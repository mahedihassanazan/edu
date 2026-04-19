"""
========================================================
  CRISPR-Cas9 + gRNA Molecular Docking Pipeline
  AutoDock Vina Integration
  Simplified Ligand Generation
========================================================
"""

import os
import sys
import subprocess
import hashlib
import logging
import requests
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ============================================================
#  LOGGING
# ============================================================

logging.basicConfig(
    level=logging.INFO,
    format="  %(levelname)s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


# ============================================================
#  CONFIGURATION (typed dataclasses — no more raw dicts)
# ============================================================

@dataclass
class GRNASequence:
    id: str
    sequence: str

    def gc_content(self) -> float:
        """Return GC content as a fraction (0.0–1.0)."""
        return sum(1 for n in self.sequence if n in "GC") / len(self.sequence)

    def stable_seed(self) -> int:
        """Deterministic seed from sequence (safe across Python versions)."""
        return int(hashlib.md5(self.sequence.encode()).hexdigest(), 16) % (2 ** 32)


@dataclass
class DockingBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


@dataclass
class PipelineConfig:
    grna_sequences: list[GRNASequence] = field(default_factory=list)
    cas9_pdb_id: str = "4ZT0"
    vina_path: str = "vina"
    docking_box: DockingBox = field(
        default_factory=lambda: DockingBox(
            center_x=10.0, center_y=-5.0, center_z=20.0,
            size_x=25.0,   size_y=25.0,   size_z=25.0,
        )
    )
    exhaustiveness: int = 8
    num_modes: int = 9
    work_dir: Path = Path("docking_workspace")
    output_csv: str = "docking_summary.csv"


# ---- Default sequences --------------------------------------------------

DEFAULT_CONFIG = PipelineConfig(
    grna_sequences=[
        GRNASequence(id="gRNA_1", sequence="CCATTGGGACCGGAGAAACC"),
        GRNASequence(id="gRNA_2", sequence="CGGAGCCGCAGTGAGCACCA"),
        GRNASequence(id="gRNA_3", sequence="GCGAGCACCCAAGTGTGCAC"),
    ]
)


# ============================================================
#  STEP 1: Download Cas9 Structure
# ============================================================

def download_cas9_pdb(cfg: PipelineConfig) -> Path:
    pdb_path = cfg.work_dir / f"{cfg.cas9_pdb_id}.pdb"
    if pdb_path.exists():
        log.info(f"[OK] Cas9 PDB already present: {pdb_path}")
        return pdb_path

    log.info(f"[>>] Downloading Cas9 structure ({cfg.cas9_pdb_id})...")
    url = f"https://files.rcsb.org/download/{cfg.cas9_pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        pdb_path.write_text(response.text)
        log.info(f"[OK] Downloaded: {pdb_path}")
        return pdb_path
    except requests.RequestException as e:
        log.warning(f"[WARN] Download failed ({e}). Using mock Cas9 PDB.")
        return _create_mock_cas9_pdb(pdb_path)


def _create_mock_cas9_pdb(path: Path) -> Path:
    """Minimal valid PDB so the pipeline can continue in offline/CI environments."""
    path.write_text(
        "REMARK Mock Cas9 — download failed\n"
        "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N\n"
        "END\n"
    )
    log.warning("[WARN] Mock PDB written — docking scores will be simulated.")
    return path


# ============================================================
#  STEP 2: Prepare Receptor (PDB ? PDBQT)
# ============================================================

def prepare_receptor_pdbqt(pdb_path: Path, cfg: PipelineConfig) -> Path:
    pdbqt_path = cfg.work_dir / "receptor.pdbqt"
    log.info("[STEP 2] Preparing Cas9 receptor...")

    # Attempt OpenBabel conversion
    try:
        result = subprocess.run(
            ["obabel", str(pdb_path), "-O", str(pdbqt_path), "-xr"],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode == 0 and pdbqt_path.exists():
            log.info(f"[OK] Receptor prepared via OpenBabel: {pdbqt_path}")
            return pdbqt_path
        log.warning(f"[WARN] obabel returned non-zero ({result.returncode}). Falling back.")
    except FileNotFoundError:
        log.warning("[WARN] obabel not found — using manual PDBQT conversion.")
    except subprocess.TimeoutExpired:
        log.warning("[WARN] obabel timed out — using manual PDBQT conversion.")

    # Manual fallback: write a simplified but column-correct PDBQT
    log.warning("[WARN] Manual PDBQT conversion — results may be less accurate.")
    out_lines: list[str] = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and len(line) >= 54:
                # PDBQT format: cols 1-66 from PDB, then partial charge + atom type
                padded = f"{line[:66]:<66}  +0.000 C\n"
                out_lines.append(padded)
            elif line.startswith(("TER", "END")):
                out_lines.append(line)

    pdbqt_path.write_text("".join(out_lines))
    log.info("[OK] Manual receptor PDBQT written.")
    return pdbqt_path


# ============================================================
#  STEP 3: Prepare Simplified gRNA Ligand
# ============================================================

def prepare_ligand_pdbqt(grna: GRNASequence, cfg: PipelineConfig) -> Path:
    """
    Build a minimal single-chain PDBQT ligand from the first 12 nucleotides.
    Each nucleotide is represented as a single Carbon atom for Vina compatibility.
    """
    ligand_path = cfg.work_dir / f"{grna.id}_ligand.pdbqt"
    box = cfg.docking_box
    lines: list[str] = [
        f"REMARK Ligand: {grna.id}\n",
        f"REMARK Sequence: {grna.sequence}\n",
    ]

    for i, _ in enumerate(grna.sequence[:12]):
        x = box.center_x + (i * 3.0)
        y = box.center_y
        z = box.center_z
        # PDBQT ATOM record — fixed-width columns matching Vina's parser
        lines.append(
            f"ATOM  {i+1:5d}  C   UNK  L {i+1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     +0.000 C\n"
        )

    lines += ["TER\n", "END\n"]
    ligand_path.write_text("".join(lines))
    log.info(f"  [OK] Ligand PDBQT written: {ligand_path.name}")
    return ligand_path


# ============================================================
#  STEP 4: Vina Configuration File
# ============================================================

def _write_vina_config(
    receptor_path: Path,
    ligand_path: Path,
    output_path: Path,
    cfg: PipelineConfig,
) -> Path:
    config_path = cfg.work_dir / "vina_config.txt"
    box = cfg.docking_box
    config_path.write_text(
        f"receptor       = {receptor_path}\n"
        f"ligand         = {ligand_path}\n"
        f"out            = {output_path}\n"
        f"center_x       = {box.center_x}\n"
        f"center_y       = {box.center_y}\n"
        f"center_z       = {box.center_z}\n"
        f"size_x         = {box.size_x}\n"
        f"size_y         = {box.size_y}\n"
        f"size_z         = {box.size_z}\n"
        f"exhaustiveness = {cfg.exhaustiveness}\n"
        f"num_modes      = {cfg.num_modes}\n"
        f"energy_range   = 4\n"
    )
    return config_path


# ============================================================
#  STEP 5: Run Vina + Scoring
# ============================================================

@dataclass
class DockingResult:
    grna_id: str
    sequence: str
    status: str                        # "success" | "simulated" | "failed"
    best_affinity: Optional[float]
    all_poses: list[float] = field(default_factory=list)
    error: Optional[str] = None


def run_vina_docking(
    grna: GRNASequence,
    receptor_path: Path,
    cfg: PipelineConfig,
) -> DockingResult:
    ligand_path = prepare_ligand_pdbqt(grna, cfg)
    output_path = cfg.work_dir / f"{grna.id}_out.pdbqt"
    config_path = _write_vina_config(receptor_path, ligand_path, output_path, cfg)

    log.info(f"  [>>] Running Vina for {grna.id}...")

    try:
        proc = subprocess.run(
            [cfg.vina_path, "--config", str(config_path)],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if proc.returncode != 0:
            log.warning(
                f"  [WARN] Vina exited with code {proc.returncode}: "
                f"{proc.stderr[:120].strip()}"
            )
            log.info("  [INFO] Falling back to simulated scoring.")
            return _simulate_docking_scores(grna)

        affinities = _parse_vina_output(proc.stdout)
        if not affinities:
            log.warning("  [WARN] Could not parse Vina output. Falling back.")
            return _simulate_docking_scores(grna)

        log.info(f"  [OK] Best affinity: {affinities[0]} kcal/mol")
        return DockingResult(
            grna_id=grna.id,
            sequence=grna.sequence,
            status="success",
            best_affinity=affinities[0],
            all_poses=affinities,
        )

    except FileNotFoundError:
        log.warning(f"  [WARN] Vina binary not found at '{cfg.vina_path}'. Simulating.")
        return _simulate_docking_scores(grna)
    except subprocess.TimeoutExpired:
        log.warning("  [WARN] Vina timed out. Simulating.")
        return _simulate_docking_scores(grna)
    except Exception as e:
        log.error(f"  [FAIL] Unexpected error: {e}")
        return DockingResult(
            grna_id=grna.id,
            sequence=grna.sequence,
            status="failed",
            best_affinity=None,
            error=str(e),
        )


def _parse_vina_output(text: str) -> list[float]:
    """Extract binding affinities from Vina stdout table."""
    affinities: list[float] = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                affinities.append(float(parts[1]))
            except ValueError:
                continue
    return affinities


def _simulate_docking_scores(grna: GRNASequence) -> DockingResult:
    """
    Fallback scoring based on GC content.
    Uses a stable MD5-derived seed so results are reproducible across Python versions.
    """
    rng = np.random.default_rng(grna.stable_seed())
    gc = grna.gc_content()
    base_affinity = -6.0 - (gc * 2.5)
    poses = sorted(
        round(base_affinity + rng.uniform(-0.5, 0.5) + i * 0.3, 1)
        for i in range(3)
    )
    log.info(f"  [SIM] Simulated best affinity: {poses[0]} kcal/mol")
    return DockingResult(
        grna_id=grna.id,
        sequence=grna.sequence,
        status="simulated",
        best_affinity=poses[0],
        all_poses=poses,
    )


# ============================================================
#  STEP 6: Analyze Results
# ============================================================

def _recommendation(affinity: float) -> str:
    if affinity <= -8.0:
        return "Recommended"
    if affinity <= -6.0:
        return "Acceptable"
    return "Low affinity"


def analyze_results(results: list[DockingResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        grna = GRNASequence(id=r.grna_id, sequence=r.sequence)
        aff = r.best_affinity if r.best_affinity is not None else float("nan")
        rows.append(
            {
                "gRNA ID": r.grna_id,
                "Sequence": r.sequence,
                "GC (%)": round(grna.gc_content() * 100, 1),
                "Affinity (kcal/mol)": aff,
                "Status": r.status,
                "Recommendation": _recommendation(aff) if not np.isnan(aff) else "N/A",
            }
        )

    df = pd.DataFrame(rows).sort_values("Affinity (kcal/mol)", na_position="last")
    return df.reset_index(drop=True)


# ============================================================
#  MAIN
# ============================================================

def run_pipeline(cfg: PipelineConfig = DEFAULT_CONFIG) -> pd.DataFrame:
    cfg.work_dir.mkdir(parents=True, exist_ok=True)
    output_path = cfg.work_dir / cfg.output_csv

    print("\n" + "=" * 60)
    print("  CRISPR DOCKING PIPELINE")
    print("=" * 60)

    # Step 1 + 2: fetch structure and prepare receptor
    pdb_path = download_cas9_pdb(cfg)
    receptor_pdbqt = prepare_receptor_pdbqt(pdb_path, cfg)

    # Steps 3–5: dock each gRNA
    print("\n[STEPS 3–5] Docking gRNA sequences...")
    results = [run_vina_docking(g, receptor_pdbqt, cfg) for g in cfg.grna_sequences]

    # Step 6: analyze
    df = analyze_results(results)
    df.to_csv(output_path, index=False)

    print("\n" + "=" * 60)
    print("  RESULTS")
    print("=" * 60)
    print(df.to_string(index=False))
    print(f"\n[OK] Saved: {output_path}")
    print("=" * 60 + "\n")

    return df


if __name__ == "__main__":
    run_pipeline()