"""
CRISPR Off-Target Analysis Tool
=================================
এই স্ক্রিপ্টটি আপনার gRNA গুলোর off-target site
ERBB2 সিকোয়েন্সের মধ্যে খুঁজে বের করবে।

ব্যবহার:
    python offtarget_analysis.py

আউটপুট:
    offtarget_report.txt  - বিস্তারিত রিপোর্ট
    offtarget_results.csv - CSV ফরম্যাট
"""

import csv
import itertools
from datetime import datetime

# ============================================================
# ধাপ ১: আপনার সেরা gRNA গুলো এখানে দিন
# (আগের স্ক্রিপ্ট থেকে পাওয়া Top 3 gRNA)
# ============================================================

GRNAS = [
    {"name": "gRNA_1", "sequence": "CCATTGGGACCGGAGAAACC"},
    {"name": "gRNA_2", "sequence": "CGGAGCCGCAGTGAGCACCA"},
    {"name": "gRNA_3", "sequence": "GCGAGCACCCAAGTGTGCAC"},
]

MAX_MISMATCHES = 3   # কতটা mismatch পর্যন্ত দেখবেন
PAM = "NGG"          # SpCas9 PAM

# ============================================================
# গুরুত্বপূর্ণ জিন — এদের মধ্যে off-target হলে বিপদ!
# ============================================================
CRITICAL_GENES = [
    "TP53", "BRCA1", "BRCA2", "RB1", "APC",
    "PTEN", "VHL", "MLH1", "MSH2", "CDKN2A"
]

# ============================================================
# Helper Functions
# ============================================================

def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def count_mismatches(seq1, seq2):
    """দুটো সিকোয়েন্সের মধ্যে mismatch গণনা"""
    if len(seq1) != len(seq2):
        return 999
    return sum(a != b for a, b in zip(seq1, seq2))


def pam_matches(pam_seq):
    """NGG PAM চেক করা (N = যেকোনো base)"""
    if len(pam_seq) < 3:
        return False
    return pam_seq[1] == "G" and pam_seq[2] == "G"


def load_fasta(filepath):
    sequence = ""
    header = ""
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
            else:
                sequence += line.upper()
    return header, sequence


# ============================================================
# মূল Off-Target Search ফাংশন
# ============================================================

def find_offtargets(grna_seq, reference_seq, max_mismatches=3):
    """
    রেফারেন্স সিকোয়েন্সে gRNA এর off-target সাইট খোঁজা।
    উভয় strand চেক করা হবে।
    """
    results = []
    grna_len = len(grna_seq)
    ref_len = len(reference_seq)

    # Forward এবং reverse strand উভয়ই চেক
    strands = [
        ("+", reference_seq),
        ("-", reverse_complement(reference_seq))
    ]

    for strand, seq in strands:
        for i in range(len(seq) - grna_len - 3):
            candidate = seq[i:i + grna_len]
            pam_seq = seq[i + grna_len:i + grna_len + 3]

            # PAM চেক
            if not pam_matches(pam_seq):
                continue

            # Mismatch গণনা
            mismatches = count_mismatches(grna_seq, candidate)

            if mismatches <= max_mismatches:
                # Seed region mismatch (শেষ 12 nt — বেশি গুরুত্বপূর্ণ)
                seed_mm = count_mismatches(grna_seq[-12:], candidate[-12:])

                # Original position
                if strand == "+":
                    orig_pos = i
                else:
                    orig_pos = ref_len - i - grna_len

                results.append({
                    "position": orig_pos,
                    "strand": strand,
                    "target_seq": candidate,
                    "pam": pam_seq,
                    "total_mismatches": mismatches,
                    "seed_mismatches": seed_mm,
                    "is_perfect": mismatches == 0,
                    "risk": classify_risk(mismatches, seed_mm),
                })

    # Mismatch অনুযায়ী সাজানো
    results.sort(key=lambda x: (x["total_mismatches"], x["seed_mismatches"]))
    return results


def classify_risk(total_mm, seed_mm):
    """Off-target এর ঝুঁকি মূল্যায়ন"""
    if total_mm == 0:
        return "ON-TARGET ✅"
    elif total_mm == 1 and seed_mm == 0:
        return "HIGH RISK ❌"
    elif total_mm == 1 and seed_mm == 1:
        return "MEDIUM RISK ⚠️"
    elif total_mm == 2 and seed_mm == 0:
        return "MEDIUM RISK ⚠️"
    elif total_mm == 2:
        return "LOW RISK 🟡"
    else:
        return "VERY LOW RISK 🟢"


def calculate_specificity_score(offtargets):
    """
    Specificity score হিসাব (0-100)
    MIT Specificity Score এর simplified version
    """
    if not offtargets:
        return 100.0

    penalty = 0
    for ot in offtargets:
        mm = ot["total_mismatches"]
        seed_mm = ot["seed_mismatches"]

        if mm == 0:
            continue  # on-target, skip
        elif mm == 1:
            penalty += 30 if seed_mm == 0 else 15
        elif mm == 2:
            penalty += 10 if seed_mm == 0 else 5
        elif mm == 3:
            penalty += 2

    score = max(0, 100 - penalty)
    return round(score, 1)


# ============================================================
# রিপোর্ট তৈরি
# ============================================================

def print_separator(char="=", width=65):
    print(char * width)


def generate_report(grna_name, grna_seq, offtargets, spec_score):
    """একটি gRNA এর বিস্তারিত রিপোর্ট"""
    on_target = [o for o in offtargets if o["total_mismatches"] == 0]
    off_1mm   = [o for o in offtargets if o["total_mismatches"] == 1]
    off_2mm   = [o for o in offtargets if o["total_mismatches"] == 2]
    off_3mm   = [o for o in offtargets if o["total_mismatches"] == 3]

    print(f"\n{'─'*65}")
    print(f"  🔬 {grna_name}: 5'-{grna_seq}-3'")
    print(f"{'─'*65}")
    print(f"  Specificity Score:  {spec_score}/100", end="  ")
    if spec_score >= 90:
        print("→ ✅ EXCELLENT — ব্যবহারযোগ্য")
    elif spec_score >= 70:
        print("→ ⚠️  GOOD — সতর্কতার সাথে ব্যবহার করুন")
    else:
        print("→ ❌ POOR — এই gRNA বাদ দিন")

    print(f"\n  Hit Summary:")
    print(f"    On-target  (0 mismatch): {len(on_target)} টি")
    print(f"    1 mismatch off-target:   {len(off_1mm)} টি  {'❌ বিপদ!' if off_1mm else '✅'}")
    print(f"    2 mismatch off-target:   {len(off_2mm)} টি  {'⚠️ সতর্ক' if off_2mm else '✅'}")
    print(f"    3 mismatch off-target:   {len(off_3mm)} টি")

    # বিস্তারিত off-target তালিকা
    dangerous = [o for o in offtargets if o["total_mismatches"] <= 2]
    if dangerous:
        print(f"\n  ⚠️  বিপজ্জনক Off-Target Sites:")
        print(f"  {'Position':<10} {'Strand':<8} {'Target Seq':<22} {'PAM':<5} {'MM':<4} {'Seed MM':<9} {'Risk'}")
        print(f"  {'─'*75}")
        for ot in dangerous[:10]:
            print(f"  {ot['position']:<10} {ot['strand']:<8} {ot['target_seq']:<22} "
                  f"{ot['pam']:<5} {ot['total_mismatches']:<4} {ot['seed_mismatches']:<9} {ot['risk']}")
    else:
        print(f"\n  ✅ কোনো বিপজ্জনক off-target site পাওয়া যায়নি!")


# ============================================================
# মেইন ফাংশন
# ============================================================

def main():
    print_separator()
    print("  CRISPR Off-Target Analysis Tool")
    print(f"  তারিখ: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print_separator()

    # সিকোয়েন্স লোড
    try:
        header, sequence = load_fasta("input_file.txt")
        print(f"\n✅ সিকোয়েন্স লোড সফল!")
        print(f"   {header[:60]}")
        print(f"   দৈর্ঘ্য: {len(sequence)} bp")
    except FileNotFoundError:
        print("❌ input_file.txt পাওয়া যায়নি!")
        return

    print(f"\n   gRNA তালিকা ({len(GRNAS)}টি):")
    for g in GRNAS:
        print(f"   • {g['name']}: {g['sequence']}")

    print(f"\n🔍 Off-target বিশ্লেষণ শুরু হচ্ছে...")
    print(f"   Max mismatches: {MAX_MISMATCHES}")
    print(f"   PAM: {PAM}\n")

    # প্রতিটি gRNA বিশ্লেষণ
    all_results = []
    summary_data = []

    for grna in GRNAS:
        name = grna["name"]
        seq  = grna["sequence"]
        print(f"  ⏳ {name} বিশ্লেষণ করা হচ্ছে...", end=" ", flush=True)

        offtargets = find_offtargets(seq, sequence, MAX_MISMATCHES)
        spec_score = calculate_specificity_score(offtargets)

        print(f"সম্পন্ন! ({len(offtargets)} hit পাওয়া গেছে)")

        generate_report(name, seq, offtargets, spec_score)

        # CSV এর জন্য ডেটা
        for ot in offtargets:
            all_results.append({
                "grna_name": name,
                "grna_sequence": seq,
                "hit_position": ot["position"],
                "strand": ot["strand"],
                "target_sequence": ot["target_seq"],
                "pam": ot["pam"],
                "total_mismatches": ot["total_mismatches"],
                "seed_mismatches": ot["seed_mismatches"],
                "risk": ot["risk"],
                "specificity_score": spec_score,
            })

        summary_data.append({
            "name": name,
            "sequence": seq,
            "spec_score": spec_score,
            "total_hits": len(offtargets),
            "mm0": sum(1 for o in offtargets if o["total_mismatches"] == 0),
            "mm1": sum(1 for o in offtargets if o["total_mismatches"] == 1),
            "mm2": sum(1 for o in offtargets if o["total_mismatches"] == 2),
            "mm3": sum(1 for o in offtargets if o["total_mismatches"] == 3),
        })

    # ফাইনাল সারসংক্ষেপ
    print(f"\n{'='*65}")
    print(f"  📊 ফাইনাল তুলনামূলক সারসংক্ষেপ")
    print(f"{'='*65}")
    print(f"  {'gRNA':<10} {'Sequence':<22} {'Spec.':<8} {'MM0':<5} {'MM1':<5} {'MM2':<5} {'MM3':<5} {'সিদ্ধান্ত'}")
    print(f"  {'─'*75}")

    best_grna = None
    best_score = -1

    for s in summary_data:
        if s["spec_score"] >= 90:
            verdict = "✅ ব্যবহার করুন"
        elif s["spec_score"] >= 70:
            verdict = "⚠️  সতর্কতা"
        else:
            verdict = "❌ বাদ দিন"

        print(f"  {s['name']:<10} {s['sequence']:<22} {s['spec_score']:<8} "
              f"{s['mm0']:<5} {s['mm1']:<5} {s['mm2']:<5} {s['mm3']:<5} {verdict}")

        if s["spec_score"] > best_score and s["mm1"] == 0:
            best_score = s["spec_score"]
            best_grna = s

    # সেরা gRNA সুপারিশ
    print(f"\n{'─'*65}")
    if best_grna:
        print(f"  🏆 সুপারিশকৃত gRNA: {best_grna['name']}")
        print(f"     Sequence: 5'-{best_grna['sequence']}-3'")
        print(f"     Specificity: {best_grna['spec_score']}/100")
        print(f"     এই gRNA টি সবচেয়ে নিরাপদ এবং কার্যকর।")
    else:
        print(f"  ⚠️  সব gRNA এ 1-mismatch off-target আছে।")
        print(f"     Cas-OFFinder দিয়ে পুরো human genome এ যাচাই করুন।")

    # CSV সেভ
    print(f"\n💾 ফাইল সেভ হচ্ছে...")
    if all_results:
        fields = ["grna_name", "grna_sequence", "hit_position", "strand",
                  "target_sequence", "pam", "total_mismatches",
                  "seed_mismatches", "risk", "specificity_score"]
        with open("offtarget_results.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(all_results)
        print(f"   ✅ offtarget_results.csv সেভ হয়েছে ({len(all_results)} rows)")

    print(f"""
{'='*65}
  ✅ Off-Target Analysis সম্পন্ন!

  পরবর্তী ধাপ:
  → Molecular Docking (AutoDock Vina দিয়ে)
    HER2 প্রোটিনে Lapatinib ড্রাগ কতটা ভালো বাইন্ড করে
    সেটা Python দিয়ে বিশ্লেষণ করব।
{'='*65}
""")


if __name__ == "__main__":
    main()
