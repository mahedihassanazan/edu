"""
CRISPR-Cas9 gRNA Design Tool for ERBB2/HER2
=============================================
এই স্ক্রিপ্টটি input_file.txt থেকে ERBB2 সিকোয়েন্স পড়বে
এবং সেরা gRNA গুলো খুঁজে বের করবে।

ব্যবহার:
    python crispr_grna_design.py

আউটপুট:
    grna_results.txt - সব gRNA রেজাল্ট
    top_grna.csv    - সেরা gRNA গুলো
"""

import re
import csv
from collections import Counter

# ============================================================
# ধাপ ১: সিকোয়েন্স লোড করা
# ============================================================

def load_fasta(filepath):
    """FASTA ফাইল থেকে সিকোয়েন্স পড়া"""
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
# ধাপ ২: PAM সাইট খোঁজা এবং gRNA এক্সট্র্যাক্ট করা
# ============================================================

def find_grnas(sequence, pam="NGG", grna_length=20):
    """
    SpCas9 এর জন্য NGG PAM সাইট খুঁজে gRNA ডিজাইন করা।
    প্রতিটি NGG এর আগের 20 nt হলো gRNA।
    """
    grnas = []
    seq_len = len(sequence)

    # Forward strand: [20nt gRNA][NGG]
    pam_pattern_fwd = re.compile(r'(?=([ACGT]{' + str(grna_length) + r'}[ACGT]GG))')
    for match in pam_pattern_fwd.finditer(sequence):
        full = match.group(1)
        grna_seq = full[:grna_length]
        pam_seq = full[grna_length:]
        position = match.start()
        grnas.append({
            "sequence": grna_seq,
            "pam": pam_seq,
            "position": position,
            "strand": "+"
        })

    # Reverse complement strand: [CCN][20nt gRNA reverse complement]
    rev_comp = reverse_complement(sequence)
    for match in pam_pattern_fwd.finditer(rev_comp):
        full = match.group(1)
        grna_seq = full[:grna_length]
        pam_seq = full[grna_length:]
        # মূল সিকোয়েন্সে অবস্থান
        rc_pos = match.start()
        orig_pos = seq_len - rc_pos - grna_length - 3
        grnas.append({
            "sequence": grna_seq,
            "pam": pam_seq,
            "position": orig_pos,
            "strand": "-"
        })

    return grnas


def reverse_complement(seq):
    """সিকোয়েন্সের রিভার্স কমপ্লিমেন্ট"""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


# ============================================================
# ধাপ ৩: gRNA স্কোরিং (Rule Set II / Doench-style)
# ============================================================

def gc_content(seq):
    """GC কন্টেন্ট শতকরা হিসাব"""
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100


def score_grna(grna_seq):
    """
    gRNA কোয়ালিটি স্কোর (0-100):
    - GC content: 40-70% হলে ভালো
    - No poly-T (TTTT) stretch: RNA Pol III terminator এড়ানো
    - No homopolymer runs
    - Seed region (3' end 12nt) বিশ্লেষণ
    """
    score = 100
    penalties = []

    gc = gc_content(grna_seq)

    # GC content penalty
    if gc < 40:
        penalty = (40 - gc) * 1.5
        score -= penalty
        penalties.append(f"Low GC ({gc:.0f}%): -{penalty:.0f}")
    elif gc > 70:
        penalty = (gc - 70) * 1.5
        score -= penalty
        penalties.append(f"High GC ({gc:.0f}%): -{penalty:.0f}")

    # Poly-T stretch (RNA Pol III terminator)
    if "TTTT" in grna_seq:
        score -= 25
        penalties.append("Poly-T (TTTT): -25")

    # Homopolymer runs (5+ same base)
    for base in "ACGT":
        if base * 5 in grna_seq:
            score -= 15
            penalties.append(f"Homopolymer ({base*5}): -15")

    # Seed region GC (last 12 nt, 3' end)
    seed = grna_seq[-12:]
    seed_gc = gc_content(seed)
    if seed_gc < 30:
        score -= 10
        penalties.append(f"Low seed GC ({seed_gc:.0f}%): -10")

    # Starting G (U6 promoter efficiency)
    if grna_seq[0] == "G":
        score += 5

    # Penalize starting with multiple Ts
    if grna_seq.startswith("TT"):
        score -= 10
        penalties.append("Starts with TT: -10")

    score = max(0, min(100, score))
    return round(score, 1), gc, penalties


def simple_offtarget_estimate(grna_seq, full_sequence):
    """
    সিম্পল অফ-টার্গেট এস্টিমেশন:
    সিকোয়েন্সে গাইডের কতটা মিলে যাচ্ছে তা দেখা।
    (সম্পূর্ণ অফ-টার্গেট বিশ্লেষণের জন্য Cas-OFFinder ব্যবহার করতে হবে)
    """
    # Seed region (last 12 nt) কতবার আছে
    seed = grna_seq[-12:]
    count = full_sequence.count(seed) + reverse_complement(full_sequence).count(seed)
    # যত কম তত ভালো (1 মানে শুধু নিজের জায়গায়)
    specificity_score = max(0, 100 - (count - 1) * 20)
    return specificity_score, count


# ============================================================
# ধাপ ৪: গুরুত্বপূর্ণ ডোমেইন চিহ্নিত করা
# ============================================================

# ERBB2/HER2 এর গুরুত্বপূর্ণ ডোমেইন (আনুমানিক CDS position)
ERBB2_DOMAINS = [
    {"name": "Signal Peptide",       "start": 1,    "end": 75},
    {"name": "Extracellular Domain I","start": 76,   "end": 630},
    {"name": "Transmembrane Domain", "start": 2533,  "end": 2598},
    {"name": "Tyrosine Kinase Domain","start": 2654,  "end": 3759},
    {"name": "C-terminal Domain",    "start": 3760,  "end": 3888},
]

def identify_domain(position):
    """gRNA এর অবস্থান কোন ডোমেইনে"""
    for domain in ERBB2_DOMAINS:
        if domain["start"] <= position <= domain["end"]:
            return domain["name"]
    return "Intergenic/UTR"


# ============================================================
# ধাপ ৫: মেইন ফাংশন
# ============================================================

def main():
    print("=" * 60)
    print("  CRISPR-Cas9 gRNA Design Tool for ERBB2/HER2")
    print("=" * 60)

    # সিকোয়েন্স লোড
    try:
        header, sequence = load_fasta("input_file.txt")
        print(f"\n✅ সিকোয়েন্স লোড সফল!")
        print(f"   Gene: {header[:60]}")
        print(f"   Length: {len(sequence)} bp")
    except FileNotFoundError:
        print("❌ input_file.txt পাওয়া যায়নি! ফাইলটি একই ফোল্ডারে রাখুন।")
        return

    # gRNA খোঁজা
    print(f"\n🔍 NGG PAM সাইট খোঁজা হচ্ছে...")
    grnas = find_grnas(sequence)
    print(f"   মোট gRNA পাওয়া গেছে: {len(grnas)}")

    # স্কোরিং
    print(f"\n📊 gRNA স্কোরিং চলছে...")
    scored_grnas = []
    for g in grnas:
        on_score, gc, penalties = score_grna(g["sequence"])
        spec_score, seed_count = simple_offtarget_estimate(g["sequence"], sequence)
        domain = identify_domain(g["position"])

        # Combined score
        combined = (on_score * 0.6) + (spec_score * 0.4)

        scored_grnas.append({
            **g,
            "on_target_score": on_score,
            "gc_content": round(gc, 1),
            "specificity_score": spec_score,
            "combined_score": round(combined, 1),
            "domain": domain,
            "penalties": "; ".join(penalties) if penalties else "None",
        })

    # সর্বোচ্চ স্কোরে সাজানো
    scored_grnas.sort(key=lambda x: x["combined_score"], reverse=True)

    # সেরা ১০টি দেখানো
    print(f"\n{'='*60}")
    print(f"  🏆 সেরা gRNA তালিকা (Top 10)")
    print(f"{'='*60}")
    print(f"{'Rank':<5} {'gRNA Sequence (20nt)':<22} {'PAM':<5} {'Strand':<7} "
          f"{'On-Score':<10} {'GC%':<6} {'Spec.':<7} {'Combined':<9} {'Domain'}")
    print("-" * 110)

    top_grnas = scored_grnas[:10]
    for i, g in enumerate(top_grnas, 1):
        print(f"{i:<5} {g['sequence']:<22} {g['pam']:<5} {g['strand']:<7} "
              f"{g['on_target_score']:<10} {g['gc_content']:<6} "
              f"{g['specificity_score']:<7} {g['combined_score']:<9} {g['domain']}")

    # সেরা ৩টি বিস্তারিত
    print(f"\n{'='*60}")
    print(f"  🎯 সেরা ৩টি gRNA - বিস্তারিত বিশ্লেষণ")
    print(f"{'='*60}")

    for i, g in enumerate(top_grnas[:3], 1):
        print(f"""
┌─ gRNA #{i} {'─'*45}
│  Sequence:        5'-{g['sequence']}-3'
│  PAM:             {g['pam']} (NGG)
│  Position:        {g['position']} bp
│  Strand:          {g['strand']}
│  GC Content:      {g['gc_content']}%
│  On-Target Score: {g['on_target_score']}/100
│  Specificity:     {g['specificity_score']}/100
│  Combined Score:  {g['combined_score']}/100
│  Domain:          {g['domain']}
│  Penalties:       {g['penalties']}
│  
│  ➡️  Cas-OFFinder এ যাচাই করুন:
│     http://www.rgenome.net/cas-offinder/
│     Sequence: {g['sequence']}NGG
└{'─'*50}""")

    # CSV সেভ করা
    print(f"\n💾 ফাইল সেভ হচ্ছে...")
    fields = ["rank", "sequence", "pam", "position", "strand",
              "on_target_score", "gc_content", "specificity_score",
              "combined_score", "domain", "penalties"]

    with open("grna_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for i, g in enumerate(scored_grnas[:50], 1):
            row = {k: g[k] for k in fields if k != "rank"}
            row["rank"] = i
            writer.writerow(row)

    print(f"   ✅ grna_results.csv সেভ হয়েছে (Top 50 gRNA)")

    # Summary
    print(f"""
{'='*60}
  📋 সারসংক্ষেপ
{'='*60}
  মোট gRNA পাওয়া গেছে:    {len(grnas)}
  Score > 70 (ভালো):       {sum(1 for g in scored_grnas if g['combined_score'] > 70)}
  Score > 80 (খুব ভালো):   {sum(1 for g in scored_grnas if g['combined_score'] > 80)}

  পরবর্তী ধাপ:
  1. উপরের সেরা ৩টি gRNA নিন
  2. Cas-OFFinder এ অফ-টার্গেট চেক করুন
  3. CRISPOR এ যাচাই করুন (যদি কাজ করে)
  
  Cas-OFFinder লিংক: http://www.rgenome.net/cas-offinder/
{'='*60}
""")


if __name__ == "__main__":
    main()
