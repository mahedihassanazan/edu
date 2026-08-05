# 🧬 DESeq2 Manual Tutorial (Glaucoma RNA-seq)

**Date:** 31 July 2026
**Tissue Analyzed:** Retina

## ১. আমাদের ডেটা (The Input)
আমরা `GSE241782` নামের একটি বিশাল RNA-seq ডেটাসেট থেকে শুধু **Retina** টিস্যুর ১২টি ইঁদুরের ডেটা আলাদা করেছি। 

আমাদের কাছে ২টি ফাইল ছিল:
1. **`Retina_Counts.csv`**: এখানে ১২টি ইঁদুরের ৫৫,০০০ জিনের Raw Count (রিসিট) দেওয়া ছিল।
2. **`Retina_Metadata.csv`**: এখানে লেখা ছিল কোন ইঁদুরটি সুস্থ (Control) আর কোনটিতে গ্লুকোমা (Glaucoma) আছে।

## ২. R-এ ডেটা লোড করা (Loading Data)
টার্মিনালে `R` লিখে কনসোল ওপেন করার পর, আমরা নিচের কমান্ডগুলো দিয়ে ফাইল দুটোকে R-এর মেমোরিতে লোড করেছি:

```R
# ফোল্ডারের লোকেশন সেট করা
setwd("/home/hmalik342/Working/Python-BioPython/Demo Project/Module 2")

# ফাইল দুটো ইমপোর্ট করা (প্রথম কলামকে জিনের নাম হিসেবে ধরা হয়েছে)
counts <- read.csv("Retina_Counts.csv", row.names=1)
colData <- read.csv("Retina_Metadata.csv", row.names=1)
```

## ৩. DESeq2 প্যাকেজ লোড এবং ডেটা ঠিক করা
R মাঝেমাঝে ফাইলের নামের মাঝখানের হাইফেন (`-`) কে ডট (`.`) বানিয়ে দেয়। তাই আমরা `colnames` কমান্ড দিয়ে সেই এররটি ফিক্স করেছি।

```R
# DESeq2 লাইব্রেরি কল করা
library(DESeq2)

# নামের হাইফেন/ডট সমস্যাটি ফিক্স করা
colnames(counts) <- rownames(colData)
```

## ৪. DESeq2 এর আসল ম্যাজিক (The Analysis)
এই ধাপে আমরা DESeq2 কে বলেছি—"তুমি মেটাডেটার `Condition` কলাম (Control vs Glaucoma) দেখে ডেটাগুলো স্ক্যান করো।"

```R
# DESeq2 কে ডেটা বুঝিয়ে দেওয়া
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = colData, 
                              design = ~ Condition)

# আসল এনালাইসিস শুরু করা
dds <- DESeq(dds)
```

## ৫. রেজাল্ট বের করা এবং সেভ করা (The Output)
স্ক্যান শেষ হওয়ার পর, আমরা `Glaucoma` এবং `Control` এর ভেতরের পার্থক্যগুলো (DEGs) বের করে এনেছি এবং সেটিকে আমাদের কম্পিউটারে একটি সুন্দর `CSV` ফাইলে সেভ করেছি।

```R
# রেজাল্ট বের করে আনা
res <- results(dds, contrast=c("Condition", "Glaucoma", "Control"))

# রেজাল্টটি কম্পিউটারে সেভ করা
write.csv(as.data.frame(res), file="Retina_DEGs_Result.csv")
```

---
**Summary:**
এই কাজগুলোর মাধ্যমে আমরা ৫৫ হাজার জিনের বিশাল জঙ্গল থেকে শুধু সেই জিনগুলোকে ফিল্টার করে এনেছি, যেগুলোর জন্য ইঁদুরের চোখে গ্লুকোমা রোগটি হয়েছে!
