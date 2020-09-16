# -*- coding: utf-8 -*-
"""
Created on Tue Sep 08 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Explore LUSC example of MutSigCV (cloned from genepattern github khttps://github.com/genepattern/MutSigCV). Compare
input files from LUSC to reference input files from http://software.broadinstitute.org/cancer/cga/mutsig_run
(MutSigCV1.41)
"""

import matplotlib.pyplot as plt
import os
import numpy  as np
import pandas as pd

plt.rc("text", usetex=True)

#### # 1. LOAD
#### # ###################################################################################################### 

folder = "./tools/MutSigCV/gpunit/example/data/"

df_maf = pd.read_csv(os.path.join(folder, "LUSC.mutations.maf"), sep='\t', low_memory=False)
df_n = df_maf[["Hugo_Symbol", "Tumor_Sample_Barcode", "is_coding", "is_silent", "categ"]].copy()
df_N = pd.read_csv(os.path.join(folder, "LUSC.coverage.txt"), sep='\t')
df_V = pd.read_csv(os.path.join(folder, "genes.covariates.txt"), sep='\t')

#### Comments
#### Mutations are either coding (exons, splice-sites, etc) or non-coding.
#### Among coding mutations, the variant may be silent or non silent.

#### Output of crosstab Variant_Classification, is_coding on LUSC
##                            is_coding 
##                             0      1
## Variant_Classification
## 3'UTR                    3491      0
## 5'Flank                  1255      0
## 5'UTR                    1634      0
## Frame_Shift_Del             0    518
## Frame_Shift_Ins             0    118
## In_Frame_Del                0     47
## In_Frame_Ins                0      3
## Intron                  63788      0
## Missense_Mutation           0  44831
## Nonsense_Mutation           0   3896
## Nonstop_Mutation            0     58
## RNA                        33      0
## Silent                      0  16261
## Splice_Site                 0   1307
## Translation_Start_Site      0    103


#### In the LUSC examples, categories are defined 
#### (i) transitions in C’s or G’s in CpG dinucleotides;  (C>T, G>C at CpG)
#### (ii) transversions in C’s or G’s in CpG dinucleotides; (C>A, C>G, G>C, G>T at CpG)
#### (iii) transitions in other C’s or G’s; (C>T, G>C not at CpG)
#### (iv) transversions in other C’s or G’s;  (C>A, C>G, G>C, G>T not at CpG)
#### (v) transitions at A’s or T’s;  (A>G or T>C)
#### (vi) transversions in A’s or T’s; (A>C, A>T, T>A, T>G and not in vii)
#### (vii) small insertions/deletions, nonsense and splice site mutations.

### rq: for multinucleotide subsitutions, consider the first nucleotide change to assign a category
### e.g AA>TT is in (vi), AA>GT is in (v)

#### REMINDER
####        Tri
####    A <-----> G   < Purines
####      |\Trv/|
####      | \ / |
#### Trv  |  \  | Trv
####      | / \ |
####      |/Trv\|
####    C <-----> T   < Pyrimidines
####        Tri
####

folder = "./tools/MutSigCV_1.41/references/"
df_ref_N = pd.read_csv(os.path.join(folder, "exome_full192.coverage.txt"), sep='\t')
df_ref_V = pd.read_csv(os.path.join(folder, "gene.covariates.txt"), sep='\t')

#### # 2. COMPARE LUSC AND REF
#### # ###################################################################################################### 

#### group by genes and compare
cols_tumor = [x for x in df_N.columns if x.startswith("LUSC")]
df_N_g = df_N.groupby("gene")[cols_tumor].sum()
df_N_g = df_N_g.mean(axis=1)

df_ref_N_g = df_ref_N.groupby("gene")["coverage"].sum()

#### plot on common genes
common_genes = list(set(df_N_g.index).intersection(set(df_ref_N_g.index)))
X = df_N_g[common_genes]
Y = df_ref_N_g[common_genes]
corr = np.corrcoef(X, Y)[0,1]

X_m = X.to_frame("coverage")
X_m.insert(0, "intercept", 1)
X_m = X_m.values
b = np.linalg.inv((X_m.T.dot(X_m))).dot(X_m.T).dot(Y.values)
r2 = 1 - np.sum(np.square((Y.values - X_m.dot(b))))/np.sum(np.square((Y.values - Y.mean())))

X_line = [min(X), max(X)]
Y_line = [b[0] + min(X)*b[1], b[0] + max(X)*b[1]]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14,8))
ax.scatter(X, Y, marker="+", color="darkred", s=100)
ax.plot(X_line, Y_line,  color="black", ls="--", lw=3, label="y = %.4g x + %.4g" % (b[1], b[0]))
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_xlim([0,50000])
ax.set_ylim([0,50000])
ax.set_xlabel("LUSC avg coverage", fontsize=20, fontweight="medium")
ax.set_ylabel("Ref coverage", fontsize=20, fontweight="medium")
ax.tick_params(axis="both", which="both", labelsize=16, length=4)
ax.legend(loc="best", fontsize=20, frameon=False)

filepath = "./example/coverage_ref_vs_lusc.pdf"
plt.savefig(filepath, bbox_inches="tight")

