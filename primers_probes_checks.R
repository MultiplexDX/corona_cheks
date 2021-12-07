setwd("/home/veve/Documents/sars2_byMonth/B.1.1.529/4thDec2021_1pm/")
library("Biostrings")
myMSA <- readDNAStringSet("Wuhan_B117cons80_Omicron450seqs.fa")
mySet <- DNAStringSet(myMSA)
Nseqs <- length(mySet)

# find positions of primers and probes
E_Sarbeco_P1 <- "ACACTAGCCATCCTTACTGCGCTTCG"
E_Sarbeco_R2 <- "TGTGTGCGTACTGCTGCAATAT"
E_Sarbeco_F7 <- "ATGTACTCATTCGTTTCGGAAGA"
E_Sarbeco_F3 <- "GTACTCATTCGTTTCGGAAGAGACAG"
E_Sarbeco_P2 <- "TAGCGTACTTCTTTTTCTTGCTTTCGTGGT"

RdRp_SARSr_P2 <- "CAGGTGGAACCTCATCAGGAGATGC"
RdRp_SARSr_R2 <- "GTTTTTAACATTTGTCAAGCTGTCACG"
RdRp_SARSr_F2 <- "GTGAAATGGTCATGTGTGGCGG"
RdRp_SARSr_P8 <- "TCAGGAGATGCCACAACTGCTTATGC"

Sgene_F1 <- "TCTTTTCCAATGTTACTTGGTTC"
Sgene_P3 <- "AGAGGTTTGATAACCCTGTCCTACCA"
Sgene_P4 <- "TTTGCTTCCACTGAGAAGTCTAACAT"
Sgene_R1_RC <- "AGATTCGAAGACCCAGTCCCTACT"

B117_F3 <- "GTTACTTGGTTCCATGCTATCTC"
B117_R36_RC <- "GCTTACCACAAAAACAACAAAAGTTG"

# create tables of start and end positions
# adjust the maximal number of mismatches
E_P1 <- as.data.frame(vmatchPattern(E_Sarbeco_P1, mySet, max.mismatch = 5, min.mismatch = 0))
E_R2 <- as.data.frame(vmatchPattern(E_Sarbeco_R2, mySet, max.mismatch = 5, min.mismatch = 0))
E_F7 <- as.data.frame(vmatchPattern(E_Sarbeco_F7, mySet, max.mismatch = 5, min.mismatch = 0))
E_F3 <- as.data.frame(vmatchPattern(E_Sarbeco_F3, mySet, max.mismatch = 5, min.mismatch = 0))
E_P2 <- as.data.frame(vmatchPattern(E_Sarbeco_P2, mySet, max.mismatch = 5, min.mismatch = 0))
  
R_P2 <- as.data.frame(vmatchPattern(RdRp_SARSr_P2, mySet, max.mismatch = 5, min.mismatch = 0))
R_R2 <- as.data.frame(vmatchPattern(RdRp_SARSr_R2, mySet, max.mismatch = 5, min.mismatch = 0))
R_F2 <- as.data.frame(vmatchPattern(RdRp_SARSr_F2, mySet, max.mismatch = 5, min.mismatch = 0))
R_P8 <- as.data.frame(vmatchPattern(RdRp_SARSr_P8, mySet, max.mismatch = 5, min.mismatch = 0))
  
S_F1 <- as.data.frame(vmatchPattern(Sgene_F1, mySet, max.mismatch = 5, min.mismatch = 0))
S_P3 <- as.data.frame(vmatchPattern(Sgene_P3, mySet, max.mismatch = 5, min.mismatch = 0))
S_P4 <- as.data.frame(vmatchPattern(Sgene_P4, mySet, max.mismatch = 5, min.mismatch = 0))
S_R1 <- as.data.frame(vmatchPattern(Sgene_R1_RC, mySet, max.mismatch = 5, min.mismatch = 0))
  
Br_F3 <- as.data.frame(vmatchPattern(B117_F3, mySet, max.mismatch = 5, min.mismatch = 0))
Br_R36 <- as.data.frame(vmatchPattern(B117_R36_RC, mySet, max.mismatch = 5, min.mismatch = 0))

# gene E P1
myLen <- dim(E_P1)[1]
E_P1_seqs <- c()
for (i in c(1:myLen)){
  E_P1_seqs <- c(E_P1_seqs, as.character(subseq(mySet[E_P1$group[i]], start = E_P1$start[i], end = E_P1$end[i])))
}

print(paste("Pattern", E_Sarbeco_P1, "(E_Sarbeco_P1) is NOT found in", length(mySet)-length(E_P1_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", E_Sarbeco_P1, "(E_Sarbeco_P1) is found in these forms (with max 5 mismatches):"))
table(sort(E_P1_seqs))

# gene E R2
myLen <- dim(E_R2)[1]
E_R2_seqs <- c()
for (i in c(1:myLen)){
  E_R2_seqs <- c(E_R2_seqs, as.character(subseq(mySet[E_R2$group[i]], start = E_R2$start[i], end = E_R2$end[i])))
}

print(paste("Pattern", E_Sarbeco_R2, "(E_Sarbeco_R2) is NOT found in", length(mySet)-length(E_R2_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", E_Sarbeco_R2, "(E_Sarbeco_R2) is found in these forms (with max 5 mismatches):"))
table(sort(E_R2_seqs))

# gene E F7
myLen <- dim(E_F7)[1]
E_F7_seqs <- c()
for (i in c(1:myLen)){
  E_F7_seqs <- c(E_F7_seqs, as.character(subseq(mySet[E_F7$group[i]], start = E_F7$start[i], end = E_F7$end[i])))
}

print(paste("Pattern", E_Sarbeco_F7, "(E_Sarbeco_F7) is NOT found in", length(mySet)-length(E_F7_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", E_Sarbeco_F7, "(E_Sarbeco_F7) is found in these forms (with max 5 mismatches):"))
table(sort(E_F7_seqs))

# gene E F3
myLen <- dim(E_F3)[1]
E_F3_seqs <- c()
for (i in c(1:myLen)){
  E_F3_seqs <- c(E_F3_seqs, as.character(subseq(mySet[E_F3$group[i]], start = E_F3$start[i], end = E_F3$end[i])))
}

print(paste("Pattern", E_Sarbeco_F3, "(E_Sarbeco_F3) is NOT found in", length(mySet)-length(E_F3_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", E_Sarbeco_F3, "(E_Sarbeco_F3) is found in these forms (with max 5 mismatches):"))
table(sort(E_F3_seqs))

# gene E P2
myLen <- dim(E_P2)[1]
E_P2_seqs <- c()
for (i in c(1:myLen)){
  E_P2_seqs <- c(E_P2_seqs, as.character(subseq(mySet[E_P2$group[i]], start = E_P2$start[i], end = E_P2$end[i])))
}

print(paste("Pattern", E_Sarbeco_P2, "(E_Sarbeco_P2) is NOT found in", length(mySet)-length(E_P2_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", E_Sarbeco_P2, "(E_Sarbeco_P2) is found in these forms (with max 5 mismatches):"))
table(sort(E_P2_seqs))


# gene RdRp P2
myLen <- dim(R_P2)[1]
R_P2_seqs <- c()
for (i in c(1:myLen)){
  R_P2_seqs <- c(R_P2_seqs, as.character(subseq(mySet[R_P2$group[i]], start = R_P2$start[i], end = R_P2$end[i])))
}

print(paste("Pattern", R_SARSr_P2, "(R_SARSr_P2) is NOT found in", length(mySet)-length(R_P2_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", R_SARSr_P2, "(R_SARSr_P2) is found in these forms (with max 5 mismatches):"))
table(sort(R_P2_seqs))

# gene RdRp R2
myLen <- dim(R_R2)[1]
R_R2_seqs <- c()
for (i in c(1:myLen)){
  R_R2_seqs <- c(R_R2_seqs, as.character(subseq(mySet[R_R2$group[i]], start = R_R2$start[i], end = R_R2$end[i])))
}

print(paste("Pattern", R_SARSr_R2, "(R_SARSr_R2) is NOT found in", length(mySet)-length(R_R2_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", R_SARSr_R2, "(R_SARSr_R2) is found in these forms (with max 5 mismatches):"))
table(sort(R_R2_seqs))

# gene RdRp F2
myLen <- dim(R_F2)[1]
R_F2_seqs <- c()
for (i in c(1:myLen)){
  R_F2_seqs <- c(R_F2_seqs, as.character(subseq(mySet[R_F2$group[i]], start = R_F2$start[i], end = R_F2$end[i])))
}

print(paste("Pattern", R_SARSr_F2, "(R_SARSr_F2) is NOT found in", length(mySet)-length(R_F2_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", R_SARSr_F2, "(R_SARSr_F2) is found in these forms (with max 5 mismatches):"))
table(sort(R_F2_seqs))

# gene RdRp P8
myLen <- dim(R_P8)[1]
R_P8_seqs <- c()
for (i in c(1:myLen)){
  R_P8_seqs <- c(R_P8_seqs, as.character(subseq(mySet[R_P8$group[i]], start = R_P8$start[i], end = R_P8$end[i])))
}

print(paste("Pattern", R_SARSr_P8, "(R_SARSr_P8) is NOT found in", length(mySet)-length(R_P8_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", R_SARSr_P8, "(R_SARSr_P8) is found in these forms (with max 5 mismatches):"))
table(sort(R_P8_seqs))


# gene S F1
myLen <- dim(S_F1)[1]
S_F1_seqs <- c()
for (i in c(1:myLen)){
  S_F1_seqs <- c(S_F1_seqs, as.character(subseq(mySet[S_F1$group[i]], start = S_F1$start[i], end = S_F1$end[i])))
}

print(paste("Pattern", Sgene_F1, "(Sgene_F1) is NOT found in", length(mySet)-length(S_F1_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", Sgene_F1, "(Sgene_F1) is found in these forms (with max 5 mismatches):"))
table(sort(S_F1_seqs))

# gene S P3
myLen <- dim(S_P3)[1]
S_P3_seqs <- c()
for (i in c(1:myLen)){
  S_P3_seqs <- c(S_P3_seqs, as.character(subseq(mySet[S_P3$group[i]], start = S_P3$start[i], end = S_P3$end[i])))
}

print(paste("Pattern", Sgene_P3, "(Sgene_P3) is NOT found in", length(mySet)-length(S_P3_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", Sgene_P3, "(Sgene_P3) is found in these forms (with max 5 mismatches):"))
table(sort(S_P3_seqs))

# gene S P4
myLen <- dim(S_P4)[1]
S_P4_seqs <- c()
for (i in c(1:myLen)){
  S_P4_seqs <- c(S_P4_seqs, as.character(subseq(mySet[S_P4$group[i]], start = S_P4$start[i], end = S_P4$end[i])))
}

print(paste("Pattern", Sgene_P4, "(Sgene_P4) is NOT found in", length(mySet)-length(S_P4_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", Sgene_P4, "(Sgene_P4) is found in these forms (with max 5 mismatches):"))
table(sort(S_P4_seqs))

# gene S R1
myLen <- dim(S_R1)[1]
S_R1_seqs <- c()
for (i in c(1:myLen)){
  S_R1_seqs <- c(S_R1_seqs, as.character(subseq(mySet[S_R1$group[i]], start = S_R1$start[i], end = S_R1$end[i])))
}

print(paste("Pattern", Sgene_R1, "(Sgene_R1) is NOT found in", length(mySet)-length(S_R1_seqs), "sequences (out of", length(mySet), ")."))
print(paste("Pattern", Sgene_R1, "(Sgene_R1) is found in these forms (with max 5 mismatches):"))
table(sort(S_R1_seqs))
