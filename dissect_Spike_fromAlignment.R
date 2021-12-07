setwd("my_path/")
library("Biostrings")
myMSA <- readDNAMultipleAlignment("Wuhan_B117cons80_Omicron450seqs_mafft.fa")
mySet <- DNAStringSet(myMSA)
WuhanStart <- matchPattern("atgtttgtttttcttgttttattgccacta", mySet$hCoV19_Wuhan_WIV06_EPI_ISL_402129)
WuhanEnd <- matchPattern("ggagtcaaattacattacacataa", mySet$hCoV19_Wuhan_WIV06_EPI_ISL_402129)

mySubset <- subseq(mySet, start = WuhanStart@ranges@start, end = WuhanEnd@ranges@start+WuhanEnd@ranges@width-1)
writeXStringSet(mySubset, filepath="./Spike.fa", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")