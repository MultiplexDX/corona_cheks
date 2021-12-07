setwd("my_path/")
# load in a report table
myTab <- read.table("Spike_ID_EPI_Ns.csv", sep="\t", header=TRUE)
dim(myTab)
mean(myTab$seqNs)
sd(myTab$seqNs)
# mean + sd
# cut-off: <10% (0.1*4209=421)
par(mar=c(1,3,2,0.5), mgp=c(1.75,0.5,0), cex.main=0.85, cex.axis=0.85, cex.lab=0.85)
boxplot(myTab$seqNs, ylab="Ns per genome", main="B.1.1.529 (omicron)\n450 sequenced genomes (5.12.2021, 1pm)")
abline(h=420, col="darksalmon", lwd=1.5, lty=2)
text(x=1,y=460,"cut-off = 10 % of the length (420 Ns)",cex=0.75, col="darksalmon")
text(x=1,y=380,"8.7 % of samples (39)",cex=0.75, col="darksalmon")
# Ns_in_samples_filteringCriteria.png

length(which(myTab$seqNs >= 421))
# 39 sequences to exclude

## check if the number of Ns differs among sequencing technologies
library(ggplot2)
# load information about Omicron sequences (from 3. to 452. line)
Nplot <- ggplot(data=myTab, aes(y=as.integer(seqNs), x=as.factor(technology), col=technology)) + geom_point() + ylab("number of unknown bases") + xlab("sequencing technology") + ggtitle("450 omicron sars2 sequences") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Nplot + geom_hline(yintercept=420, linetype="dashed", color = "red", size=2)
# Ns_per_technology.png

## check if the number of Ns differs among countries
NsCplot <- ggplot(data=myTab[c(3:452),], aes(y=as.integer(seqNs), x=as.factor(country), col=country)) + geom_point() + ylab("number of unknown bases") + xlab("country") + ggtitle("450 omicron sars2 sequences") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
NsCplot + geom_hline(yintercept=420, linetype="dashed", color = "red", size=2) + theme(legend.title = element_blank()) +  theme(legend.position = "none")
# UnknownBases_omicron_perCountry.png

ExcludeEPIs <- myTab$EPI[which(myTab$seqNs >= 420)]
# create a list of samples that need to be excluded
write(ExcludeEPIs, "samples_to_exclude.list", sep="\t")

