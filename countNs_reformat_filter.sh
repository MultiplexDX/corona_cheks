# count sequences
count=$(grep -c ">" Spike.fa)

# 2 blocks below can be combined into one

# to get a report table with 5 columns
# prepare the header of the table
echo -e "ID\tEPI\tseqNs\ttechnology\tcountry" > Spike_ID_EPI_Ns.csv
for((i=2; i<$(($count+2)); i++))
do
ID=$(awk -v i=$i 'BEGIN{RS=">"; FS="\n"}; NR==i { print $1 }' Spike.fa | cut -d"/" -f3)
EPI=$(awk -v i=$i 'BEGIN{RS=">"; FS="\n"}; NR==i { print $1 }' Spike.fa | cut -d"/" -f4 | cut -d"|" -f2)
# count ambiguous bases (Ns)
# check if your aligning tool uses capital N or lower case n
# here, we count the lower case
seqNs=$(awk -v i=$i 'BEGIN{RS=">"}; NR==i { print ">"$0 }' Spike.fa | grep -v ">" | tr -d '\n' | grep -o "n" | wc -l)
# which sequencing technology is used? Usually reported in the 9th columns of the tsv gisaid annotation file
tech=$(grep -P $EPI"\t" gisaid_hcov-19.tsv | cut -f9)
# where was a sample collected?
Loc=$(grep -P $EPI"\t" gisaid_hcov-19.tsv | cut -f4 | cut -d"/" -f2)
echo -e $ID"\t"$EPI"\t"$seqNs"\t"$tech"\t"$Loc >> Spike_ID_EPI_Ns.csv
done

# redo the fasta file from block-like format into 2-lines format
for ((i=2; i<$(($count+2)); i++))
do
ID=$(grep ">" Spike.fa | sed -n ''$(($i-1))'p' | cut -d">" -f2)
echo ">"$ID >> Spike_2lines.fa
awk -v i=$i ' BEGIN{RS=">"; FS="\n"}; NR==i { print ">"$0 } ' Spike.fa | grep -v ">" | tr -d '\n' >> Spike_2lines.fa
echo "" >> Spike_2lines.fa
done


# filtering: if more 10 % of unknown bases (420 bases) in Spike (Wuhan Spike: 4209 bases = 1403 AAs)

# add Wuhan reference
awk ' BEGIN{RS=">"}; NR==2 { print ">"$0 } ' Spike.fa >> Spike_filtered.fa
# add Alpha (B.1.1.7; consensus 80 %)
awk ' BEGIN{RS=">"}; NR==3 { print ">"$0 } ' Spike.fa >> Spike_filtered.fa

for ((i=3; i<$(($count+2)); i++))
do
EPI=$(grep ">" Spike.fa | sed -n ''$i'p' | cut -d"/" -f4 | cut -d"|" -f2)
# use a list of sequences (EPI IDs) which you want to exclude (created in script plot_Ns.R)
myHit=$(grep -c $EPI samples_to_exclude.list | cut -d" " -f1)
if [[ $myHit == 0 ]]
then
awk ' BEGIN{RS=">"}; /'$EPI'\|/ { print ">"$0 } ' Spike.fa >> Spike_filtered.fa
fi
done


