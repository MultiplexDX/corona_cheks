# corona_cheks
Scripts can be used to filter and align sequences, find primers and probes (with particular focus on the Spike protein), and design of new diagnostic kits.

## Framework (only for small datasets, the parallelization is a must for bigger ones)
### Input
a) Omicron sequences from Gisaid.org in fasta format:<br/>
headers: ">virus/place/ID|EPI|date"<br/>
- to get the overview of countries sorted by the number of occurrences:
```bash
grep ">" gisaid_hcov-19_2021_12_04_11.fasta | cut -d"/" -f2 | sort | uniq -c | awk ' { print $1"\t"$2 } ' | sort -n -k1,1
```

b) A tsv annotation file of Omicron sequences (a file with information about sequencing technology etc.) 

### the whole genome alignment
1.) concatenate the reference sequence (Wuhan) to Omicron fasta file or any other sars-CoV-2 sequences:

```bash
cat Wuhan_reference.fa B.1.1.7_Consensus80.fa gisaid_hcov-19.fasta > Wuhan_B117cons80_omicron.fa
```

2.) use any multiple sequence alignment tool (mafft, clustal, ...):

```bash
mafft --thread 8 --auto Wuhan_B117cons80_omicron.fa > Wuhan_B117cons80_omicron_mafft.fa
```

3.) check MultiplexDX diagnostic kits
- all primers and probes in use: MultiplexDX_primers_probes_in_use_sars2_diagnostics.list
- use R script primers_probes_checks.R

### focus on Spike
4.) cut out the locus of Spike protein (see script dissect_Spike_fromAlignment.R) <br/>
In Wuhan reference, Spike starts "atgtttgtttttcttgttttattgccacta" and ends "ggagtcaaattacattacacataa".

5.) quality check
To count Ns (ambiguous bases) in every sequence, to create a 2-lines alignment format or to filter sequences, see:
countNs_reformat_filter.sh <br/>
To plot an overview of found Ns per sequence,  the N-counts dependence on sequencing technology or collection information, see:
plot_Ns.R

6.) a novel diagnostic kit or PAML analysis - an alignment of Spike sequences without codon splits (the fastest way for small datasets I know):

To realigned your filtered Spike sequences, use mafft:
```bash
mafft --thread 8 --auto Spike_filtered.fa > Spike_filtered_mafft.fa
```

- again create the 2-line format of your alignment (see the script countNs_reformat_filter.sh)
- load your block alignment into SeaView tool, switch on colours based on codons
- input for sed commands (see below): Spike_mafft_2lines.fa
- mark all false gaps and codon splits, fix them using sed commands (example: sed -i s/gtta------tctct/gttatc------tct/)
- SeaView also allows manual editing



