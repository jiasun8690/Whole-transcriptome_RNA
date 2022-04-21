# Table of contents

- [LncRNA Filtering of Raw Reads](#L0)
  - [Alignment with Ribosome RNA (rRNA)](#L1)
  - [Alignment with reference genome](#L2)
  - [Transcripts Reconstruction](#L3)
  - [Novel Transcripts](#L4)
  - [Quantification of Transcripts](#L5)
  - [New lncRNAs identification](#L6)
  - [Differentially expressed transcripts (DEGs) Analysis](#L7)
  - [Enrichment Analysis](#L8)
- [miRNA Filtering of Raw Reads](#M0)
  - [Exist mirna statistics](#M1)
  - [Removing ncRNAs](#M2)
  - [Known mirna identification](#M3)
  - [Alignment with Reference Genome](#M4)
  - [Novel mirna identification](#M5)
  - [Differentially expressed transcripts (DEGs) Analysis](#M6)
  - [Target gene prediction](#M7)
  - [Enrichment Analysis](#M8)

# <span id='L0'>LncRNA Filtering of Raw Reads</span>

Fastp v0.18.0 (Chen, et al. 2018) was used for quality control of raw reads off the machine, and clean reads were obtained after filtering low quality data. Filtration parameters were as follows:

- the adapter for read1: -a AGATCGGAAGAGC 

- the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified： -q 20

- how many percents of bases are allowed to be unqualified (0~100): -u 50

- reads shorter than length_required will be discarded: -l 50

- if one read's number of N base is >n_base_limit, then this read/pair is discarded: -n 15

```shell
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i $input_1 -I $input_2 -o $output_1 -O $output_2 -j $output.json -h $output.html
```

## <span id='L1'>Alignment with Ribosome RNA (rRNA)</span>

Short reads alignment tool bowtie2 v2.28 (Langmead et al. 2012) was used for mapping reads to ribosome RNA (rRNA) database. The rRNA mapped reads were then removed. The reserved unmapped reads were used for subsequent transcriptome analysis. The comparison parameters were as follows:

- adopting local comparison method: --local
- Outputting reads that do not match:--un-conc-gz

```shell
bowtie2 --local -x $rRNA_index -1 $input_1 -2 $input_2 --un-conc-gz $output
```

## <span id='L2'>Alignment with reference genome</span>

HISAT2 v2.1.0 (Kim et al. 2015) was used for mapping analysis based on reference genome. The main parameters were as follows:

- Chain specific ratio:--rna-strandness RF 
- Clipping information in reference: --known-splicesite-infile $genome.splices_sites

```shell
hisat2 --rna-strandness RF -x $genome_index -q -1 $input_1 -2 $Input_2 --known-splicesite-infile $genome.splices_sites --dta --summary-file $output.stat --new-summary >$output.sam
```

## <span id='L3'>Transcripts Reconstruction</span>

The reconstruction of transcripts was carried out with software stringtie v1.3.4 (Pertea et al. 2015). The parameters were as follows:

- library type: stranded library fr-firststrand(--rf)
- minimum isoform fraction: -f 0.3
- The details of the other parameters can be found in: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual.

```shell
stringtie --rf -f 0.3 $input.bam -G $genome.gtf -o $output.gtf
```

stringtie merge

- minimum isoform fraction: -f 0.3
- minimum input transcript length to include in the merge: -m 200

```shell
stringtie --merge -f 0.3 -m 200 -G $genome.gtf -o $output.merge.gtf $input.gtf.list
```

## <span id='L4'>Novel Transcripts</span>

Cuffcompare was used to obtain the new gene type of Stringtie, and the self-written Perl script was used to extract class_code (u, i, j, x, c, e, o), and the transcript with exon number greater than 1 was the new transcript.
cuffcompare class_code: http://www.omicsclass.com/article/234

```shell
cuffcompare -r $genome.gtf -o $output $input.merge.gtf
```

## <span id='L5'>Quantification of Transcripts</span>

RSEM v1.2.19 (Li et al. 2011) software was called to quantify the transcript with script align_and_estimate_abundance.pl in Trinity v2.0.6. The parameters were as follows:

- library type: stranded library fr-firststrand(--SS_lib_type RF)
- Quantitative method:--est_method RSEM 
- Comparison methods：--aln_method bowtie2

```shell
perl align_and_estimate_abundance.pl --transcripts $exon.fa --gene_trans_map $exon.map --est_method RSEM --aln_method bowtie2 --seqType fq --left $input_1 --right $input_2 --output_dir $output --quality --phred33-quals --SS_lib_type RF
```

## <span id='L6'>New lncRNAs identification</span>

After obtaining the reconstructed transcript of Stringtie, the new lncRNA was predicted:

- Two softwares CPC2 (Yu-Jian et al. 2017) and CNCI (Liang et al. 2013) were used to assess the protein-coding potential of novel transcripts by default parameters. The intersection of both non protein-coding potential results were chosen as long non-coding RNAs.

CNCI parameters

- model types： -m v

```shell
## cpc
python CPC2 -i $input.fa --ORF -o $output.CPC2.result

## CNCI
python CNCI.py -f $input.fa -m ve -o $output
```

## <span id='L7'>Differentially expressed transcripts (DEGs) Analysis</span>

DESeq v1.20.0 (Anders et al, 2014) was used to performe lncRNAs differential expression analysis. BH FDR correction. The genes with the parameter of |log2(FC)| > 1,FDR< 0.05 were considered differentially expressed genes.

```shell
dds <- DESeqDataSetFromMatrix(count, colData, design= ~ condition)
dds <- DESeq(dds)
res= results(dds)
```

## <span id='L8'>Enrichment Analysis</span>

Enrichment analysis was performed using omicshare tools.
kegg: https://www.omicshare.com/tools/Home/Soft/pathwaygseasenior
go: https://www.omicshare.com/tools/Home/Soft/gogseasenior

# <span id='M0'>miRNA Filtering of Raw Reads</span>

Fastp v0.18.0 (Chen, et al. 2018) was used for quality control of raw reads off the machine, and clean reads were obtained after filtering low quality data. Filtration parameters were as follows:

- the adapter for read1: -a TGGAATTCTCGG 

- the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified： -q 20

- how many percents of bases are allowed to be unqualified (0~100): -u 50

- reads shorter than length_required will be discarded: -l 10

```shell
fastp -a TGGAATTCTCGG -q 20 -u 50 -n 15 -l 10 -i $input_1 -I $input_2 -o $output_1 -O $output_2 -j $output.json -h $output.html
```

## <span id='M1'>Exist mirna statistics</span>

Tags were compared with the precursors of mirbase miRNA sequences of this species using bowtie v1.1.2 software. The parameters were as follows:

- report end-to-end hits w/ <=v mismatches：-v 0 
- hits guaranteed best stratum; ties broken by quality: --best
- hits in sub-optimal strata aren't reported: --strata
- report all alignments per read: -a 
- do not align to reverse-complement reference strand: --norc

```shell
bowtie $mirna_hairpin -f $input -v 0 --best --strata -a --norc -S >$output.sam
```

## <span id='M2'>Removing ncRNAs</span>

- blastall v2.2.26 was used for mapping reads to Rfam database. The parameter was as follows: 
  - Expectation value: -e 0.01

```shell
blastall -p blastn -i $input -d $Rfam_db -e 0.01 -o $output
```

- blastall v2.2.26 was used for mapping reads to ncRNAs sequence of this species. The parameter was as follows:
  - Expectation value: -e 0.01

```shell
blastall -p blastn -i $input -d $Rfam_db -e 0.01 -o $output
```

## <span id='M3'>Known mirna identification</span>

Selecting to remove the sequence alignment with exist mirna, and then selecting the tags sequence after removing ncRNA. Bowtie v1.1.2 was used to align with the precursor sequences of mirbase mirna. The parameters were as follows:

- report end-to-end hits w/ <=v mismatches：-v 2 
- hits guaranteed best stratum; ties broken by quality: --best
- hits in sub-optimal strata aren't reported: --strata
- report all alignments per read: -a 
- do not align to reverse-complement reference strand: --norc

```shell
bowtie $mirna_hairpin -f $input -v 2 --best --strata -a --norc -S >$output.sam
```

## <span id='M4'>Alignment with Reference Genome</span>

Bowtie v1.1.2 was used to alignment with the tags of reference genome. The parameters were as follows:

- report end-to-end hits w/ <=v mismatches：-v 0 
- hits guaranteed best stratum; ties broken by quality: --best
- hits in sub-optimal strata aren't reported: --strata
- report up to <int> good alignments per read: -k 500

```shell
bowtie $genome -f $input -v 0 --best --strata -k 500 -S >$output.sam
```

## <span id='M5'>Novel mirna identification</span>

The tags sequences of exist mirna, known mirna were removed and aligned with the tags of genome. MiRDeep2 v2.0.0.7 was used to identify the novel mirnas. The parameters were as follows:

- maximum number of precursors to analyze when automatic excision gearing is used.: -g 50000

```shell
perl miRDeep2.pl $input $genome $aln.txt none none none $output -v -g 50000
```

## <span id='M6'>Differentially expressed transcripts (DEGs) Analysis</span>

edgeR v3.12.1 (Robinson et al, 2010) was used to performe miRNAs differential expression analysis. The parameter of |log2(FC)| > 1,FDR< 0.05 were considered differentially expressed mirnas.

```shell
Diff_list = estimateCommonDisp(Diff_list)
Diff_list = estimateCommonDisp(Diff_list)
Diff_list = estimateTagwiseDisp(Diff_list)
etest = exactTest(Diff_list)
```

## <span id='M7'>Target gene prediction</span>

Target genes were obtained by intersection of targetscan v7.0 and miranda v3.3a. 

miranda parameters：

- Set energy threshold to -E kcal/mol: -en -10 
- Demand strict 5' seed pairing: -strict

```shell
# targetscan
perl targetscan_70.pl $input_mirna $input_mrna $output

# miranda
miranda $input_mirna $input_mrna -en -10 -strict >$output
```

## <span id='M8'>Enrichment Analysis</span>

Enrichment analysis was performed using omicshare tools.
kegg: https://www.omicshare.com/tools/Home/Soft/pathwaygseasenior
go: https://www.omicshare.com/tools/Home/Soft/gogseasenior
