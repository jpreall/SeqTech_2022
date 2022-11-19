# SeqTech_2022
## CSHL Advanced Sequencing Technologies Course, Single Cell Dry Lab, 2022


### [Lecture Slides](https://www.dropbox.com/s/t1u9mogpsmjggjh/Preall_SeqTech_2022.pptx?dl=0) (112MB .pptx file)
### [Lecture Slides: Intro to scRNA-seq Analysis](https://www.dropbox.com/s/edkr5lgsbtscylp/Intro_to_scRNAseq.pptx?dl=0) (36MB .pptx file)
-------

This tutorial will be a guide through the first few steps of primary data analysis:
1. [FASTQ generation with `cellranger mkfastq`](#section1)
2. [Making a custom genome reference with `cellranger mkref`](#section2)
3. [Mapping and count matrix generation with `cellranger count`](#section3)
4. [Combining two samples into a shared, normalized matrix with `cellranger aggr`](#section4)

-------
## Single Cell Lab, Experiment A: Single-Cell RNAseq with Cell Hashing
[Cell Hashing Description](https://cite-seq.com/cell-hashing/)

<img src="https://citeseq.files.wordpress.com/2018/02/cell_hashing.png" width="700">

[Wet Lab Protocol](https://www.dropbox.com/s/mitbrqaxtgbavgo/SeqTech_2022_SingCell_protocol.docx?dl=0)  
For the scRNA-seq portion, we ran a simple experiment that will explore changes in gene expression as a function of cell culture density. Two different human cervical cancer lines, HeLa and Siha, were plated in 6-well format on Friday spanning a wide range of densities.  On Tuesday morning before the lab, we inspected these cultures and identified wells that roughly corresponded to a desirable range of confluencies (33%, 66%, 100%).  During lecture, we trypsinized these wells to produce a single-cell suspension.


| Cell Line | Cells seeded    | Confluency% | Hash Tag |
| --------- | :-------------: | :--------:  | :------: |
|HeLa       | 25k             | 33          | 1        |
|HeLa       | 50k             | 66          | 2        |
|HeLa       | 150k            | 100         | 3        |
|Siha       | 25k             | 33          | 4        |
|Siha       | 75k             | 66          | 5        |
|Siha       | 200k            | 100         | 6        |


Libraries were prepared according to the 10X Genomics User Guide and loaded onto a NextSeq2000 P3 flow cell with the following read lengths:  

| Read1 | Index1 | Index2 | Read2 |
|---|---|---|---|
|28bp|10bp|10bp|90bp|

Hash tags (a.k.a. Antibody-derived tags, or ADTs) were amplified in parallel, quantified, and spiked in at 5% molarity compared with the GEX library.

**Congratulations!  You didn't screw up an experiment.  Now you might have data.**

To demulitiplex, we are going to use 10X's built-in pipeline, designed for their cholesterol-modified oligo hash tagging system called <a href=https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cellplex>Cellplex</a>

Cell Ranger `count` doesn't support this natively, so we need to use a newer Cell Ranger function called `multi`.  This pipeline also runs all the functionality of `count` to generate a GEX matrix, but also demuxes.  To configure it, you need to build a configuration file pointing to all the necessary paths for FASTQs, references, and hash tag information.

```bash
cellranger multi \
  --id=SeqTech22_RNA_10k \
  --csv=multiplex_config.csv
```

**multiplex_config.csv**
```bash
[gene-expression],,,,,
reference,/seq/CellRanger/references/refdata-gex-GRCh38-2020-A,,,,
cmo-set,/mnt/grid/scc/data/Preall/SeqTech22/gex/hash_feature_barcodes.csv,,,,
expect-cells,12000,,,,
,,,,,
[libraries],,,,,
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
SeqTech22_RNA_10k,/mnt/grid/ngs/data/Elzar_Illumina/Illumina_runs/221114_VH00553_137_AAATLNYHV/AAATLNYHV/outs/fastq_path/SeqTech22/311810,any,SeqTech22_RNA_10k,gene expression,
SeqTech22_RNA_10k_ADT,/mnt/grid/ngs/data/Elzar_Illumina/Illumina_runs/221114_VH00553_137_AAATLNYHV/Data/Intensities/BaseCalls/311815/311813/,any,SeqTech22_RNA_10k_ADT,Multiplexing Capture,
SeqTech22_RNA_10k_ADT,/mnt/grid/ngs/data/Elzar_Illumina/Preall_Lab/221114_NB551387_0784_AHCYJJBGXN/311817/,any,SeqTech22_RNA_10k_ADT,Multiplexing Capture,

,,,,,
[samples],,,,,
sample_id,cmo_ids,description,,,
SeqTech22_10k_Hash1,Hash_1
SeqTech22_10k_Hash2,Hash_2
SeqTech22_10k_Hash3,Hash_3
SeqTech22_10k_Hash4,Hash_4
SeqTech22_10k_Hash5,Hash_5
SeqTech22_Hash6,Hash_6
```

**hash_feature_barcodes.csv**
```bash
id,name,read,pattern,sequence,feature_type
Hash_1,Hash_1,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GTCAACTCTTTAGCG,Multiplexing Capture
Hash_2,Hash_2,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGATGGCCTATTGGG,Multiplexing Capture
Hash_3,Hash_3,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TTCCGCCTCTCTTTG,Multiplexing Capture
Hash_4,Hash_4,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AGTAAGTTCAGCGTA,Multiplexing Capture
Hash_5,Hash_5,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AAGTATCGTTTCGCA,Multiplexing Capture
Hash_6,Hash_6,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GGTTGCCAGATGTCA,Multiplexing Capture
Hash_7,Hash_7,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGTCTTTCCTGCCAG,Multiplexing Capture
Hash_8,Hash_8,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CTCCTCTGCAATTAC,Multiplexing Capture
Hash_9,Hash_9,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CAGTAGTCACGGTCA,Multiplexing Capture
Hash_10,Hash_10,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,ATTGACCCGCGTTAG,Multiplexing Capture
```

### Demultiplexing
<img src=https://support.10xgenomics.com/img/multi_config_csv_expt_diagrams/multi_config_csv_gex_cmo.png align=left width=400>

