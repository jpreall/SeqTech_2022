# SeqTech_2022
## CSHL Advanced Sequencing Technologies Course, Single Cell Dry Lab, 2022

### [Link to pre-baked data](DUMMY LINK) (XXXMB .tar.gz file)

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


XXXXXXXXX  

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

Let's dive in:

### Illumina sequencing output
*(This is taken care of for you this year.  But here is some useful information about this step anyway, in case it's ever your responsibility to do the FASTQ generation step.):*

**Congratulations!  You didn't screw up an experiment.  Now you might have data.**
