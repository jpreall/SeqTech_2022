### CSHL 2022 Advanced Sequencing Technologies  Course Materials

# Pandas Exercise
[Pancreas Atlas metadata](https://www.dropbox.com/s/jm1kg2x5u87w11e/metadata.csv.gz?dl=0)  

# Loupe interactive demo
[Mediocre PBMC data from SeqTech2017](https://www.dropbox.com/sh/qksaunln69yrqd1/AAAKLZ4E-yyfhb5-eYSnvnnZa?dl=0)  
[CRISPR KO A549 Dataset from 10X](https://www.dropbox.com/sh/z0h8nszrxcgjigx/AAD4Mgm_4-XNunVgT8tUdoBma?dl=0)

# Configuring your lab computer for Scanpy 
-------

### Running Scanpy in your AWS instance

In JupyerLab: File -> New -> Terminal
```bash
pip install 'scanpy[leiden]'
pip install harmony-pytorch gtfparse scrublet

cd workspace
```

Download the Jupyter Notebook we will be walking through:
```bash
wget http://34.224.25.10/workspace/SeqTech22_Data_Exploration.ipynb
```
Download the data we will be exploring:
```bash
wget http://34.224.25.10/workspace/data.tar.gz
tar -zxvf data.tar.gz
```
Relevant XKCD:  

<img src='https://imgs.xkcd.com/comics/tar.png'>
