### CSHL 2022 Advanced Sequencing Technologies  Course Materials
# Configuring your cloud instance for Scanpy 
-------

### 1. Create a new Conda environment for Scanpy analysis

Open Terminal:
```bash
conda create --name scanpy -c conda-forge python=3.8 scanpy python-igraph 
conda activate scanpy
conda install jupyter
conda install -c bioconda gtfparse harmony-pytorch
```

### 2. Grab the data folder from Dropbox:
[Data Folder]([https:FIXME](https://www.dropbox.com/t/0nf8d4Wx48pCQpj7))  
Save this to your desktop

### 3. Launch the notebook
```bash
cd ~/Desktop/
jupyter notebook SeqTech22_Data_Exploration.ipynb
```
