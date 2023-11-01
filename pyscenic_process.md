# Pyscenic  

* **Step 00 Create new environmnet**

```python
# download conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh 

# change env
#sudo apt install vim
vim ~/.bashrc 
export PATH="/home/data/t070406/miniconda3/bin"
source ~/.bashrc

# NO RUN
vim /home/data/t070406/miniconda3/envs/pyscenic/lib/python3.7/site.py
USER_SITE = "/home/data/t070406/miniconda3/envs/pyscenic/lib/python3.7/site-packages"
USER_BASE = "/home/data/t070406/miniconda3/envs/pyscenic"

python -m site

# conda configure
conda config --add channels r 
conda config --add channels conda-forge 
conda config --add channels bioconda
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/
conda config --set show_channel_urls yes

# create conda env
conda create -n pyscenic python=3.7 
conda activate pyscenic 

# install package
conda install pip
pip install scanpy -i https://mirrors.aliyun.com/pypi/simple/
pip install loompy -i https://mirrors.aliyun.com/pypi/simple/
pip install pyscenic -i https://mirrors.aliyun.com/pypi/simple/
```

* **Step 00 Extract scRNA-seq profile**

```R
library(Seurat)
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
AvailableData()
# InstallData("pbmc3k")
data("pbmc3k")
pbmc3k
#pbmc3k = readRDS("./pbmc3k.test.seurat.Rds")
#pbmc3k
write.csv(t(as.matrix(pbmc3k@assays$RNA@counts)),file = "for.scenic.data.csv")
```



* **Step 01 prepare file**

```
cat >01_cvs_to_loom.py
```

```python
# output the sample.loom file [input filename "pbmc3k.csv"]
import os,sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("pbmc3k.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs);
# ctrl+D
```



* **Step 02  fast pyscenic**

```python
cat >02_pyscenic_steps.bash
```

```python
# change to loom
python 01_cvs_to_loom.py 

# for human
dir=/home/data/t070406/Rdata/Case00_scanpy_test/pyscenic_reference #work directory
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl 
input_loom=./sample.loom

# step 01 grn
pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$dir/hs_hgnc_tfs.txt 

# step 02 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20  \
--mask_dropouts

# step 03 AUCcell
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 3
```

```
nohup bash 02_pyscenic_steps.bash 1>pySCENIC.log 2>&1 &
```

