# Nuc-PF
##Introduction

Nuc-PF is a collection of tools for training nucleosome footprint with multiple histone modification datasets to find potential pioneer factor regulations. We have used Nuc-PF successfully to find GATA2 as pioneer factor in LNCap cell lines. If you wish to focus on PF-associated states (downstream of Nuc-PF), we recommend that only the PF of interest be used, although Nuc-PF supports multiple PFs.

##Dependencies

Nuc-PF is written in python as several scripts and should run on any x86 64-bit Linux system, however a computer with several gigabytes of RAM is strongly recommended to avoid an out-of-memory error.

Required libraries

##Quick Start (Corresponding to Figure 1)
Step 1: Nucleosome positioning and spacing detection.
MNase-seq data with MNase-ChIP-seq datasets were preferred for this pipeline.










##  Border Calling Workflow (Corresponding to Figure 4 & Figure 5)
#### Alignment by bowtie2
$ bowtie2 -v 3 -k 2 -m 1 -p 15 --fr -I 20 -X 400 -S /data/reference/hg19 -1 /data/ChIP-ePEST/GATA2_veh_ChipePEST_L005_R1.fastq -2 /data/ChIP-ePEST/GATA2_veh_ChipePEST_L005_R2.fastq /data/ChIP-ePEST/vehLNCaP_GATA2_ChipePEST.R1R2.Paired.sam

#### Preparation for ePEST input
$samtools view -bhS -q 30 /data/ChIP-ePEST/vehLNCaP_GATA2_ChipePEST.R1R2.Paired.sam -o /data/ChIP-ePEST/vehLNCaP_GATA2_ChipePEST.R1R2.Paired.bam 

#### ChIP-ePEST border calling by ePEST
$python ePEST.py  -D True -p 1e-8 -R 25  -t 12 -c 0.05 -k 2.0 -o ePEST_vehFOXA1 /data/yez/Projects/ChIP-exo/vehLNCaP_FoxA1_ChipEXO.R1R2.Paired.Align.bam

#### The results columns of ePEST could be explained by following:
chrid 	start 	end	bordername 	depth	strand	chernoff	peakid	compid	pnb





### Contact Us
Tianbao Li: tianbaoli89@gmail.com; Victor Jin: jinv@uthscsa.edu
