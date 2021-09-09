# Nuc-PF
## Introduction

Nuc-PF is a collection of tools for training nucleosome footprint with multiple histone modification datasets to find potential pioneer factor regulations. We have used Nuc-PF successfully to find GATA2 as pioneer factor in LNCap cell lines. If you wish to focus on PF-associated states (downstream of Nuc-PF), we recommend that only the PF of interest be used, although Nuc-PF supports multiple PFs.

![Pipeline_V1](https://user-images.githubusercontent.com/17072186/111445001-eb978100-8745-11eb-9a9c-0afa8e3da56c.png)

## Quick Start & Dependencies

Nuc-PF is written in python as several scripts and should run on any x86 64-bit Linux system, however a computer with several gigabytes of RAM is strongly recommended to avoid an out-of-memory error.

__Required libraries:__  
  * **Numpy**:      http://www.numpy.org/  
  * **Scipy**:      http://www.scipy.org/  
  * **Pysam**:      https://github.com/pysam-developers/pysam  
  * **Cement**:     http://cement.readthedocs.org/  
  * **Treelib**:    https://github.com/caesar0301/treelib  
  * **Networkx**:   https://networkx.github.io/  

## Step 1: Nucleosome positioning and spacing detection.
MNase-seq data with MNase-ChIP-seq datasets were preferred for this pipeline (Corresponding to **Figure 1**).
We did the genome-wide detection of nucleosome positions from MNase-seq data, as well as identifies sharper nucleosome profiles for nucleosome spacing detection.
####
$ python pearson_correlation.py replicate_1 replicate_2 > Pearson_results

$ python NucDetect.py -i <input.bed> -t <S/P> -o <output_path/filename>  
 


Arguments |   ..  
 ---- | -----   
-i, --input | inputfile in sorted BED format.  
-o, --output | outputfile path and file name.  
-t, --type | S: single-end, P: pair-end, Default as "S".    

#### The nucleosome results columns are as following:
Chrid 	Start 	End	Index_No Length Height AUC Shape -log10(p-value)

## Step 2: Nucleosome grouping and Nuc_State Assignments.
MNase-ChIP-seq datasets were generated for nucleosome footprint detection and nucleosome states (Corresponding to **Figure 2**).

$ python NucGroup.py -g <Nuc_position_information.bed> -m [<histone_marker_1> <histone_marker_2>...] -o <Nuc_group_information.bed>.

$ python NucState.py -p <thread_number> -i <Nuc_group_information.bed> -o <State_assign_results>.

## Step 3: Transcript factor associated with dynamic nucleosome re-organization.
Potential pioneer factors with activation of treatment condition for detecing functional nucleosome regulators (Corresponding to **Figure 3**).

$ python 
$ python Plots.py <State_assign_results>

##  Step 4: ChIP-ePENS Border Calling Workflow 
Pioneer factor ChIP-ePENS analysis (Corresponding to **Figure 4 & Figure 5**)
#### Alignment by bowtie2
$ bowtie2 -v 3 -k 2 -m 1 -p 15 --fr -I 20 -X 400 -S /data/reference/hg19 -1 /data/ChIP-ePENs/GATA2_veh_ChIP-ePENs_R1.fastq -2 /data/ChIP-ePENs/GATA2_veh_ChIP-ePENs_R2.fastq /data/ChIP-ePENs/GATA2_veh_ChIP-ePENs.sam

#### Preparation for ePENs input
$ samtools view -bhS -q 30 /data/ChIP-ePENs/GATA2_veh_ChIP-ePENs.sam -o /data/ChIP-ePENs/GATA2_veh_ChIP-ePENs.bam 

#### ChIP-ePENS border calling by ePENs
$ python ePENS.py -D True -p 1e-8 -R 25  -t 12 -c 0.05 -k 2.0 -o ChIP-ePENs_results_GATA2_veh /data/ChIP-ePENs/GATA2_veh_ChIP-ePENs.bam

#### The results columns of ePENs could be explained by following:
Chrid 	Start 	End	Bordername 	Depth	Strand	Chernoff	Peakid	Compid	Pnb





### Contact Us
Tianbao Li: tianbaoli89@gmail.com; Victor Jin: jinv@uthscsa.edu
