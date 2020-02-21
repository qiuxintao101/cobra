Welcome to the documentation of CoBRA
=======================
CoBRA - Containerized Bioinformatics workflow for Reproducible ChIP/ATAC-seq Analysis
=======================
Introduction to CoBRA:
CoBRA is a comprehensive ChIP/ATAC‐seq analysis tool built using snakemake and Docker which allows for escalable, reproducible, portable and easy-to-use workflows.

CoBRA combines the use of several dozen ChIP/ATAC‐seq tools, suites, and packages to create a complete pipeline that takes ChIP/ATAC‐seq analysis from unsupervised analyses, differential peak calling, and downstream pathway analysis. In addition, CoBRA has been outfitted with several recently published tools that allow for better normalziation and CNV correction. The results are compiled in a simple and highly visual report containing the key figures to explain the analysis, and then compiles all of the relevant files, tables, and pictures into an easy to navigate folder.

Table of Contents

*Quick Start and Installation

*Pipeline Details


* System Requirements
* Anatomy of a CoBRA PROJECT
* Getting Started
* * [Installing Docker]
* * [Running the CoBRA docker container]
* * [Downloading the example data]
* * [Test dry run on the workflow]
* * [Run the workflow]
* Setting up a Project for a CoBRA Run
* Set up data folder
* The Config File
* The Metasheet
* Running CoBRA
* Appendix A: Dana-Farber Members
* Appendix B: Specific CoBRA Commands for Replotting
* Setting up CoBRA for a group of users or server

System requirements:
Some of the tools that CoBRA uses, e.g. samtools and homer are memory intensive programs. Therefore we recommend the following system requirements for CoBRA:

Minimal system requirements:
We recommend that you run CoBRA on a server that has at least 30GB of ram. This will allow for a single-threaded run (on human samples).

Recommended system requirements:
We recommend that you have at least 128GB of ram and at least a 4-core CPU if you want to run CoBRA in multi-threaded mode (which will speedup the workflow significantly). Our own servers have 256GB of ram and 32 cores.

Anatomy of a CoBRA project:
All work in CoBRA is done in a PROJECT directory, which is simply a directory to contain a single CoBRA analysis run. PROJECT directories can be named anything (and they usually start with a simple mkdir command, e.g. mkdir CoBRA_for_thesis), but what is CRITICAL about a PROJECT directory is that you fill them with the following core components: (We first lay out the directory structure and explain each element below)

PROJECT/
Snakefile/
scipts/
config.yaml
metasheet.csv 
The 'scipts' directory contains all of the code. We'll explain how to download that directory below.

The config.yaml and metasheet.csv are configurations for your CoBRA run (also explained below).

After a successful CoBRA run, another 'analysis' folder is generated which contains all of the resulting output files.

Although included in this README are step-by-step instructions, it is assumed that the user has a basic understanding of the nix command line interface.

Getting started
Installing Docker:
We will be using the Docker to manage all of the software packages that CoBRA is dependent on.

Please follow the instructions on the following website to install docker Community Edition (CE): https://docs.docker.com/install/

Running the CoBRA docker container:
We are now ready to use docker command to run the CoBRA.

>   docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && bash /init.sh"
NOTE: The "/directory_on_your_system/data" refers to the folder in your system, "/data" refers the folder inside the container

NOTE: If you have not run CoBRA before, this command will download the CoBRA container and run it. After running the command, we should have all codes that needed for this workflow.

Downloading the example data:
>   docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && bash /example.sh"
NOTE: This will download an example data (ENCODE GR ChIP-seq from different dexamethasone concentration treatment)for the test run.

Test dry run on the workflow:
>   docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all --cores 2 -np"
NOTE: If all needed code and files are ready, you should see the workflow printing the jobs that is planed by the workflow.

Run the workflow:
>   docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all --cores 2"
#Setting up a Project for a CoBRA Run

Setting up your data directory:
As mentioned above, we highly recommend that you pool your raw data into a directory called 'bam_bed_bw' (within PROJECT).

IF all of your data is already centrally stored in a directory called '/some/path/to/my/data/', then symbollically link your data to your PROJECTS folder.
cd PROJECT
ln -s /some/path/to/my/data/ ./bam_bed_bw
The 'ln -s' command creates a symbolic link from '/some/path/to/my/data/' and names it 'data'. Symbolic Links have been shown to be adequate when testing, if not, then please try one of the other solutions.

IF your files are not centrally stored:
- make a directory within PROJECT called 'data' and copy your raw files into data

What your PROJECT directory should look like (up to now):
PROJECT/
CoBRA/
bam_bed_bw/ config.yaml
metasheet.csv 
NOTE: you will have to setup your PROJECT directory for each CoBRA run.

Configuring the META files: config.yaml
The config.yaml file has three main sections. PATHS, PARAMS, SAMPLES:

PARAMS:
This section holds parameters specific to your project design

######Project Name ######Use in pca, sample-sample, sample-feature plot. Please use "_" to seperate different words project: GR_ChIP_seq

######enhancer option enhancer/promoter/all enhancer: all

######Location of metasheet metasheet: metasheet.csv ref: "scripts/ref.yaml"

######Assembly is needed when seperate enhancer/promoter, motif finding, nearyby gene assembly: hg19

######At least mini_num_sample should have RPKM > rpkm_threshold rpkm_threshold: 1 mini_num_sample: 0

######Scale method for the nomalize counts among samples ######z- z-score ######q- quantile-normalize ######l- log-transform scale: q

######Fliter metric in feature selection ######sd- Standard deviation ######cov- Coefficient of Variation ######av- mean filter-opt: cov

######top percent cutoff filter-percent: 100

######limited of peaks to use for plot SSpeaks: 20000000

SFpeaks: 20000000

num_kmeans_clust: 6

######DEseq_cut_off - Padj/LG2FC Padj: 0.05 LG2FC: 0

######Motif analysis - true/false motif: 'true'

######BAM files sorted? true/false bam_sort: 'true'

######CNV correction? true/false CNV_correction: 'false'

######unchanged heatmap unchanged_heatmap: 'false' 
SAMPLES:
In this section of the configuration file, you specify the NAMES of each sample, and the PATHS to the sample's raw data. Raw data files can either be fastq, fastq.gz, or bam formated files.

As recommended above, if all of your raw data are located in PROJECTS/data, then each path will simply start like:
'bam_bed_bw/XX1.bed'

If you did not follow the recommended best practice then you will have to specify the full paths here.

Each sample should be given a NAME (arbitrary text) and a PATH

EXAMPLE:

bed:
   sample1: ./XX1.bed
   sample2: ./XX2.bed


samples:
   sample1: ./XX1.bam
   sample2: ./XX2.bam


bigwig:
   sample1: ./XX1.bw
   sample2: ./XX2.bw
IMPORTANT: You cannot mix Paired-end and Single-end samples within the same CoBRA run as this will cause an ERROR. If necessary, run all of one type of data first, followed by the net type of data after.

Configuring the META files: metasheet.csv
Make the metasheet file in excel, and save it as a .txt or .csv, It doesn’t matter what it is named as long as it is called in the config in the spot marked “metasheet,” see the config section if confused. The format should be something like the following:

Sample	Cell	Condition	Treatment	Replicates	comp_M7_DOX_over_NoDox	comp_T47D__DOX_over_NoDox
A1	MCF7	Full_Media	NoDOX	1	1	
A2	MCF7	Full_Media	NoDOX	2	1	
B1	MCF7	Full_Media	DOX	1	2	
B2	MCF7	Full_Media	DOX	2	2	
C1	T47D	Full_Media	NoDOX	1		1
C2	T47D	Full_Media	NoDOX	2		1
D1	T47D	Full_Media	DOX	1		2
D2	T47D	Full_Media	DOX	2		2
The first column should always be sample names that exactly match the sample names used in config.yaml (see SAMPLES just above)
The samples that you want to perform a Differential Peak Calling (DE) on using limma and deseq should be marked by the “comp” columns more on this below
This is important! The “control” should be marked with a 1, and the “treatment” should be marked with a 2.
It is recommended that if you should have a “replicates” column to denote different samples, it is a good idea to not only have each of the sample names be unique, but also make sure that the associated metadata is unique as well to each sample.
The rest of the metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many “comp_XXXX” columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates)

Again, make this in excel so that all of the spacing is done correctly and save it out as a .txt or .csv file. This is the most common bug, so please follow this.

Common Problems with metasheet
Characters to avoid: ("-", "(", ")", " ", "/", "$") To avoid bugs, the only punctuation that should be used is the underscore “_”. Dashes, periods, etc, could cause a bug because there is a lot of table formatting and manipulation, or they are invalid characters in R. NOTE: CoBRA parses the meta file and will convert MOST of these invalid characters into '.'--dollarsigns will just be dropped. The CoBRA parser will also convert between dos/mac files to unix format.
It is very important that you know that samples A is what you mark with 1, and samples B is what you mark with a 2. You should name your output following this format as well "comp_cond_AvB” This will let the reader know what the output DE files refer to.
Deseq: ”baseMeanA” refers to samples A, which follows condition 1 and “baseMeanB” refers to samples B which follows condition 2. logfc is B/A
Limma: Logfc refers to B/A
Running CoBRA:
Now that we have setup our PROJECT directory (downloading the 'CoBRA' code directory, creating our 'data' directory, and configuring our config.yaml and metasheet.csv), we are (finally!) ready to run CoBRA.

Next we will perform a DRY-RUN to make sure that we setup the CoBRA PROJECT directory correctly. In your PROJECT folder run the following command:

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all --cores 2 -np"

This will return a large output which basically outlines what CoBRA is about to do. If no errors come back, then you will mostly see GREEN and YELLOW print outs. If there are errors, then you will see some RED print outs.

If there are no errors, then use the following command to run CoBRA:

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all --cores 2"

If there are errors, try to see what the error is about. Was it a mistyped path? Etc. If all else fails, email the CoBRA team (email address needed)

APPENDIX B: Specific Replotting:
After you have run CoBRA in its entirety, you may want to go back and tweak your outputs. Maybe adding or subtracting metadata columns, differential expression columns, or maybe just doing a subset of your data. Below is a list of snakemake commands to run CoBRA to rerun some specifics for further downstream analysis tweaking. Note that CoBRA is built to automatically rerun downstream analysis if you adjust the config or the metasheet.

To learn about how snakemake works, and some of the specifics of the following commands and others, look into the snakemake documentation

The following are some useful commands for rerunning and adding to the download analysis without having to rerun the whole pipeline:

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake -np"

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all --cores 2"

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all heatmapSF_plot -f "

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all heatmapSS_plot -f "

docker run --rm -v /directory_on_your_system/data:/data -t cfce/cobra /bin/bash -c "cd data && source activate CoBRA && snakemake all pca_plot -f "

Adding comp columns will automatically make it generate new differential expressions analysis and adjust figures accordingly.
"Touching" the metasheet will have CoBRA rerun all downstream output.

touch metasheet.csv

.. toctree::
   :maxdepth: 2
   :caption: Contents

   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* :ref:`Q&A`
