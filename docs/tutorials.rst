
.. _docs-tutorial:

********
Case Studies
********

The three case studies below gives you a taste of different capabilities of our CoBRA workflow. The data sets you will be exploring are, 
1) a glucocorticoid receptor (GR) ChIP-seq data set from the ENCODE project, 
2) a set of H3K27ac ChIP-seq data from colon cancer cell lines, and
3) an ATAC-seq experiment on HL-60 promyelocytes differentiating into macrophages. 

Each of the case studies wiil demonstrate a subfield of expertise of our CoBRA pipeline. 

.. note::  Our *CoBRA* pipeline is implemented using the *Snakemake* workflow management system. We recommend getting to know the basics of *Snakemake* prior to trying this tutorial, as it helps with basic troubleshooting and solving common errors associated with running *CoBRA*. The *Snakemake* documentation and tutorial page can be directed through `here <https://snakemake.readthedocs.io/en/stable/index.html>`_.


Setup
=====

In preparation for the tutorials, please use the following steps to set up the cobra environment and retrieve the latest version of our pipeline:

1. **Initiate Docker Container**: 
  
  Use the following command to start the cobra container in an empty working folder:
  
  .. code-block:: Bash

    docker run --rm -v $PWD:/cobra -it cfce/cobra:2.0

2. **Retrieve the Latest Version of Cobra Pipeline:**

  If installed using docker, run the following command to change the working directory. Otherwise, skip to next command:
   
  .. code-block:: Bash
   
     cd cobra && source activate cobra 
   
  pull the latest version of CoBRA using ``git clone`` :

  .. code-block:: Bash

     git clone https://bitbucket.org/cfce/cobra.git .

  If you receive an error, *Git* may not be installed on your system. Please consult the internet on how to best install Git for your system.



Case Study 1: GR ChIP-set Data Set
================

Background
**********
This tutorial makes use of a publicly available glucocorticoid receptor (GR) ChIP-seq data (`GSE32465 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32465>`_) from a lung adenocarcinoma cell line (A549) at 3 different concentrations of dexamethasone, a potent GR agonist. In the experiment, samples were treated with 0.5nM, 5nM, or 50nM of dexamethasone. 


Download and set-up for running the GR_ChIP sample dataset
**********************************************************

  This dataset is of moderate size and may take 5-10 minutes to download. It contains the data files, as well as the config files - ``config.yaml`` and ``metasheet.csv`` - filled with correct parameters. One advantage of *CoBRA* is that it accommodate the use of different input files; it may take either 1) the set of ``.bam``, ``.bed``, and ``.bigwig`` files for each sample, or 2) a single ``.fastq.gz`` file for each sample. We have two separate downloads prepared this specific dataset, one containing the bam, bed and bigwig files as input, the other containing fastq files as input. **Choose only one to download.** The analysis process is exactly the same.

  .. code-block:: Bash
   
     snakemake download_example_GR_ChIP
                   #OR
     snakemake download_example_GR_ChIP_fastq 
  
  When the data set is downloaded, we can proceed to set up for the run. Usually for running CoBRA on a new experiment, the two config files ``config.yaml`` and ``metasheet.csv`` would need to be set up acccordingly. In this tutorial, they have been filled already. 
  
  .. note::  In ``config.yaml``, the parameter `motif` has been set as `true` to perform motif enrichement and clustering analysis. The DEseq normalize method parameter `norm_method` was set as `depth` to opt for normlization by the sequence depth of each sample.

  To check if the setup is correct, begin a dry run via the following command:
  
  .. code-block:: Bash

     snakemake all -np

  As seen below, the ``-np`` command of *Snakemake* outputs the execution plan of the run instead of actually perform the steps. It produces a job count list, that is, a list of all the snakemake rules that will be run to achieve the outputs, and a summary for each snakemake rule including the rule name, input, and output. 
  
  .. code-block:: shell-session            
     
     #Sample Job Summary 
     $ snakemake all -np
     Job 81: ALIGN: Running BWA mem for alignment
     bwa mem -t 8 ref_files/hg19/bwa_indices/hg19/hg19.fa /mnt/cfce-stor1/home/xq08/Projects/Diff_Peak_Methods_Investigation/FASTQ_files_GR_ENCSR989EXF/dexamethasone_at_500pM/ENCFF000NBL.fastq.gz | samtools view -Sb - > analysis/preprocessed_files/align/0.5nM_Dex_1/0.5nM_Dex_1.bam 2>>analysis/logs/align.log
     
  .. note:: Once a run is initiated on the Docker container, please DO NOT exit the container while the run is still ongoing. This would result in **interruption of the current `CoBRA` run**. 


Quick One-Step Analysis
**********************************************************

  Once the dry run completes without errors, run the pipeline using the following command (using 6 cores).

  .. code-block:: Bash

     snakemake all --cores 6

  Then wait for the result to come out in a few hours. It is plain and simple!


Step-By-Step Analysis
**********************************************************

  While the CoBRA pipeline is designed to be fast and efficient, easily-excuetable with just a few lines of commands, it is possible to produce the analysis in a step-wise fashion by running specific parts of the pipeline.

1. **Unsupervised Analysis - PCA Plot**: 

    .. code-block:: Bash

       snakemake pca_plot -f
  
  This command produces the ``pca_plot_100_percent.pdf`` file located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. The first page of the file is a color-coded Principal component analysis (PCA) plot that depicts how samples are separated in the first two principal components (those with the largest variance). The second page includes a scree plot indicating the percentage of variance captured by each principal component.


  .. figure:: ./tutorial_figures/1_pca.png
      :scale: 20 %
      :alt: case 1 pca plot
      :align: center
      
  As illustrated in the PCA plot, PC1 separates the samples with different treatment concentration of dexamethasone, while PC2 further separates the sample replicates.
 
  .. figure:: ./tutorial_figures/1_pca_scree.png
      :scale: 20 %
      :alt: tutorial 1 pca scree
      :align: center

  As illustrated in the PCA plot and scree plot above, PC1 (capturing 40.8% of variance explained) separates the samples with different treatment concentration of dexamethasone - namely 0.5nM from 5nM and 50nM, while PC2 (18.7% variance) further separates the sample replicates.


2. **Unsupervised Analysis - Sample-Sample Correlation Plot**: 

    .. code-block:: Bash

       snakemake heatmapSS_plot -f
  
  This command produces the ``heatmapSS_plot_100_percent.pdf`` file located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. It provides information on the clustering result based on the Pearson correlation coefficient, and illustrates the similarity between all samples in a pairwise fashion.
  
  .. figure:: ./tutorial_figures/1_SS.png
      :scale: 20 %
      :alt: case 1 ss heatmap
      :align: center
      
  As illustrated in the Sample-Sample correlation plot, samples replicates cluster tightly together (r > 0.6). And samples treated with 0.5nM of dexamethasone exhibited to be far different from samples treated with 5nM or 50nM dexamethasone.


3. **Supervised Analysis - DeSeq2 Differential Peak Analysis**: 

  The key inquiry to be satisfied for any ChIP-seq/ATAC-seq analysis is what the differential sites are between sample groups of interest. In *CoBRA*, this analysis is done by incorporating differential peak callin gby DESeq2 while using sequencing depth as a scale factor, and thus significantly reducing false positive differential peak-calling.
  
    .. code-block:: Bash

       snakemake run_limma_and_deseq -f
  
  This command produces a series of files located in the ``analysis_result/differential_peaks/c50nm_vs_0.5nm`` folder, including the following:
    - ``c50nm_vs_0.5nm.deseq.csv``: a differentail peaks analysis table produced by DESeq2
    - ``c50nm_vs_0.5nm.deseq.Padj0.05.LG2FC.0.up.bed`` and ``c50nm_vs_0.5nm.deseq.Padj0.05.LG2FC.-0.down.bed``: bed files of peaks that are differentially up- and down-regulated, respectively
    - ``c50nm_vs_0.5nm.deseq.sum.csv``: a table including total number of differential peaks under different thresholds
    - ``c50nm_vs_0.5nm.t.test.csv``: a t-test table of the differential peaks
    - ``MA_plot.pdf``: a MA plot comparing the two treatment samples
  
  DEseq2 by default normalizes all samples by total reads in the read count table. In contrast, in the GR ChIP-seq experiment, samples treated with 50nM dexamethasone exhibit much greater GR binding and the FRiP score is higher than samples treated with 0.5nM (9.3 vs 0.9). Therefore, DESeq2’s normalization method decreases the peak intensity in the 0.5nM treated samples because the FRiP scores are higher in the 50nM sample resulting in false positive differential peaks. In *CoBRA*, we use a scaling factor dependent on the sequencing depth of each sample. This eliminates the false positive downregulated peaks called by DESeq2 using the default scaling factor. 
  
  However, normalizaiton by default DESeq2 method is still included as an option in our pipeline, see :ref:`parameter_norm_method` for detail.
  
  .. figure:: ./tutorial_figures/1_maplot.png
      :scale: 35 %
      :alt: case 1 ma plot
      :align: center
  
  The MA plot above shows that the 50nM treatment samples have significant numbers of upregulated peaks called by DESeq2 and no downregulated peaks.
  
  Intensity measurement of the differnetial peaks can be done using the following command
  
    .. code-block:: Bash

       snakemake run_deeptools_diff_peaks -f
  
  It produces ``c50nm_vs_0.5nm.deseq.Padj0.05.LG2FC.0.pdf`` which illustrates the peak intensity of the differentially up and downregulated peaks. 

  .. figure:: ./tutorial_figures/1_peaks.png
      :scale: 30 %
      :alt: case 1 diff peats
      :align: center
       
  The peak-intensity heatmap above further illustrates that there only exist differentially upregulated peaks in 50nM treatment samples as compared to 0.5 nM dexamethasone treated samples, and intensity goes as high as 1.75.


4. **Comparison of Up and Down-regulated Site: Cistrome Toolkit**: 

  *CoBRA* has a built-in feature that compares up and down-regulated sites to a comprehesnive database of ChIP/ATAC and DNase data, and outline a series of most similar samples in terms of genomic interval overlaps with the differential sites located in the (`Cistrome database <http://cistrome.org/db/#/>`_). This feature allows researchers to pin-point those similar data set of interest and download for further investigation. It can provide unique insight into gained or lost sites such as identifying which transcription factor potentially binds to a differential peak set after a perturbation and in investigating similar cellular systems.
  
    .. code-block:: Bash

       snakemake run_cistrome_toolkit -f
  
  Using the command above, *CoBRA* outputs a series of files located in the ``analysis_result/differential_peaks/c50nm_vs_0.5nm/cistrome_toolkit`` folder, including:
    - a plot of most similar samples ranked by their giggle score, and
    - two tables of cistrome toolkit result, each include a list of GEO accession numbers corresponding to all ChIP-seq data with similarity to the differential peak set (up or down-regulated)
    
  .. figure:: ./tutorial_figures/1_cistrome_geo.png
      :scale: 30 %
      :alt: case 1 cistrome GEO accession table
      :align: center
      
      The Cistrome Toolkit result table would include Cistrome DB sample ID, GEO accession number (GSM) and key information about the data set, i.e. factor name, cell line, cell type, giggle score. The entries are ranked by their giggle score.
  
  .. figure:: ./tutorial_figures/1_cistrome.png
      :scale: 25 %
      :alt: case 1 cistrome result
      :align: center

  As show in the plot above, for the gained GR binding sites in the dexamethasone treatment, the NR3C1 factor in Lung is the most similar ChIP-seq in the Cistrome database to this GR data set.


Case Study 2: MSS and MSI Colorectal Cancers ChIP-seq Data Set
================

Background
**********
This tutorial makes use six samples from several experiments: three Microsatellite Instable (MSI) samples and three Microsatellite Stable (MSS) samples (Tak et al. 2016; Piunti et al. 2017; Piunti et al. 2017; Maurano et al. 2015; McCleland et al. 2016; Rahnamoun et al. 2018). Microsatellite Instable (MSI) and Microsatellite Stable (MSS) are two classses used to characterize colorectal cancers. MSS tumors are one of the most highly mutated tumor types (Taieb et al. 2017) and exhibit a high copy number variations. Without adjustment, a differential peak caller will rank peak loci with high copy number gain in MSS as being the most differential compared to MSI. To observe differential peaks between the MSI and MSS samples, *CoBRA* allows for **copy number variation adjustment** during the supervised analysis.


Download and set-up for running the MSS_MSI sample dataset
**********************************************************

  Please use the following command to download the MSS_MSI ChIP-seq sample dataset. 

  .. code-block:: Bash
   
     snakemake download_example_MSS_MSI
  
  When the data set is downloaded, we can proceed to set up for the run. 
  
  .. note::  In ``config.yaml``, the parameter `cnv` has laid out a path for **CNV files** (usually in ``.igv`` format) corresponding to each sample. For this dataset, the copy number was called on the ChIP-seq data itself using CopywriteR (Kuilman et al. 2015) with IP but can also be done using qDNAseq (Scheinin et al. 2014) with the input control if available. Any other source of CNV data can also be used when put in a standard format. See details in :ref:`section_cnv` for how to prepare the files for CNV analysis to be listed in the ``config.yaml``.
  
  .. note::  In ``metasheet.csv``, there is column ``comp_MSS_vs_MSI`` that compare all 3 samples of MSS to all 3 samples of MSI, and column ``comp_MSS2_vs_MSI2`` that only compares only 2 samples from each group. This demonstrates *CoBRA*'s feature of setting multiple independent differential expression analysis. See details in :ref:`section_metadata`.

  To check if the setup is correct, begin a dry run via the following command:

  .. code-block:: Bash

     snakemake all -np


Quick One-Step Analysis
**********************************************************

  Once the dry run completes without errors, run the pipeline using the following command (using 6 cores).

  .. code-block:: Bash

     snakemake all --cores 6

  
Step-By-Step Analysis
**********************************************************

1. **Unsupervised Analysis - PCA Plot, Sample-Sample Correlation Plot, etc.**: 

    .. code-block:: Bash

       snakemake pca_plot -f
       snakemake heatmapSS_plot -f
  
  As demonstrated in the previous case study, these command produces the pca plot and the heatmaps located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. 

  .. figure:: ./tutorial_figures/2_pca.png
      :scale: 20 %
      :alt: case 2 pca plot
      :align: center

  As illustrated in the PCA plot above, PC1 (capturing 44.5% of variance explained) clearly separates the MSS samples (colored in turquois) and MSI samples (colored in pink).

  
  .. figure:: ./tutorial_figures/2_SS.png
      :scale: 20 %
      :alt: case 2 ss heatmap
      :align: center
 
  The Sample-Sample Correlation shows clearly that the MSS samples cluster together, and the same applies to the MSI samples. And the two sample groups exhibit little correlation. 


2. **Supervised Analysis - DeSeq2 Differential Peak Analysis**: 

    .. code-block:: Bash

       snakemake run_limma_and_deseq -f
       snakemake run_deeptools_diff_peaks -f
  
  As demonstrated in Case Study 1, these command produces a series of differential peak analysis results located in the ``analysis_result/differential_peaks/MSS_vs_MSI`` folder, including a MA plot and a peak intensity plot. Applying copy number variation adjustment eliminates false positive peaks that would otherwise be called as differential due to their significant copy number difference between the two sample groups MSI and MSS.

  .. figure:: ./tutorial_figures/2_peaks.png
      :scale: 30 %
      :alt: case 2 diff peaks
      :align: center
      
      Peaks Intensity Plot with CNV Adjustment
     
  .. figure:: ./tutorial_figures/2_peaks_nocnv.png
      :scale: 30 %
      :alt: tutorial 2 diff peaks no cnv
      :align: center
      
      Peaks Intensity Plot with No CNV Adjustment
  
  Comparing the two peak intensity heatmaps above, differential peaks in the plot generated with CNV adjustment generally shows in general higher intensity.
  
  
3. **GSEA**: 

  *CoBRA* has built-in features to do the Gene Set Enrichment analysis, which is performed on the ranked list of genes produced by the pipeline.
  
    .. code-block:: Bash

       snakemake run_GSEA -f
  
  Using the command above, *CoBRA* outputs a series of GSEA analysis results in ``analysis_result/differential_peaks/MSS_vs_MSI/GSEA`` folder, including:
    - ``index.html``: summary report for the GSEA
    - ``gsea_report_for_na_neg`` and ``gsea_report_for_na_pos``: summary report including all ranked genes sets and their statistics 
    - ``neg_snapshot.html`` and ``pos_snapshot.html``: snapshots of all enrichment plots of enriched gene sets curated
    - ``enplot_{Gene_Set}``: individual enrichment plots of an enriched gene set
    - ``{Gene_Set}.html`` and ``{Gene_Set}.xls``: individual GSEA Results Summary of an enriched gene set
  
  .. figure:: ./tutorial_figures/2_gsea_nocnv.png
      :scale: 60 %
      :alt: case 2 GSEA
      :align: center
      
      Enrichment Plot with No CNV Adjustment
      
  .. figure:: ./tutorial_figures/2_gsea_cnv.png
      :scale: 60 %
      :alt: case 2 GSEA
      :align: center
      
      Enrichment Plot with CNV Adjustment
  
  Without CNV adjustment, GSEA will indicate greatest enrichment in gene sets solely related to amplification. As a result, it is challenging to assess the true epigenetic differences between the two colorectal cancer types. MSS vs MSI type tumors presents an especially challenging scenario. The MSS tumors exhibits large scale high copy number variations across the genome, including the 8q arm. However, the MSI tumors exhibits a focal amplification directly at 8q12-q22 region, making it very difficult for regular DE pipelines to assess the difference between these two types of amplifications. *CoBRA* is able to distinguish that difference by CNV adjustment and demonstrate in the GSEA result.
  
  The gene set NIKOLSKY_BREAST_CANCER_8Q12_Q22_AMPLICON includes genes up-regulated in non-metastatic breast cancer tumors with amplification in the 8q22 region. Without adjustment for copy number variation, this gene set is significantly enriched in MSS samples, with a normalized enrichment score of -2.84 and an adjusted p-value less than 0.0001. With CNV adjustment, this gene set is considered far less enriched, with a normalized enrichment score of -1.32 and an adjusted p-value of 1.
  

Case Study 3: ATAC-seq from HL-60 promyelocytes differentiating into macrophages
================

Background
**********
This tutorial makes use of ATAC-seq from HL-60 promyelocytes differentiating into macrophages (`GSE79019 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79019>`_). The samples were taken utilized a five-day time course (0hr, 3hr, 24hr, 96hr, and 120hr) to profile accessible chromatin of HL-60 promyelocytes differentiating into macrophages. Here *CoBRA* results shows investigation of the differentiation of macrophages through changes in the landscape of accessible chromatin. 


Download and set-up for running the Macrophage_atac sample dataset
**********************************************************

  Please use the following command to download the Macrophage ATAC-seq sample dataset. 

  .. code-block:: Bash
   
     snakemake download_example_Macrophage_atac
  
  When the data set is downloaded, we can proceed to set up for the run. 
  
  .. note::  In the ``metadata.csv``, a couple of different comparison columns were set up in order to do pair-wise comparison of samples taken from different time point. This is another efficent feature of *CoBRA* - allowing for multiple differential expression analysis done separately. For each comparison, a complete set of supervised analysis results (motif analysis, cistrome toolkit, GSEA) will be completed in the respective subfolder under ``analysis_result/differential_peaks``. See details in :ref:`section_metadata` for how to prepare ``metadata.csv`` for multiple comparisons.
  
  .. note::  In the ``config.yaml``, the parameter `percent` has been set to 10, indicating that only top 10% peaks will be used in the unsupervised analysis and clustering analysis. The `rpkm_threshold` and `mini_num_sample` can also be adjusted accordingly to different data sets. See details in :ref:`configurationFile` for how to set those parameters. 


  To check if the setup is correct, begin a dry run via the following command:

  .. code-block:: Bash

     snakemake all -np

Quick One-Step Analysis
**********************************************************

  Once the dry run completes without errors, run the pipeline using the following command (using 6 cores).

  .. code-block:: Bash

     snakemake all --cores 6

  
Step-By-Step Analysis
**********************************************************

1. **Unsupervised Analysis - PCA Plot, Sample-Sample Correlation Plot, Sample-Feature Heatmap, etc.**: 

    .. code-block:: Bash

       snakemake pca_plot -f
       snakemake heatmapSS_plot -f
  
  Like illustrated in Case Study 1, this command produces the pca plot and the heatmaps located in the ```analysis_result/clustering_analysis/rpkm.3_num_sample.2_scale.q_fliter.cov.10/plots`` folder. 

  .. figure:: ./tutorial_figures/3_pca.png
      :scale: 20 %
      :alt: case 3 pca plot
      :align: center

  As illustrated in the PCA plot and scree plot above, PC1 (capturing 57=0.7% of variance explained) clearly separates the samples by their time frame
  
  .. figure:: ./tutorial_figures/3_SS.png
      :scale: 20 %
      :alt: case 3 ss heatmap
      :align: center


2. **Unsupervised Analysis - Sample-Feature Heatmap**: 

    .. code-block:: Bash

       snakemake heatmapSF_plot -f
  
  This command produces the ``heatmapSF_plot_10_percent.pdf`` file located in the ``analysis_result/clustering_analysis/rpkm.3_num_sample.2_scale.q_fliter.cov.10/plots`` folder. It illustrates clustering of samples based on correlation on the horizontal axis and clustering of peaks on the vertical axis. It presents patterns of peaks (by k-means clustering) across samples and identifies the clusters that are enriched in a subset of samples.
  
  .. figure:: ./tutorial_figures/3_SF.png
      :scale: 20 %
      :alt: case 3 sf heatmap
      :align: center
 
  The Sample-Sample Correlation shows clearly that the samples collected at different time frame cluster together. In addition, samples collected closer time points (for instance, 0h and 3h) appears to be more similar. We observe three clusters that show clear differences in open chromatin between the early (cluster 3 - 0h and 3h), intermediate (cluster 2 - 24h), and late stage (cluster 1 - 96h and 120h) time points.
  
  .. note::  Samples in the clustering tree are ordered by what is given in the :ref:`section_metadata`. Simply switch sample order in the Metasheet if you want the Sample-Feature clusters to look more neat. 
  

 3. **Cluster Analysis - Motif and Cistrome Analysis**: 
 
 Following the Sample-Feature heatmap, *CoBRA* is implemented to run a cluster analysis focusing on each cluster of the peaks differentiated by the sample-feature heatmap. 
  
    .. code-block:: Bash

       snakemake cluster_analysis -f
 
 Using the command above, *CoBRA* outputs three additional subfolders in the ``analysis_result/clustering_analysis/rpkm.3_num_sample.2_scale.q_fliter.cov.10`` folder:
  - ``cluster``: includes the peak information in each cluster (bed file and a table containing genes associated with each peak) 
  - ``cistrome_toolkit``: cistrome toolkit analysis giggle plot for each of the cluster 
  - ``motif``: motif analysis result fo reach of the cluster
 
 In the previous part, cluster 3 exhibits to be the peaks differentially upregulated in the 96h and 120h samples. The motifs significantly enriched in these peaks are shown below:
 
 .. figure:: ./tutorial_figures/3_cluster_motif_120.png
      :scale: 20 %
      :alt: case 3 cluster motif
      :align: center


4. **Supervised Analysis - DeSeq2 Differential Peak Analysis**: 

    .. code-block:: Bash

       snakemake run_limma_and_deseq -f
       snakemake run_deeptools_diff_peaks -f
  
  As demonstrated in Case Study 1, these command produces a series of differential peak analysis results located in the ``analysis_result/differential_peaks/{your_comparison}`` folder, including a MA plot and a peak intensity plot. 
  
  .. figure:: ./tutorial_figures/3_maplot.png
      :scale: 35 %
      :alt: case 3 ma plot
      :align: center
  
  .. figure:: ./tutorial_figures/3_peaks.png
      :scale: 30 %
      :alt: case 3 diff peaks
      :align: center
      
  The above MA plot and peak intensity plot are for comparing the 0hr and 120hr samples, and exhibits very robust results. 
  
     
 5. **Pilot Feature - RNA-seq result Intergration**: 
 
  A pilot feature of *CoBRA* that is not implemented in its main snakemake workflow is that it may intergrate differential expression analysis result of the data set's corresponding RNA-seq to create an annotated volcano plot that perfectly illustrated all the differential genes of interest. 
    
    .. code-block:: Bash

       Rscript scripts/volcano_plot.R RNA_seq/120h_over_0h.deseq.csv ChIP_seq/120h_over_0h.deseq.with.Nearby.Gene.csv ref_files/hg19/refGene.hg19.id.bed vol.pdf
       
  
    .. figure:: ./tutorial_figures/3_vol.png
      :scale: 25 %
      :alt: case 3 Volcano Plot
      :align: center

  The command consists of the following inputs:
 
    - ``scripts/volcano_plot.R``: the R script that produces the volcano plot, located in the ``scripts`` folder
    - ``ChIP_seq/120h_over_0h.deseq.with.Nearby.Gene.csv``: DESeq analysis result output by *CoBRA*
    - ``RNA_seq/120h_over_0h.deseq.csv``: the differential expression gene output from RNA-seq result with the same comparison of interest
    - ``ref_files/hg19/refGene.hg19.id.bed``: reference genome file, located in the ``ref_files/{your_genome}`` folder
    - ``vol.pdf``: the pdf file of which the figure will be saved to

  The R script includes the following parameter available for alteration:
 
    - ``Max_Peak_Gene_Distance``: The maximum distance (in bp) considered when pairing a peak with its nearby gene
    - ``Gene_FC_cutoff`` and ``Gene_P_cutoff``: cutoff (fold change and adjusted p-value) for including differentially expressed genes in RNA-seq DE result, default set at 2 and 0.01, respectively
    - ``Peak_FC_cutoff`` and ``Peak_P_cutoff``: cutoff (fold change and adjusted p-value) for including differentially expressed genes in ChIP-seq/ATAC-seq DE result, default set at 2 and 0.01, respectively
    - ``Min_padj``: y_axis limit, any adjusted p-value smaller than this value will be set to this value
    - ``X_axis_limit``: x_axis limit, any fold change value larger than this value will be set to this value
    - ``Genes_To_Label``: the top N genes that will be labeled on this graph
    - ``Transparancy``: transparancy level on genes below the cutoff, default set at 0 (i.e. not shown on the graph)

