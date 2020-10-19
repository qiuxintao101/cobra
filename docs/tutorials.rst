
.. _docs-quickstart:

********
Case Study
********

The three case studies below gives you a taste of different capabilities of our CoBRA workflow. The data sets you will be exploring are, 
1) a glucocorticoid receptor (GR) ChIP-seq data set from the ENCODE project, 
2) a set of H3K27ac ChIP-seq data from colon cancer cell lines, and
3) an ATAC-seq experiment on HL-60 promyelocytes differentiating into macrophages. 

Each of the case studies wiil demonstrate a subfield of expertise of our CoBRA pipeline. 


Setup
=====

In preparation for the tutorials, please use the following steps to set up the cobra environment and retrieve the latest version of our pipeline:

1. **Initiate Environment**: 
  
  Use the following command to activate the cobra environment:
  
  .. code-block:: Bash

    source activate cobra

2. **Retrieve the Latest Version of Cobra:**

  If installed using docker, run the following command to change the working directory. Otherwise, skip to next command:
   
  .. code-block:: Bash
   
     cd cobra
   
  Create a new directory for the tutorial (ex: gr_chip), change to the new working directory and pull the latest version of CoBRA using ``git clone`` :

  .. code-block:: Bash

     mkdir gr_chip
     git clone https://bitbucket.org/cfce/cobra.git .

  If you receive an error, *Git* may not be installed on your system. Please consult the internet on how to best install Git for your system.


Case Study 1: GR ChIP-set Data Set
================

Background
**********
This tutorial makes use of a publicly available glucocorticoid receptor (GR) ChIP-seq data (`GSE32465 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32465>`_) from a lung adenocarcinoma cell line (A549) at 3 different concentrations of dexamethasone, a potent GR agonist. In the experiment, samples were treated with 0.5nM, 5nM, or 50nM of dexamethasone. 


Download and set-up for running the GR_ChIP sample dataset
**********************************************************

  Please use the following command to download the GR_ChIP sample dataset. This dataset is of moderate size (3.9 G) and may take 5-10 minutes to download. 

  .. code-block:: Bash
   
     snakemake download_example_GR_ChIP
  
  When the data set is downloaded, we can proceed to set up for the run. Usually for running CoBRA on a new experiment, the two config files ``config.yaml`` and ``metasheet.csv`` would need to be set up acccordingly. In this tutorial, they have been filled already. Note in ``config.yaml``, the parameter `motif` has been set as true to perform motif enrichement and clustering analysis.

  To check if the setup is correct, begin a dry run via the following command:

  .. code-block:: Bash

     snakemake all -np


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
      :scale: 25%
      :alt: case 1 pca plot
      :align: center

  .. figure:: ./tutorial_figures/1_pca_scree.png
      :scale: 25 %
      :alt: case 1 pca scree
      :align: center

  .. figure:: ./tutorial_figures/1_pca.png
      :scale: 30 %
      :alt: tutorial 1 pca plot
      :align: center
      
  As illustrated in the PCA plot, PC1 separates the samples with different treatment concentration of dexamethasone, while PC2 further    separates the sample replicates.
 
  .. figure:: ./tutorial_figures/1_pca_scree.png
      :scale: 30 %
      :alt: tutorial 1 pca scree
      :align: center

  As illustrated in the PCA plot and scree plot above, PC1 (capturing 40.8% of variance explained) separates the samples with different treatment concentration of dexamethasone - namely 0.5nM from 5nM and 50nM, while PC2 (18.7% variance) further separates the sample replicates.


2. **Unsupervised Analysis - Sample-Sample Correlation Plot**: 

    .. code-block:: Bash

       snakemake heatmapSS_plot -f
  
  This command produces the ``heatmapSS_plot_100_percent.pdf`` file located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. It provides information on the clustering result based on the Pearson correlation coefficient, and illustrates the similarity between all samples in a pairwise fashion.
  
  .. figure:: ./tutorial_figures/1_SS.png
      :scale: 28 %
      :alt: case 1 ss heatmap
      :align: center
      
  As illustrated in the Sample-Sample correlation plot, samples replicates cluster tightly together (r > 0.6). And samples treated with 0.5nM of dexamethasone exhibited to be far different from samples treated with 5nM or 50nM dexamethasone.


3. **Unsupervised Analysis - Sample-Feature Heatmap**: 

    .. code-block:: Bash

       snakemake heatmapSF_plot -f
  
  This command produces the ``heatmapSF_plot_100_percent.pdf`` file located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. It illustrates clustering of samples based on correlation on the horizontal axis and clustering of peaks on the vertical axis. It presents patterns of peaks (by k-means clustering) across samples and identifies the clusters that are enriched in a subset of samples.
  
  .. figure:: ./tutorial_figures/1_SF.png
      :scale: 28 %
      :alt: case 1 sf heatmap
      :align: center
 

4. **Supervised Analysis - Limma/DeSeq2 Differential Peak Analysis**: 

  The key inquiry to be satisfied for any ChIP-seq/ATAC-seq analysis is what the differential sites are between sample groups of interest. In *CoBRA*, this analysis is done by incorporating differential peak callin gby DESeq2 while using sequencing depth as a scale factor, and thus significantly reducing false positive differential peak-calling.
  
    .. code-block:: Bash

       snakemake limma_and_deseq -f
  
  This command produces a series of files located in the ``analysis_result/differential_peaks/c50nm_vs_0.5nm`` folder, including the following:
  - ``c50nm_vs_0.5nm.deseq.csv``: a differentail peaks analysis table produced by DESeq2
  - ``c50nm_vs_0.5nm.deseq.Padj0.05.LG2FC.0.up.bed`` and ``c50nm_vs_0.5nm.deseq.Padj0.05.LG2FC.-0.down.bed``: bed files of peaks that are differentially up- and down-regulated, respectively
  - ``c50nm_vs_0.5nm.deseq.sum.csv``: a table including total number of differential peaks under different thresholds
  - ``c50nm_vs_0.5nm.t.test.csv``: a t-test table of the differential peaks
  - ``MA_plot.pdf``: a MA plot comparing the two treatment samples
  
  .. figure:: ./tutorial_figures/1_maplot.png
      :scale: 50 %
      :alt: case 1 ma plot
      :align: center
  
  The MA plot above shows that the 50nM treatment samples have significant numbers of upregulated peaks called by DESeq2 and no downregulated peaks.
  
  Intensity measurement of the differnetial peaks can be done using the following command
  
    .. code-block:: Bash

       snakemake deeptools_diff_peaks -f
  
  It produces ``c50nm_vs_0.5nm.deseq.Padj0.05.LG2FC.0.pdf`` which illustrates the peak intensity of the differentially up and downregulated peaks. 

  .. figure:: ./tutorial_figures/1_peaks.png
      :scale: 50 %
      :alt: case 1 diff peats
      :align: center
       
  The peak-intensity heatmap above further illustrates taht there only exist differentially upregulated peaks in 50nM treatment samples as compared to 0.5 nM dexamethasone treated samples, and intensity goes as high as 1.75.


5. **Comparison of Up and Down-regulated Site: Cistrome Toolkit**: 

  *CoBRA* has a built-in feature that compares up and down-regulated sites to a comprehesnive database of ChIP/ATAC and DNase data, and outline a series of most similar samples in terms of genomic interval overlaps with the differential sites located in the Cistrome database. This feature allows researchers to pin-point those similar data set of interest and download for further investigation. It can provide unique insight into gained or lost sites such as identifying which transcription factor potentially binds to a differential peak set after a perturbation and in investigating similar cellular systems.
  
    .. code-block:: Bash

       snakemake cistrome_tookit -f
  
  Using the command above, *CoBRA* outputs a series of files located in the ``analysis_result/differential_peaks/c50nm_vs_0.5nm/cistrome_toolkit`` folder, including:
  - a plot of most similar samples ranked by their giggle score, and
  - two tables of cistrome toolkit result, each include a list of GEO accession numbers corresponding to all ChIP-seq data with similarity to the differential peak set (up or down-regulated)
  
  .. figure:: ./tutorial_figures/1_cistrome.png
      :scale: 40 %
      :alt: case 1 cistrome result
      :align: center

  As show in the plot above, for the gained GR binding sites in the dexamethasone treatment, the NR3C1 factor in Lung is the most similar ChIP-seq in the Cistrome database to this GR data set.


Case Study 2: H3K27ac ChIP-seq Data Set
================

Background
**********
This tutorial makes use six samples from several experiments: three Microsatellite Instable (MSI) samples and three Microsatellite Stable (MSS) samples (Tak et al. 2016; Piunti et al. 2017; Piunti et al. 2017; Maurano et al. 2015; McCleland et al. 2016; Rahnamoun et al. 2018). Microsatellite Instable (MSI) and Microsatellite Stable (MSS) are two classses used to characterize colorectal cancers. MSS tumors are one of the most highly mutated tumor types (Taieb et al. 2017) and exhibit a high number of copy number variations. Without adjustment, a differential peak caller will rank peak loci with high copy number gain in MSS as being the most differential compared to MSI. To observe differential peaks between the MSI and MSS samples, *CoBRA* allows for **copy number variation adjustment** during the supervised analysis.


Download and set-up for running the H3K27ac ChIP-seq sample dataset
**********************************************************

  Please use the following command to download the H3K27ac ChIP-seq sample dataset. 

  .. code-block:: Bash
   
     snakemake download_example_H3K27ac_ChIP
  
  When the data set is downloaded, we can proceed to set up for the run. Note in ``config.yaml``, the parameter `cnv` has laid out a path for **CNV bam files*** corresponding to each sample.

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
       snakemake heatmapSF_plot -f
  
  As demonstrated in the previous case study, these command produces the pca plot and the heatmaps located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. 

  .. figure:: ./tutorial_figures/2_pca.png
      :scale: 28 %
      :alt: tutorial 2 pca plot
      :align: center
      
  .. figure:: ./tutorial_figures/2_pca_scree.png
      :scale: 28 %
      :alt: tutorial 2 pca scree
      :align: center

  As illustrated in the PCA plot and scree plot above, PC1 (capturing 44.5% of variance explained) clearly separates the MSS samples (colored in turquois) and MSI samples (colored in pink).

  
  .. figure:: ./tutorial_figures/2_SS.png
      :scale: 28 %
      :alt: tutorial 2 ss heatmap
      :align: center

  .. figure:: ./tutorial_figures/2_SF.png
      :scale: 28 %
      :alt: tutorial 2 sf heatmap
      :align: center
 
  The Sample-Sample Correlation shows clearly that the MSS samples cluster together, and the same applies to the MSI samples. And the two sample groups exhibit little correlation. 


2. **Supervised Analysis - Limma/DeSeq2 Differential Peak Analysis**: 

    .. code-block:: Bash

       snakemake limma_and_deseq -f
       snakemake deeptools_diff_peaks -f
  
  As demonstrated in Case Study 1, these command produces a series of differential peak analysis results located in the ``analysis_result/differential_peaks/MSS_vs_MSI`` folder, including a MA plot and a peak intensity plot. Applying copy number variation adjustment eliminates false positive peaks that would otherwise be called as differential due to their significant copy number difference between the two sample groups MSI and MSS.
  
  
  .. figure:: ./tutorial_figures/2_maplot.png
      :scale: 50 %
      :alt: tutorial 2 ma plot
      :align: center
      
      MA Plot with CNV Adjustment
  
  .. figure:: ./tutorial_figures/2_maplot_nocnv.png
      :scale: 50 %
      :alt: tutorial 2 ma plot no cnv
      :align: center
      
      MA Plot with No CNV Adjustment
  
  Comparing the two MA Plots above, differential peaks in the MA Plot generated with CNV adjustment exhibits less significant log fold change. 
  
  .. figure:: ./tutorial_figures/2_peaks.png
      :scale: 50 %
      :alt: tutorial 2 diff peaks
      :align: center
      
      Peaks Intensity Plot with CNV Adjustment
     
  .. figure:: ./tutorial_figures/2_peaks_nocnv.png
      :scale: 50 %
      :alt: tutorial 2 diff peaks no cnv
      :align: center
      
      Peaks Intensity Plot with No CNV Adjustment
  
  Comparing the two peak intensity heatmaps above, differential peaks in the plot generated with CNV adjustment generally shows in general higher intensity.
  
  
3. **GSEA**: 

  *CoBRA* has built-in features to do the Gene Set Enrichment analysis, which is performed on the ranked list of genes produced by the pipeline.
  
    .. code-block:: Bash

       snakemake GSEA -f
  
  Using the command above, *CoBRA* outputs a series of GSEA analysis results in ``analysis_result/differential_peaks/MSS_vs_MSI/GSEA`` folder, including:
    - ``gsea_report_for_na_neg`` and ``gsea_report_for_na_pos``: summary report including all ranked genes sets and their statistics 
    - ``neg_snapshot.html`` and ``pos_snapshot.html``: snapshots of all enrichment plots of enriched gene sets curated
    - ``enplot_{Gene_Set}``: individual enrichment plots of an enriched gene set
    - ``{Gene_Set}.html`` and ``{Gene_Set}.xls``: individual GSEA Results Summary of an enriched gene set
  
  .. figure:: ./tutorial_figures/2_gsea_farmer1.png
      :scale: 50 %
      :alt: tutorial 2 GSEA
      :align: center
      
      An Example Enrichment Plot
  
  Without CNV adjustment, GSEA will indicate greatest enrichment in gene sets solely related to amplification. As a result, it is challenging to assess the true epigenetic differences between the two colorectal cancer types. For instance, the gene set NIKOLSKY_BREAST_CANCER_8Q12_Q22_AMPLICON includes genes up-regulated in non-metastatic breast cancer tumors with amplification in the 8q22 region. Without adjustment for copy number variation, this gene set is significantly enriched and ranked 6th, with a normalized enrichment score of -1.91 and an adjusted p-value less than 0.0001. With CNV adjustment, this gene set is far less enriched and ranked 55th, and has a normalized enrichment score of -1.69 and an adjusted p-value of 0.076.


Case Study 3: ATAC-seq from HL-60 promyelocytes differentiating into macrophages
================

Background
**********
This tutorial makes use of ATAC-seq from HL-60 promyelocytes differentiating into macrophages. The samples were taken utilized a five-day time course (0hr, 3hr, 24hr, 96hr, and 120hr) to profile accessible chromatin of HL-60 promyelocytes ddifferentiating into macrophages. Here *CoBRA* results shows investigation of the differentiation of macrophages through changes
in the landscape of accessible chromatin. 

Download and set-up for running the Macrophage_atac sample dataset
**********************************************************

  Please use the following command to download the Macrophage ATAC-seq sample dataset. 

  .. code-block:: Bash
   
     snakemake download_example_Macrophage_atac
  
  When the data set is downloaded, we can proceed to set up for the run. 


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

  While the CoBRA pipeline is designed to be fast and efficient, easily-excuetable with just a few lines of commands, it is possible to produce the analysis in a step-wise fashion by running specific parts of the pipeline.

1. **Unsupervised Analysis - PCA Plot, Sample-Sample Correlation Plot, Sample-Feature Heatmap, etc.**: 

    .. code-block:: Bash

       snakemake pca_plot -f
       snakemake heatmapSS_plot -f
       snakemake heatmapSF_plot -f
  
  Like illustrated in Case Study 1, this command produces the pca plot and the heatmaps located in the ``analysis_result/clustering_analysis/rpkm.1_num_sample.0_scale.q_fliter.cov.100/plots`` folder. 

  .. figure:: ./tutorial_figures/3_pca.png
      :scale: 28 %
      :alt: tutorial 3 pca plot
      :align: center
      
  .. figure:: ./tutorial_figures/3_pca_scree.png
      :scale: 28 %
      :alt: tutorial 3 pca scree
      :align: center

  As illustrated in the PCA plot and scree plot above, PC1 (capturing 57=0.7% of variance explained) clearly separates the samplesby their time frame
  
  .. figure:: ./tutorial_figures/3_SS.png
      :scale: 28 %
      :alt: tutorial 3 ss heatmap
      :align: center

  .. figure:: ./tutorial_figures/3_SF.png
      :scale: 28 %
      :alt: tutorial 3 sf heatmap
      :align: center
 
  The Sample-Sample Correlation shows clearly that the samples collected at different time frame cluster together. In addition, samples collected closer time points (for instance, 0h and 3h) appears to be more similar. We observe three clusters that show clear differences in open chromatin between the early (cluster 1), intermediate (cluster 2), and late stage (cluster 3) time points.


2. **Supervised Analysis - Limma/DeSeq2 Differential Peak Analysis**: 

    .. code-block:: Bash

       snakemake limma_and_deseq -f
       snakemake deeptools_diff_peaks -f
  
  As demonstrated in Case Study 1, these command produces a series of differential peak analysis results located in the ``analysis_result/differential_peaks/{your_comparison}`` folder, including a MA plot and a peak intensity plot. 
  
  .. figure:: ./tutorial_figures/3_maplot.png
      :scale: 50 %
      :alt: tutorial 2 ma plot
      :align: center
      
  
  .. figure:: ./tutorial_figures/3_peaks.png
      :scale: 50 %
      :alt: tutorial 2 diff peaks
      :align: center
     
 3. **Supervised Analysis - Limma/DeSeq2 Differential Peak Analysis**: 
 
    .. code-block:: Bash

       snakemake deseq_motif -f
       
  The command above runs de novo motif analysis on each cluster of accessible sites across all 3 clusters automatically to identify potential transcriptional regulators enriched in differentially accessible chromatin elements. The results are located in the ``analysis_result/differential_peaks/{your_comparison}``, including two different subfolder ``analysis_result/differential_peaks/{your_comparison}/{your_comparison}.{thresholds}.up.bed_motifs``, and ``analysis_result/differential_peaks/{your_comparison}/{your_comparison}.{thresholds}.down.bed_motifs``.
  

 
