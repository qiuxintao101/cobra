.. _workflow:

Pipeline Overview
************************************************************
For more details, please see out pre-print (linked here: :ref:`citation`).

*CoBRA's* conceptual idea and schematic is illustrated here. The first figure introduces the biological motivation. 


   .. figure:: workflow.png
         :scale: 30 %
         :alt: CoBRA schematics
         :align: center

         Conceptual idea and schematic of *CoBRA*, the input and the output

 
 The second figure shows a technical schematic of the *CoBRA* workflow:


   .. figure:: Cobra_workflow.png
      :scale: 16 %
      :alt: Schematic of the CoBRA workflow and the GC binning
      :align: center

      Detailed workflow for bioinformatics data processing

The following is a directed acyclic graph that shows sequence of processes executed by *Snakemake*:
         
   .. figure:: dag.png
         :scale: 20 %
         :alt: Directed acyclic graph of an example workflow
         :align: center
         
         A directed acyclcic graph (DAG). Each *Snakemake* rule is represented by a box, and each dependency is represented by an arrow.


*Snakemake* was used to implement *CoBRA*. For documentation detailing *Snakemake*, see Section :ref:`workingWithPipeline`. As illustrated in the DAG above, *CoBRA* consists of the following rules: 

- ``merged_bed``: Using the bedops package to get the union set of peaks of all samples
- ``bed_enhancer_promoter``:  Filter the union set of peaks within ±1kb of TSS to get enhancer and promoter sites
- ``bedtools_intersect``: Count reads in peak regions for all bam files
- ``pca_plot``: R script that plots all samples in a two dimensional space including PC1 and PC2
- ``heatmapSS_plot``: R script that calculates pair-wise sample-sample correlation and plots all samples with hierarchical clustering
- ``heatmapSF_plot``: R script that performs k-means/hierarchical clustering for the most variable peaks
- ``limma_and_deseq``: R script that performs a differential peak analysis using deseq and limma
- ``deseq_motif``: HOMER performs the known as well as de novo motif enrichment analysis with GC content mathced back ground
- ``GSEA``: GESA analysis on the preranked peak lists. Peaks are assigned to the nearest Gene.
- ``cistrome_tookit``: Calcuate the giggle score to compare the differential peaks to the existing TF ChIP-seq data avaiable on Cistrome DB



Input
************************************************************


Outline
==============================

The following files are needed to run *CoBRA* on your own experiment:

- *Fastq* files with reads for each sample (see :ref:`parameter_FastqFile`)

**OR**

- *BAM* file with aligned reads for each sample (see :ref:`parameter_BamFile`)
- *BED* file with called peaks for each sample (see :ref:`parameter_BedFile`)
- *BIGWIG* file with compressed, indexed, binary format for genome-wide signal data for calculations (see :ref:`parameter_BigwigFile`)

- Optionally: corresponding CNV data (see :ref:`section_cnv`)

Additionally, reference files are all precompiled and are automatically downloaded as needed. These include hg19, hg38,mm9, mm10.

- genome & genome_dict (see :ref:`parameter_RefGenome`)
- refseqGenes (see :ref:`parameter_RefGene`)
- lift chain files (see :ref:`parameter_LiftChain`)
- Cistrome DB in giggle format (see :ref:`parameter_CistromeGiggle`)

Metadata and config files must be filled out by the user to run *CoBRA* on your own experiment:

- a configuration file (:ref:`configurationFile`)
- a metadata file for the samples (:ref:`section_metadata`)


.. _configurationFile:

Configuration file
==============================

A configuration file that defines various parametrs is needed to run *CoBRA*.

.. note:: Please pay attention to the following requirements:

  - Header names should not be changed
  - Absolute and relative paths are acceptable in the config file. When using *Docker*, all input files must be mounted in the container. Please refer to section :ref:`docs-DockerNotes`.
  
All parameters are organized by section. See the following for details:

SECTION ``par_general``
--------------------------------------------

.. _parameter_Project_Name:


``projectName``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default "ChIP_seq". The name will be use for pca, sample-sample, and sample-feature plot titles.

Details
  Please use "_" to seperate different words, as spaces are not allowed.


``enhancer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Enhancer option: enhancer / promoter / all (default). 

Details
  Enhancer options to filter the union set of peaks, which will be used in all analysis in the workflow.


``metasheet``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Location of metasheet, default is metasheet.csv.

Details
  Specifies the location of metasheet that will be used.
  

``ref``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default ""scripts/ref.yaml".

Details
  Specifies the location of ref.yaml that will be used. Most of reference files that will not need to be changed are in the ref.yaml.


``assembly``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default hg19. hg38 / mm9 / mm10 are avaiable.

Details
  Specifies the assembly that the input files are aligned to, all options need to be listed in the ref.yaml.


``rpkm_threshold``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Number. Default 1. This provide a threshold that can be applied to filter the union peak set for all downstream unsupervised analysis.
  
Details
  At least ``mini_num_sample`` should have RPKM > ``rpkm_threshold``


``mini_num_sample``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Number. Default 1. This paramter toghter with rpkm_threshold provide threshold that can apply to filter for all unspuervised analysis.
  
Details
  At least ``mini_num_sample`` should have RPKM > ``rpkm_threshold``


``scale``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default q. The scale method used to nomalize counts for downstream unsupervised analysis.

Details
  The scale method for the normalization options: z- z-score, q- quantile-normalize, l- log-transform


``filter-opt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default cov. Fliter metric in feature selection.

Details
  Metric in feature selection options: sd- Standard deviation, cov- Coefficient of Variation, av- mean


``filter-percent``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >=  0. Default 100. Top percent cutoff that is aplied with ``filter-opt``.

Details
  Top ``filter-percent`` of ``filter-opt`` peaks will be use for the unsupervised analysis.


``SSpeaks``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 20000000. 

Details
  This parameter sets the Maxium number of peaks can be used for the Sample-Sample correlation plot.
  

``SFpeaks``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 20000000. 

Details
  This parameter sets the Maxium number of peaks can be used for the Sample-Feature plot.


``num_kmeans_clust``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 6. 

Details
  This parameter sets the number of clusters that will be used in the k-means clustering for Sample-Feature plot.


``cor_method``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default pearson. Correlation method used for sample-sample and sample-feature plot
  
Details
  The correlation method options: pearson, spearson


``dis_method``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default euclidean. Distance method used for sample-sample and sample-feature plot
  
Details
  Distance measurement options: euclidean, manhattan, canberra, binary, maximum, or minkowski


``Padj``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 0.05. 

Details
  This parameter sets the cut-off for DEseq differential peak calling.


``LG2FC``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >= 0. Default 0. 

Details
  This parameter sets the cut-off for DEseq differential peak calling.

.. _parameter_norm_method:

``norm_method``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
   String. Default depth. DESeq normalization method used for differential expression analysis

Details
  This parameter sets the DESeq normalization method, options: def- normlize by default setting of DEseq2, depth- normlize by the sequence depth of each sample


``motif``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, default 'false'.

Details
  This parameter is use to determine if motif enrichement and clustering analysis is performed.


``bam_sort``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, default 'true'.

Details
  This parameter is needed to flag if the bam files provieded input are sorted or not. If set to 'false', *CoBRA* will automatically sort and reorder the bam files.


``CNV_correction``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, default 'false'.

Details
  This parameter is required to flag if CNV correction should be perfomed or not.


``unchanged_heatmap``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, default 'false'.

Details
  This parameter is required to flag if heatmap change should be perfomed or not.
  

``fastq_in``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, default 'true'.

Details
  This parameter is required to indicate types of file used as input. If `true`, only fastq files for each sample will be used. If `false`, then bed, bam, bigwig will need to be provided


``thread``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >= 0. Default 8. 

Details
  Number of threads used in bwa mem alignment. If run on a local PC, use 1 thread.


SECTION ``samples``
--------------------------------------------

.. _parameter_FastqFile:

``fastq``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the fastq files.

Details
  Path to a fastq file that summarizes the peaks for each sample. The following is an example:
  
  .. code-block:: Bash
  
     bed:
       sample1: ./XX1.fastq
       sample2: ./XX2.fastq

.. _parameter_BedFile:

``bed``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the bed files.

Details
  Path to a bed file that summarizes the called peaks for each sample. The following is an example:
  
  .. code-block:: Bash
  
     bed:
       sample1: ./XX1.bed
       sample2: ./XX2.bed

.. _parameter_BamFile:

``bam``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the bam files.

Details
  Path to a bam file for each sample. The following is an example:
  
  .. code-block:: Bash
  
     bam:
       sample1: ./XX1.bam
       sample2: ./XX2.bam

.. _parameter_BigwigFile:

``bigwig``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the bigwig files.

Details
  Path to a bigwig file for each sample. The following is an example:
  
  .. code-block:: Bash
  
     bigwig:
       sample1: ./XX1.bw
       sample2: ./XX2.bw


.. _section_cnv:

SECTION ``CNV``
--------------------------------------------


``cnv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the igv files for CNV analysis.

Details
  Path to an igv file for each sample. The following is an example:
  
  .. code-block:: Bash
  
     cnv:
       sample1: ./XX1.igv
       sample2: ./XX2.igv
       
  If a file is provided, it must be a valid *igv* file with at least 5 columns:

  - tab-separated columns
  - column names in the first row
  - Columns 1 to 5:

     1. Chromosome
     2. Start
     3. End
     4. Identifier (will be made unique for each if this is not the case already)
     5. log2CNV


SECTION ``additionalInputFiles``
--------------------------------------------


.. _parameter_refGenome_fasta:


``genome``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default hg19.fasta. Path to the reference genome *fasta* file.

Details
  For user convenience, CoBRA will automatic download this file if it has not been downloaded. However, you may also manually create this file to run *CoBRA* on a new species.

  .. Warning:: Chromosome order must correspond to the following files :download:`chr_order.txt <chr_order.txt>` file.. 


``TSS.plus.minus.1kb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path where the refGene plus minus 1kb bed file are stored.

Details
  Each file must be a valid *BED* file with 5 columns, as follows:

  1. chromosome
  2. start
  3. end
  4. strand
  5. Gene_ID

  For user convenience, CoBRA will automatic download this file if it has not been downloaded. However, you may also manually create this file to apply to new species.


``refseqGenes``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Details
  Each file must be a valid *BED* file with 5 columns, as follows:

  1. chromosome
  2. start
  3. end
  4. Gene_ID
  5. Gene_Name

  For user convenience, CoBRA will automatic download this file if it has not been downloaded. However, you may also manually create this file to apply to new species.



``lift.chain``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the lift.chain.gz.

Details
  For user convenience, CoBRA will automatic download this file if it has not been downloaded. This file is used for hg19 and mm9 analysis. It can lift-over coordinates to hg38 and mm10.


``giggle``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the giggle.tar.gz that can be use for cistrome toolkit analysis for finding similar ChIP-seq data that compare to the peaks of interest.

Details
  For user convenience, CoBRA will automatic download this file if it has not been downloaded. It can also be downloaded `here <http://cistrome.org/~chenfei/MAESTRO/giggle.tar.gz>`__.

 
 
.. _section_metadata:


Metadata
=============================================

  
  The metadata file is a comma separated file that contains the annotation and differential comparisson information. The sample names must match those in the configuration file. *CoBRA* can perform as many differential peak analyses as are indicated in the metadata file.
  
  .. warning:: Make sure that end of line characters match default of the operating system. Please convert all line endings to unix format. Please see `here <https://blog.shvetsov.com/2012/04/covert-unix-windows-mac-line-endings.html>`__ .

  Make the  ``metasheet`` file in excel, and save it as a .csv, It doesn’t matter what it is named as long as it is called in the  ``config`` in the section marked  ``metasheet``. See the  ``config`` section for details. The format should be something like the following:

  +--------+------+------------+-----------+------------+--------------------------+
  | Sample | Cell | Condition  | Treatment | Replicates | comp_MCF7_DOX_over_NoDox | 
  +--------+------+------------+-----------+------------+--------------------------+
  | A1     | MCF7 | Full_Media | NoDOX     | 1          | 1                        |
  +--------+------+------------+-----------+------------+--------------------------+
  | A2     | MCF7 | Full_Media | NoDOX     | 2          | 1                        |
  +--------+------+------------+-----------+------------+--------------------------+
  | B1     | MCF7 | Full_Media | DOX       | 1          | 2                        |
  +--------+------+------------+-----------+------------+--------------------------+
  | B2     | MCF7 | Full_Media | DOX       | 2          | 2                        |
  +--------+------+------------+-----------+------------+--------------------------+



  The first column should always contain the sample names that exactly match the sample names used in the config.yaml file.
  The samples that you want to perform a Differential Peak Calling (DE) on using limma and deseq should be marked by the  ``comp`` columns. More on this below.

  .. warning:: This is important! The  ``control`` should be marked with a 1, and the  ``treatment`` should be marked with a 2.

  The remaining metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many ``comp_XXXX`` columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates).

  Again, make this in excel so that all of the spacing is done correctly and save it out as a .csv file. This is the most common bug, so please follow this.
  
  .. warning:: Common Problems with  ``metasheet`` Characters to avoid: ("-", "(", ")", " ", "/", "$"). To avoid bugs, the only punctuation that should be used is the underscore “_”. Dashes, periods, etc, could cause a bug because there is a lot of table formatting and manipulation, or they are invalid characters in R. 
  
  .. note:: CoBRA parses the metadata file and will convert MOST of these invalid characters into '.'--dollarsigns will just be dropped.  The CoBRA parser will also convert between dos/mac files to unix format.
  
  .. note:: It is very important that you know that samples ``A`` is what you mark with 1, and samples ``B`` is what you mark with a 2. You should name your output following this format as well  ``comp_B_over_A`` This will let the reader know what the output DE files refer to. Deseq:  ``baseMeanControl`` refers to samples ``A``, which follows condition 1 and ``baseMeanTreatment`` refers to samples ``B`` which follows condition 2. logfc is ``B/A``

  .. warning:: Do not change the samples data after you started an analysis. You may introduce inconsistencies that will result in error messages. If you need to alter the sample data, we strongly advise you to rerun all steps in the pipeline.


Output
************************************************************

*CoBRA* generates output files that are produced after each of step of the pipeline.

.. note:: Some output folder names are dependent on parameters and comparisons set by the user in the metasheet and config file. Major output filetype and folder structure is described below. 

Common output files can be found in the following formats:

- .bed (bed file)
- .csv (file with comma as column separators)
- .png (PNG format)
- .pdf (PDF format)
- .log (text format)

FOLDER ``Analysis``
=============================================

The final output results are stored here.

Sub-folder ``preprocessed_files``
----------------------------------------------

Stores results related to bam, bed, bigwig, read counts.

.. note:: Output files in this folder do not need to be examined unless itermediate output files are of interest to the user.


Sub-folder ``clustering_analysis``
----------------------------------------------

Stores results related to Principal Component Analysis (PCA) plot, Sample-sample correlation and Sample-Feature clustering plot.


Sub-folder ``differential_peaks``
----------------------------------------------

Stores results related to differential peak calling, motif enrichment, GSEA and cistrome toolkit analyses.


Sub-folder ``logs``
----------------------------------------------

Stores all log files that are created by the pipeline. Each log file is produced by the corresponding rule and contains debugging information as well as warnings and errors.


FOLDER ``preprocessed_files``
=============================================

Stores temporary and intermediate files. Since they are usually not relevant for the user, they are explained in brief.

Sub-folder ``bam``
------------------------------

Stores sorted versions of the *BAMs* that are optimized for fast count.


Sub-folder ``bed``
----------------------------------------------

Stores all original and union bed files, the union peaks are seperated by enhancer and promoter bed files.


Sub-folder ``bigwig``
------------------------------

Stores bigwig files for all samples.


Sub-folder ``read_counts``
------------------------------

Stores sample-peak counts for each sample and merged sample-peak count matrix.


FOLDER ``clustering_analysis``
=============================================


Sub-folder ``rpkm.{}_num_sample.{}_scale.{}_fliter.cov.{}``
------------------------------

Stores unsupervised anlaysis results. Paramaters used for filtering the read counts file is indicated in the folder name.

For example, the folder name 'rpkm.2_num_sample.3_scale.q_fliter.cov.2' means that the unsupervised analysis under this folder is filter by the following criteria:

- ``rpkm.2_num_sample.3`` - at least three samples in the data set have minmal rpkm 2 

- ``scale.q_fliter.cov.2`` - the normalization method is quantile-normalized, fliter metric in feature selection is Coefficient of Variation, the top 2 percent of peaks are being selected.

FILES ``plots/pca_plot.pdf``
----------------------------------------------------------------------------------------------

Details
  Produced in rule ``pca_plot``. PCA is mostly used as a tool in exploratory data analysis. It is often used to visualize distance and relatedness between samples. 

FILES ``plots/heatmapSS_plot.pdf``
----------------------------------------------------------------------------------------------

Details
  Produced in rule ``heatmapSS_plot``. Sample similarity as determined by hierarchical clustering based on the Spearman correlation between samples. 

FILES ``plots/heatmapSF_plot.pdf``
----------------------------------------------------------------------------------------------

Details
  Produced in rule ``heatmapSF_plot``. Peaks from all study samples were merged to create a union set of sites. Each column is a sample, and each row is a peak. K-means clustering is applied to the peak sets. Cluster information can be found in the file "heatmapSF_plot.txt".



FOLDER ``differential_peaks``
=============================================

Sub-folders ``{comparsion_defined_in_metasheet}``
------------------------------

Stores differential anlaysis results that was defined by user in the metasheet. The following are files that can be found in the folder:

- ``{comparsion_defined_in_metasheet}.deseq.csv`` - differential peaks list based on the union peaks. In the file, log2FoldChange and padj for the comparisson of each peak can be found. 

- ``{comparsion_defined_in_metasheet}.deseq.with.Nearby.Gene.csv`` - in addition to the differential peak list, the nearby gene is annotated for each peak.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.up.bed`` - treatment enriched peaks based on the Padj and log2FoldChange cutoff defined in the config file.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.-{}.down.bed`` - control enriched peaks based on the Padj and log2FoldChange cutoff defined in the config file.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.pdf`` - heatmap showing the differential peaks between the treatment and control groups.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.up.bed_motif`` - motif enrichment results for treatment enriched motifs. Both known and de novo results are included.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.-down.bed_motif`` - motif enrichment results for control enriched motifs. Both known and de novo results are included.

- ``GSEA`` - GSEA analysis result based on the log2FoldChange for each nearby gene in the differential peak list.

- ``cistrome_toolkit`` - cistrome_toolkit analysis result based on the treatment and control enriched differential peaks.

- ``DEseq.normalized.counts.csv`` - DEseq normalized counts for each sample and each peak.


FOLDER ``logs``
=============================================

Folder contains log files with errors for each step of the pipeline.

- ``*.log`` A log file is produced for each rule. They contain warnings, errors, and debugging information.

  - ``clean_bam`` logs for picard bam clean 
  - ``remove_duplicates`` logs for picard remove duplicate 
  - ``reorder`` logs for reorder the bam files
  - ``read_counts`` for bedtools intersect to get the sample-peak count matrix

.. _workingWithPipeline:

Running *CoBRA*
******************

General notes
==============================

We present a new pipeline, Containerized workflows for ChIP/ATAC‐seq Experiments (*CoBRA*), that is fast, efficient, portable, customizable and reproducible. The workflow builds upon the ongoing effort to make computational research reproducible using Docker containers. *CoBRA* allows users of varying levels of technical skill to quickly process and analyze new data from ChIP-seq and ATAC-seq experiments. It is the authors’ hope that *CoBRA* can be a starting point for others to build upon and improve *CoBRA* as a tool and extend its ability to analyze the cistrome. 

The *CoBRA* workflow is implemented into a snakemake workflow management system (Köster and Rahmann 2012). Workflows are described via a human-readable, Python-based language. It can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. For ChIP-seq and ATAC-seq experiments, *CoBRA* provides both unsupervised and supervised analyses. 

Further, to make *CoBRA* more easily deployable on any system, it is distributed as a Docker container, which can be used on any machine as long as Docker is installed. Docker containers provide a tool for packaging bioinformatics software. It encapsulates all of the supporting software and libraries, eliminates the possibility of conflicting dependencies, and facilitates the installation of required software. As a result, *CoBRA* is reproducible, portable and easy to deploy.



.. _timeMemoryRequirements:

Running *CoBRA* - Computation time and memory usage
--------------------

*CoBRA* can be computationaly intensive if ``Bam`` files are not sorted. Analyses with a larger sample size (100+ samples) and peak number (10,0000+) generally take longer.


Running time
--------------------

Details about total time consumption:

- the running time is based on the number of samples and the number of peaks.
- for typical analyses in which the sample size is less than 15, running times are roughly 30 minutes with 2 cores for sorted bam files.
- for a large number of samples, running time will be up to 2 hrs or so when executed on a cluster machine.
- if motif analysis is turned on, add 1 additional hour to the running time listed above.



.. _FAQs:

Frequently asked questions (FAQs)
****************************************

The following are commonly asked questions:

1. Why does *CoBRA* need to use a config file and metasheet file to setup the run? Why not just simply use the command to setup the run?

  The unsupervised and supervised anlaysis of ChIP/ATAC-seq experiment requires many paramaters, and could vary from one experiement to another. The config and metasheet files allow the user to save all paramaters that have been used in this run and allow others to reproduce the analysis when needed.

2. Have a problem running docker?

  Please go to https://docs.docker.com/toolbox/faqs/troubleshoot/ to get docker running.

3. How can I rerun a specific part of the pipeline?

  This can be accomplished by running *Snakemake* with the rule name of interest. For example, to produce a new PCA plot or sample-sample heatmap, the following commands can be invoked:
  
     .. code-block:: Bash

        snakemake pca_plot -f
         
        snakemake heatmapSS_plot -f

     ..
   
4. How can I modify the workflow?

  The Snakefile can be modified to change current rules or to accomodate additional ones.


.. _docs-errors:


Troubleshooting
************************************************************

If an issue running *CoBRA* is encountered and you do not find a solution here, please post an issue on our `Bitbucket Issue Tracker <https://bitbucket.org/cfce/cobra/issues>`_ .


Common errors
================

Here are some common errors that users have encountered and reported. 

1. Error in rule ``bedtools_intersect``

  .. code-block:: Bash

    Error in rule bedtools_intersect:

    jobid: 86

    output: ananlysis/preprocessed_files/sample_counts/sample1.total_count,
    ananlysis/preprocessed_files/read_counts/sample_counts/sample1.count

    log: analysis/logs/read_coutns/samle1.log
    RuleException:
         CalledProcessError in line 154 of Snakefile:
  ..

  .. note:: This particular message is normally encountered when the user indicates in the config file that the bam files are sorted when they are not. CoBRA requires that bam files and bed files have the same sorting order. To solve the problem, set the ``bam_sort`` option in the  ``config`` file to ``false``.


2. KeyError in ``metasheet_setup.py``

  .. code-block:: Bash

    *** KeyError in line 9 of Snakefile ***
    File "Snakefile", line 9, in <module>
    File "metasheet_setup.py", line 19, in updateMeta
    File "metasheet_setup.py", line 19, in <dictcomp>
  ..

  .. note:: This particular message appears when a mismatch between the sample names in the ``config`` and ``metasheet`` files exists.


  Simply check if the names are matched to solve this error.


3. rule ``heatmapSS_plot`` duplicate 'row.names' are not allowed

 
  .. code-block:: Bash

     Rscript --default-packages=methods,utils scripts/heatmapSS_plot.R
     analysis/rpkm.1_num_sample.10_scale.q_fliter.cov.100/read_counts/read.counts.rpkm.threshold.scale.fliter.csv
     metasheet.csv 20000000 analysis/rpkm.1_num_sample.10_scale.q_fliter.cov.100/plots/heatmapSS_plot_100_percent.pdf
     analysis/rpkm.1_num_sample.10_scale.q_fliter.cov.100/plots/heatmapSS_100_percent.txt ChIP_seq
     analysis/rpkm.1_num_sample.10_scale.q_fliter.cov.100/plots/images/heatmapSS_plot_100_percent/
     There were 24 warnings (use warnings() to see them)
     Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
     duplicate 'row.names' are not allowed
     Calls: heatmapSS_plot -> read.csv -> read.table
     Execution halted

  This error is normally encountered when you have duplicate sample names in the metasheet.csv. *CoBRA* does not allow duplicate sample names in the ``config`` and ``metasheet`` files.

   

Bug solutions
==============================

When an error is encountered, see the log file that corresponds to the failing *Snakemake* rule. Do a dry run to assess which command must be run. Running the command outside of the workflow will provide a more detailed error message. It is also recommended to check the intermediate files (such as the input and output files of the rule) ensure that they are correct.

Resuming *Snakemake* run
----------------------

After debugging, run *Snakemake* again. It will automatically continue from the rule at which the error occured.


Rerun Incomplete Files 
----------------------

Sometimes when a *Snakemake* run is exited by force (this may also include Docker container exited in the middle of the run, see :ref:`docs-DockerReminder`) and not because of an error, output files would be left incomplete. Then the following error message may appear when trying to resume the *Snakemake* run:

.. code-block:: Bash
   IncompleteFilesException:
   The files below seem to be incomplete. If you are sure that certain files are not incomplete, mark them as complete with

      snakemake --cleanup-metadata <filenames>

   To re-generate the files rerun your command with the --rerun-incomplete flag.
   Incomplete files:
   analysis/preprocessed_files/bam/sorted_reads/453_DMSO_2.bam
   analysis/preprocessed_files/bam/sorted_reads/453_co_NER_1.bam
   analysis/preprocessed_files/bam/sorted_reads/453_co_DMSO_1.bam
   analysis/preprocessed_files/bam/sorted_reads/453_co_NER_2.bam
   analysis/preprocessed_files/bam/sorted_reads/453_NER_1.bam
   analysis/preprocessed_files/bam/sorted_reads/453_co_DMSO_2.bam
..
  
  When this error message appears, simply follow the instruction and add the ``--rerun-incomplete`` flag next to the rule that need to be re-run, for instance:
  
.. code-block:: Bash
   snakemake all --cores 6 --rerun-incomplete
..
  

If you do encounter an error and are unable to find a solution in the FAQ, post an Issue in the `Bitbucket Issue Tracker <https://bitbucket.org/cfce/cobra/issues>`_ tracker.


Customized analysis
****************************************

*CoBRA* is capable of performing unsupervised analyses, differential peak calling, and downstream pathway analysis for ChIP/ATAC‐seq experiemnt. Running *CoBRA* with the default setup is helpful. However, sometimes you may want to further customize the analysis. 


Summary
=================

  - Regardless of whether or not differntial analysis was conducted, we recommend that you first check the pca_plot and heatmapSS_plot pdf.
  - If a specific question is not addressed here, feel free to contact us.
  - *CoBRA* calls differential peaks, which may need different cutoffs for significance for different experiements.



Specifics for the unsupervised analysis
=================================

The following steps are a good starting point:
  - Start with default paramaters
  - Handling the paramaters in the following way:
    - (a) Adjust the ``filter-percent`` to 20 or 100, this will change the percent of the most variable peaks that will go into the unsupervised analysis.
    - (b) Adjust ``num_kmeans_clust`` to change the heatmapSF_plot result for observations in different group of clustering.
    - (c) Quantile-normalize for ``Scale method`` and Coefficient of Variation for ``filter-opt`` is often recommended.


Specifics for the supervised analysis
==========================================

The following steps are a good starting point:
  - Start with the default adjusted p-value threshold (0.05)
  - The following are examples of common occurences when conducting differential analysis:
    - (a) There are very few or 0 significant differential peaks: You may use a less stringent adjusted p-value threshold. You may check the GSEA result even if the there are very few differential peaks, the GSEA analysis provides the enrichment of all nearby genes using the ranking of log2fold change in all peaks. Sometimes the subtle change in the peaks may not reach the significant threshold, but the overal ranking of the peaks may help identify pathways that reflect changes accross the perturbation. 
    - (B) A lot peaks are significant (say more than 10000): You may use a more stringent adjusted p-value and log2fold change threshold. Check the deeptools heatmap to see if called differential peaks are truly differential. 


