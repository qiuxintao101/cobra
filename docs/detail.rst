.. _workflow:

Workflow
************************************************************
For full information, please see the latest publication (linked here: :ref:`citation`).


The workflow and conceptual idea behind *CoBRA* is illustrated by the following three Figures. First, we give a high-level conceptual overview and a biological motivation:

   .. figure:: workflow.png
         :scale: 30 %
         :alt: CoBRA schematics
         :align: center

         Conceptual idea and workflow of *CoBRA*, the input and the output


 Next, we show a schematic of the *CoBRA* workflow from a more technical perspective by showing the actual steps that are performed:


   .. figure:: Cobra_workflow.png
      :scale: 16 %
      :alt: Schematic of the CoBRA workflow and the GC binning
      :align: center

      Summary workflows for data processing and methodological details for *CoBRA* (for more details, see Suppl. Figure 1 in the publication).


We now show which rules are executed by *Snakemake* for a specific example (see the caption of the image):
         
   .. figure:: dag.png
         :scale: 20 %
         :alt: Directed acyclic graph of an example workflow
         :align: center
         
         Exact workflow (a so-called directed acyclic graph, or DAG) that is executed. Each node represents a rule name as defined in the Snakefile, and each arrow a dependency.

*CoBRA* is currently implemented as a *Snakemake* pipeline. For a gentle introduction about *Snakemake*, see Section :ref:`workingWithPipeline`. As you can see, the workflow consists of the following steps or *rules*:

- ``merge_bed``: Using the bedops to get the union set of peaks that is exist in sample sets
- ``bed_enhancer_promoter``:  Filter the union peaks with ±1kb of TSS to get enhancer and promoter sites
- ``bedtools_intersect``: Count all reads for peak regions across all bam files
- ``pca_plot``: R script that plot all samples in a two dimensional space with PC1 and PC2
- ``heatmapSS_plot``: R script that calculate pair-wise sample-sample correlation and plot all samples with hierarchical clustering
- ``heatmapSF_plot``: R script that performs k-means/hierarchical clustering for the most variable peaks
- ``limma_and_deseq``: R script that performs a differential accessibility analysis for the peak regions based on deseq and limma
- ``deseq_motif``: HOMER to performs the know as well as de novo motif enrichment analysis with GC content mathced back ground
- ``GSEA``: GESA analysis on the preranked peak lists, gene are assigned to the nearist peaks
- ``cistrome_tookit``: Calcuate the giggle socre compare the differential peaks to the existing TF ChIP-seq data avaiable on Cistrome DB



Input
************************************************************


Summary
==============================

As input for *CoBRA* for your own analysis, the following data are needed:

- *BAM* file with aligned reads for each sample (see :ref:`parameter_BamFile`)
- *BED* file with called peaks for each sample (see :ref:`parameter_BedFile`)
- *BIGWIG* file with compressed, indexed, binary format for genome-wide signal data for calculations (see :ref:`parameter_BigwigFile`)
- Optionally: corresponding CNV data (see :ref:`parameter_CNVFile`)

In addition, the following files are need, all of which we provide already for human hg19, hg38 and mouse mm9, mm10:

- genome & genome_dict (see :ref:`parameter_RefGenome`)
- refseqGenes (see :ref:`parameter_RefGene`)
- lift chain files (see :ref:`parameter_LiftChain`)
- Cistrome DB in giggle format (see :ref:`parameter_CistromeGiggle`)


Lastly, some metadata files are needed that specify *CoBRA*-specific and Snakemake-specific parameters. They are explained in detail in the next sections. If this sounds complicated, don't worry, just take the example analysis, and you will understand within a few minutes what these files are:

- a general configuration file (:ref:`configurationFile`)
- a metadata file for the samples (:ref:`section_metadata`)


.. _configurationFile:

General configuration file
==============================

To run the pipeline, a configuration file that defines various parameters of the pipeline is required.

.. note:: Please note the following important points:

  - the name of this file is irrelevant, but it must be in the right format (JSON) and it must be referenced correctly when calling *Snakemake* (via the ``--configfile`` parameter). We recommend naming it ``config.yaml``
  - neither section nor parameter names must be changed.
  - For parameters that specify a path, both absolute and relative paths are possible.  We recommend specifying an absolute path. Relative paths must be specified relative to the *Snakemake* working directory.
  - For parameters that specify a directory, there should be no trailing slash.

In the following, we explain all parameters in detail, organized by section names.


SECTION ``par_general``
--------------------------------------------

.. _parameter_Project_Name:


``projectName``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default "ChIP_seq". The name will be use for pca, sample-sample, sample-feature plot titles.

Details
  Please use "_" to seperate different words.



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
  Specifies the location of ref.yaml that will be used. Most of reference files that will not need to be changed commonly are in the ref.yaml.


``assembly``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default hg19. hg38 / mm9 / mm10 are avaiable.

Details
  Specifies the assembly that the input files are aligned to, all options need to be listed in the ref.yaml.


``rpkm_threshold``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Number. Default 1. This provide threshold that can apply to filter for all unspuervised analysis.
  
Details
  At least mini_num_sample should have RPKM > rpkm_threshold


``mini_num_sample``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Number. Default 1. This paramter toghter with rpkm_threshold provide threshold that can apply to filter for all unspuervised analysis.
  
Details
  At least mini_num_sample should have RPKM > rpkm_threshold


``scale``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default q. The scale method for the nomalize counts among samples.

Details
  The scale method for the normaliztion: z- z-score, q- quantile-normalize, l- log-transform


``filter-opt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Default cov. Fliter metric in feature selection.

Details
  Metric in feature selection: sd- Standard deviation, cov- Coefficient of Variation, av- mean


``filter-percent``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer >=  0. Default 100. Top percent cutoff that will apply with the filter metric.

Details
  Top filter-percent of filter-opt peaks will be use for the un-supervised analysis.


``SSpeaks``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 20000000. 

Details
  This parameter sets the Maxium peaks can be used for Sample-Sample correlation plot.
  

``SFpeaks``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 20000000. 

Details
  This parameter sets the Maxium peaks can be used for Sample-Feature plot.


``num_kmeans_clust``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  Integer > 0. Default 6. 

Details
  This parameter sets the number of clusters that will be used in the k-means clustering.


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


``motif``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, deafult 'false'.

Details
  This parameter is use to decide the on and off for the motif enrichement and clustering analysis.


``bam_sort``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, deafult 'true'.

Details
  This parameter is to flag if the bam files provied for input are sorted or not.


``CNV_correction``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String, deafult 'false'.

Details
  This parameter is to flag if the CNV correction should be perfomed or not.



SECTION ``samples``
--------------------------------------------

.. _parameter_summaryFile:


``bed``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the bed files.

Details
  Path to a bed file that summarizes the peak information for the data. Following is an example:
  
  .. code-block:: Bash
  
     bed:
       sample1: ./XX1.bed
       sample2: ./XX2.bed


``samples``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the bam files.

Details
  Path to a bam file for each sample. Following is an example:
  
  .. code-block:: Bash
  
     bam:
       sample1: ./XX1.bam
       sample2: ./XX2.bam


``bigwig``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the bigwig files.

Details
  Path to a bigwig file for each sample. Following is an example:
  
  .. code-block:: Bash
  
     bigwig:
       sample1: ./XX1.bw
       sample2: ./XX2.bw




SECTION ``CNV``
--------------------------------------------


``cnv``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Summary
  Paths to the igv files for CNV analysis.

Details
  Path to a igv file for each sample. Following is an example:
  
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
  For user convenience, CoBRA will automatic download this file if it is not been downloaded. However, you may also manually create this file to apply to new species.
  .. warning:: You need write access to the directory in which the *fasta* file is stored, make sure this is the case or copy the *fasta* file to a different directory. The reason is that the pipeline produces a *fasta* index file, which is put in the same directory as the corresponding *fasta* file. This is a limitation of *samtools faidx* and not our pipeline.

  .. note:: This file has to be in concordance with the input data; that is, the exact same genome assembly version must be used. In the first step of the pipeline, this is checked explicitly, and any mismatches will result in an error.


``TSS.plus.minus.1kb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the refGene plus minus 1kb bed file are stored.

Details
  Each file must be a valid *BED* file with 5 columns, as follows:

  1. chromosome
  2. start
  3. end
  4. strand
  5. Gene_ID

  For user convenience, CoBRA will automatic download this file if it is not been downloaded. However, you may also manually create this file to apply to new species.


``refseqGenes``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Details
  Each file must be a valid *BED* file with 5 columns, as follows:

  1. chromosome
  2. start
  3. end
  4. Gene_ID
  5. Gene_Name

  For user convenience, CoBRA will automatic download this file if it is not been downloaded. However, you may also manually create this file to apply to new species.



``lift.chain``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the lift.chain.gz.

Details
  For user convenience, CoBRA will automatic download this file if it is not been downloaded.


``giggle``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Summary
  String. Path to the giggle.tar.gz, that can be use for cistrome toolkit analysis for finding similar ChIP-seq data that compare to the peaks of interested.

Details
  For user convenience, CoBRA will automatic download this file if it is not been downloaded. Can also be downloaded `here <http://cistrome.org/~chenfei/MAESTRO/giggle.tar.gz>`__.

 
 
.. _section_metadata:


Input metadata
=============================================

  This file summarizes the data and corresponding available metadata  that should be used for the analysis. The format is flexible and may contain additional columns that are ignored by the pipeline, so it can be used to capture all available information in a single place. Importantly, the file must be saved as comma-separated, the exact name does not matter as long as it is correctly specified in the configuration file.

  .. warning:: Make sure that the line endings are correct. Different operating systems use different characters to mark the end of line, and the line ending character must be compatible with the operating system in which you run *CoBRA*. For example, if you created the file in MAC, but you run it in a Linux environment (e.g., a cluster system), you may have to convert line endings to make them compatible with Linux. For more information, see `here <https://blog.shvetsov.com/2012/04/covert-unix-windows-mac-line-endings.html>`__ .

  Make the  ``metasheet`` file in excel, and save it as a .csv, It doesn’t matter what it is named as long as it is called in the  ``config`` in the spot marked  ``metasheet`` see the  ``config`` section if confused. The format should be something like the following:

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



  The first column should always be sample names that exactly match the sample names used in config.yaml
  The samples that you want to perform a Differential Peak Calling (DE) on using limma and deseq should be marked by the  ``comp`` columns more on this below

  .. warning:: This is important! The  ``control`` should be marked with a 1, and the  ``treatment`` should be marked with a 2.

  It is recommended that if you should have a “replicates” column to denote different samples, it is a good idea to not only have each of the sample names be unique, but also make sure that the associated metadata is unique as well to each sample.
  The rest of the  metadata columns are up to the user to write. Sample must always be first, and you are allowed to have as many ``comp_XXXX`` columns as you want at the end. All of the middle columns are your metadata (for this example, this is cell, condition, treatment, replicates)

  Again, make this in excel so that all of the spacing is done correctly and save it out as a .txt or .csv file. This is the most common bug, so please follow this.
  
  .. warning:: Common Problems with  ``metasheet`` Characters to avoid: ("-", "(", ")", " ", "/", "$") To avoid bugs, the only punctuation that should be used is the underscore “_”. Dashes, periods, etc, could cause a bug because there is a lot of table formatting and manipulation, or they are invalid characters in R. 
  
  .. note:: CoBRA parses the meta file and will convert MOST of these invalid characters into '.'--dollarsigns will just be dropped.  The CoBRA parser will also convert between dos/mac files to unix format.
  
  .. note:: It is very important that you know that samples ``A`` is what you mark with 1, and samples ``B`` is what you mark with a 2. You should name your output following this format as well  ``comp_B_over_A`` This will let the reader know what the output DE files refer to. Deseq:  ``baseMeanA`` refers to samples ``A``, which follows condition 1 and ``baseMeanB`` refers to samples ``B`` which follows condition 2. logfc is ``B/A``

  .. warning:: Do not change the samples data after you started an analysis. You may introduce inconsistencies that will result in error messages. If you need to alter the sample data, we strongly advise to recalculate all steps in the pipeline.


Output
************************************************************

The pipeline produces quite a large number of output files, only some of which are however relevant for the regular user.

.. note:: In the following, the directory structure and the files are briefly outlined. As some directory or file names depend on specific parameters in the configuration file, curly brackets will be used to denote that the filename depends on a particular parameter or name. For example, ``{comparisonType}`` and ``{regionExtension}`` refer to ``comparisonType`` (:ref:`parameter_comparisonType`) and ``regionExtension`` ( :ref:`parameter_regionExtension`) as specified in the configuration file.

Most files have one of the following file formats:

- .bed (bed file)
- .csv (file with comma as column separators)
- .png (PNG format)
- .pdf (PDF format)
- .log (text format)

FOLDER ``Analysis``
=============================================

In this folder, the final output files are stored. Most users want to examine the files in here for further analysis.

Sub-folder ``preprocessed_files``
----------------------------------------------

Stores results related to bam, bed, bigwig, read counts.

.. note:: Output files in this folder do not need to examine unless itermidate output files are interested to the user.


Sub-folder ``clustering_analysis``
----------------------------------------------

Stores results related to Principal Component Analysis (PCA) plot, Sample-sample correlation and Sample-Feature clustering plot.


Sub-folder ``differential_peaks``
----------------------------------------------

Stores results related to differential peaks calling, motif, GSEA and cistrome toolkit analyses.


Sub-folder ``logs``
----------------------------------------------

Stores results related to all logs that is created by the pipeline .Each log file is produced by the corresponding rule and contains debugging information as well as warnings and errors:


FOLDER ``preprocessed_files``
=============================================

Stores temporary and intermediate files. Since they are usually not relevant for the user, they are explained only very briefly here.

Sub-folder ``bam``
------------------------------

Stores sorted versions of the *BAMs* that are optimized for fast count.


Sub-folder ``bed``
----------------------------------------------

Stores all original and union bed files, the union peaks are seperate by enahcner and promoter bed files.


Sub-folder ``bigwig``
------------------------------

Stores all bigwig files for all samples.


Sub-folder ``read_counts``
------------------------------

Stores sample-peak count for each sample and merged sample-peak count matrix.


FOLDER ``clustering_analysis``
=============================================


Sub-folder ``rpkm.{}_num_sample.{}_scale.{}_fliter.cov.{}``
------------------------------

Stores unsupervised anlaysis results that paramaters used for filter the read counts is indicated in the folder name.

For example, the folder name 'rpkm.2_num_sample.3_scale.q_fliter.cov.2' means that the unsupervised analysis under this folder is filter by the following criteria.

- ``rpkm.2_num_sample.3`` - at least three samples in the data set have minmal rpkm 2 

- ``scale.q_fliter.cov.2`` - the normalization method is quantile-normalized, fliter metric in feature selection is Coefficient of Variation, the top 2 percent of peaks are being selected.

FILES ``plots/pca_plot.pdf``
----------------------------------------------------------------------------------------------

Details
  Produced in rule ``pca_plot``. PCA is mostly used as a tool in exploratory data analysis. It is often used to visualize  distance and relatedness between samples. 

FILES ``plots/heatmapSS_plot.pdf``
----------------------------------------------------------------------------------------------

Details
  Produced in rule ``heatmapSS_plot``. Sample similarity was determined by hierarchical clustering using the Spearman correlation between samples. 

FILES ``plots/heatmapSF_plot.pdf``
----------------------------------------------------------------------------------------------

Details
  Produced in rule ``heatmapSF_plot``. Peaks from all study samples were merged to create a union set of sites. Each column is a sample, and each row is a peak. K-means clustering were applied to the peak sets, cluster information can be locate in the file "heatmapSF_plot.txt".



FOLDER ``differential_peaks``
=============================================

Sub-folders ``{comparsion_defined_in_metasheet}``
------------------------------

Stores differential anlaysis results that was defined by user in the metasheet.csv. Following are files inside the folder:

- ``{comparsion_defined_in_metasheet}.deseq.csv`` - differential peaks list based on the union peaks, log2FoldChange, padj for each peak is in the list. 

- ``{comparsion_defined_in_metasheet}.deseq.with.Nearby.Gene.csv`` - in addition to the differential peak list, the nearby gene is annotate for each peak.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.up.bed`` - treatment enriched peaks based on the Padj and log2FoldChange cut off defined in the config.yaml.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.-{}.down.bed`` - control enriched peaks based on the Padj and log2FoldChange cut off defined in the config.yaml.

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.pdf`` - heatmap that showing the differential peaks between the treatment and control group

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.up.bed_motif`` - motif enrichedment results for treatment enriched motifs, both know and denovo results are included

- ``{comparsion_defined_in_metasheet}.deseq.Padj{}.LG2FC.{}.-down.bed_motif`` - motif enrichedment results for control enriched motifs, both know and denovo results are included

- ``GSEA`` - GSEA analysis result based on the log2FoldChange for each nearby genes in the differential peak list

- ``cistrome_toolkit`` - cistrome_toolkit analysis result based on the treatment and control enriched differential peaks

- ``DEseq.normalized.counts.csv`` - DEseq normalized counts for each sample and each peak


FOLDER ``logs``
=============================================

Stores various log and error files.

- ``*.log`` files from bam file sorting and cleaing: Each log file is produced by the corresponding rule and contains debugging information as well as warnings and errors:

  - ``clean_bam`` logs for picard bam clean 
  - ``remove_duplicates`` logs for picard remove duplicate 
  - ``reorder`` logs for reorder the bam files
  - ``read_coutns`` for bedtools intersect to get the sample-peak count matrix

.. _workingWithPipeline:

Running *CoBRA*
******************

General remarks
==============================

*CoBRA* is programmed as a *Snakemake* pipeline. *Snakemake* is a bioinformatics workflow manager that uses workflows that are described via a human readable, Python based language. It offers many advantages to the user because each step can easily be modified, parts of the pipeline can be rerun, and workflows can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition or only minimal modifications. However, with great flexibility comes a price: the learning curve to work with the pipeline might be a bit higher, especially if you have no *Snakemake* experience. For a deeper understanding and troubleshooting errors, some knowledge of *Snakemake* is invaluable.

Simply put, *Snakemake* executes various *rules*. Each *rule* can be thought of as a single *recipe* or task such as sorting a file, running an R script, etc. Each rule has, among other features, a name, an input, an output, and the command that is executed. You can see in the ``Snakefile`` what these rules are and what they do. During the execution, the rule name is displayed, so you know exactly at which step the pipeline is at the given moment. Different rules are connected through their input and output files, so that the output of one rule becomes the input for a subsequent rule, thereby creating *dependencies*, which ultimately leads to the directed acyclic graph (*DAG*) that describes the whole workflow. You have seen such a graph in Section :ref:`workflow`.


.. _timeMemoryRequirements:

Executing *CoBRA* - Running times and memory requirements
===============================================================

*CoBRA* can be computationally demanding depending on the sample size and the number of peaks. In the following, we discuss various issues related to time and memory requirements and we provide some general guidelines that worked well for us.

.. warning:: We generally advise to run *CoBRA* in a cluster environment. For small analysis, a local analysis on your machine might work just fine (see the example analysis in the Git repository), but running time increases substantially due to limited amount of available cores.

Analysis size
---------------

We now provide a *very rough* classification into small, medium and large with respect to the sample size and the number of peaks:

- Small: Fewer than 10-15 samples, number of peaks not exceeding 50,000-80,000, normal read depth per sample
- Large: Number of samples larger than say 20 or number of peaks clearly exceeds 100,000, or very high read depth per sample
- Medium: Anything between small and large


Total running time
--------------------

Some notes regarding the total running time:

- the total running time is based on the number of samples, their read depth, the number of peaks, and the number of TF included in the analysis
- for small analysis such as the example analysis in the Git repository, running times are roughly 30 minutes with 2 cores for sorted bam files
- for large analysis, running time will be up to 2 hrs or so when executed on a cluster machine
- for the motif anslysis, running time will add 1 hr in additional to the running time above



.. _FAQs:

Frequently asked questions (FAQs)
****************************************

Here a few typical use cases, which we will extend regularly in the future if the need arises:

1. Why the CoBRA need to use an config.yaml, metasheet.csv to setup the run? Why not just simply use the command to setup the run?

  The unsupervised and supervised anlaysis of ChIP/ATAC-seq experiment need to have many paramaters, and could be very different from different experiement. The config.yaml and metasheet.csv allow user to save all paramaters that have been used in this run and allow others to reproduce the analysis when needed.

2. I have a problem running docker?

  Please try to go https://docs.docker.com/toolbox/faqs/troubleshoot/ to get the docker running.

3. I want to rerun a specific part of the pipeline only.

  This common scenario is also easy to solve: Just invoke *Snakemake* with ``--forcerun {rulename}``, where ``{rulename}`` is the name of the rule as defined in the Snakefile. *Snakemake* will then rerun the  specified run and all parts downstream of the rule. If you want to avoid rerunning downstream parts (think carefully about it, as there might be changes from the rerunning that might have consequences for downstream parts also), you can combine ``--forcerun`` with ``--until`` and specify the same rule name for both.

4. I want to modify the workflow.

  Simply add or modify rules to the Snakefile, it is as easy as that.


.. _docs-errors:


Handling errors
************************************************************


We tried to cover all cases for which *CoBRA* may fail, so please post an issue on our `Bitbucket Issue Tracker <https://bitbucket.org/cfce/cobra/issues>`_ if you believe you found a new problem.


Common errors
================

We here provide a list of some of the errors that can happen and that users reported to us. This list will be extended whenever a new problem has been reported.

1. KeyError in ``metasheet_setup.py``

  .. code-block:: Bash

    *** KeyError in line 9 of Snakefile ***
    File "Snakefile", line 9, in <module>
    File "metasheet_setup.py", line 19, in updateMeta
    File "metasheet_setup.py", line 19, in <dictcomp>
  ..

  .. note:: This particular message normally related to an mismatch between the sample names in ``config.yaml`` and ``metasheet.csv``.


  Simply check if the names are matched would solve this error.


2. rule ``heatmapSS_plot`` duplicate 'row.names' are not allowed

 
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

  This error normally happens when you have duplicate sample names in the metasheet.csv, the CoBRA do not allow duplicate smples names in the config.yaml and metasheet.csv.

  If you do not know what the error is, post an Issue in the `Bitbucket Issue Tracker <https://bitbucket.org/chrarnold/CoBRA>`_ tracker and we are hopefully able to help you quickly.

3. Data-specific errors

  Errors can also be due to incompatible data. For example, if a BAM file contains both single-end and paired-end reads (which is unusual, lots of programs may exit with errors for such data) and in the configuration file the parameter *pairedEnd* is set to true, the *repair* step from *Subread* will fail with an error message. In such a case, the single-end reads should either be removed from the BAM file (this is the preferred option, unless the majority of reads are single-end) or *pairedEnd* is set to false, the latter of which then treats all reads to be single-end (with the consequence that then, not fragments are counted, but just individual reads, which may result in different results due to altered number of counts).


Fixing the error
==============================

General guidelines
--------------------
After locating the error, fix it accordingly. We here provide some guidelines of different error types that may help you fixing the errors you receive:

- Errors related to erroneous input: These errors are easy to fix, and the error message should be indicative. If not, please let us know, and we improve the error message in the pipeline.
- Errors of technical nature: Errors related to memory, missing programs, R libraries etc can be fixed easily by making sure the necessary tools are installed and by executing the pipeline in an environment that provides the required technical requirements. For example, if you receive a memory-related error, try to increase the available memory. In a cluster setting, adjust the cluster configuration file accordingly by either increasing the default memory or (preferably) or by overriding the default values for the specific rule.
- Errors related to *Snakemake*: In rare cases, the error can be due to *Snakemake* (corrupt metadata, missing files, etc). If you suspect this to be the case, you may delete the hidden ``.snakemake`` directory in the folder from which you started the analysis. *Snakemake* will regenerate it the next time you invoke it then.
- Errors related to the input data: Error messages that indicate the problem might be located in the data are more difficult to fix, and we cannot provide guidelines here. Feel free to contact us.

Debugging R scripts to identify the cause of an error
--------------------------------------------------------------------
If an R script fails with a technical error such as ``caught segfault`` (a segmentation fault), you may want to identify the library or function call that causes the message in order to figure out which library to reinstall. To do so, open the R script that fails in *RStudio*, and execute the script line by line until you identify the line that causes the issue. Importantly, read the instructions in the section at the beginning of the script that is called ``SAVE SNAKEMAKE S4 OBJECT THAT IS PASSED ALONG FOR DEBUGGING PURPOSES``. Briefly, you simply have to make the *snakemake* object available in your R workspace, which contains all necessary information to execute the R script properly. Normally, *Snakemake* automatically loads that when executing a script. To do so, simply execute the line that is pasted there in R, it is something like this:

.. code-block:: R

  snakemake = readRDS("{outputFolder}/LOGS_AND_BENCHMARKS/checkParameters.R.rds")

Replace ``{outputFolder}`` by the folder you used for the analysis, and adjust the ``checkParameters`` part also accordingly. Essentially, you just have to provide the path to the corresponding file that is located in the ``LOGS_AND_BENCHMARKS`` subdirectly within the specified output directory.

Rerunning *Snakemake*
----------------------
After fixing the error, rerun *Snakemake*. *Snakemake* will continue at the point at which the error message occurred, without rerunning already successfully computed previous steps (unless specified otherwise).



Understanding and interpreting results
****************************************

Having results is exciting; however, as with most software, now the maybe even harder part starts: Understanding and interpreting the results. Let's first remind ourselves: The main goal of *CoBRA* is to aid in formulating testable hypotheses and ultimately improve the understanding of regulatory mechanisms that are driving the differences on a system-wide scale.

General notes
=================

  - Irrespective of whether or not you also used the classification mode, we recommend that the first thing to check is the Volcano plot PDF.
  - If a specific question is not addressed here, feel free to contact us. We ill then add it here.
  - Note that *CoBRA* captures differential accessibility, which does not necessarily imply a functional difference. See the publication for more discussion and details.
  - the significance as calculated by the empirical or analytical approach should not be over-interpreted from our point of view. We find the TF activity to be the more important measure.


Specifics for the basic mode
=================================

The following procedure may be useful as a rough guideline:
  - Start with the most stringent adjusted p-value threshold (0.01)
  - Categorize into one of the 3 following cases:
    - (a) There are no or almost none TFs significant: You may simply use a less stringent adjusted p-value threshold. If the least stringent adjusted p-value threshold (0.2) does also not have any or only very few significant TFs,  see the :ref:`FAQs` for possible explanations. In such a (rare) case, it might be worthwhile then to check the raw p-values instead of the adjusted ones.
    - (b) A few TFs are significant: You hit the sweet spot! Try to characterize and understand the TFs and whether they make biological sense for you. See also the notes for (c) below.
    - (c) A lot or the majority of TFs are significant (say more than 50 to 100): See the :ref:`FAQs` for possible explanations and how to best proceed.



Specifics for the classification mode
==========================================

- be aware of the limitations, see below

Limitations
-------------

As written in the publication, we note that *CoBRA* is prone to mis-classifying TFs that (1) act bifunctionally as activators and repressors in different genomic contexts or along with different co-factors, (2) are heavily regulated post-translationally, or (3) show little variation in RNA expression across the samples. Some of these mis-classifications may represent interesting subjects for future investigations.
Furthermore, if two TFs have similar motifs, which makes it difficult to distinguish them, *CoBRA* may have difficulties in classifying them correctly. Thus, for distinguishing the functional roles of TFs from the same motif-family, further biochemical experiments are needed.
