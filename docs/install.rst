.. _docs-quickstart:

Try it out now!
============================================================

*CoBRA* runs on Linux and macOS and is even independent on the operating system if combined with ``Docker``. The following quick start briefly summarizes the necessary steps to install and use it.

Principally, there are two ways of installing *CoBRA* and the proper tools:

1a. **The "easy" way**: Using ``Docker`` and our preconfigured *CoBRA* containers that contain all necessary tools, R, and R libraries

  You only need to install ``Docker``. *CoBRA* supports Docker in Versions >=1.8. You can check whether you already have ``Docker`` installed by simply typing

  .. code-block:: Bash

    docker --version

  If docker is not installed, please install `docker <https://docs.docker.com/install/>`_ on your platform first. Once docker is installed, you can simply type following command to pull the CoBRA image.
  
   .. code-block:: Bash

    docker pull cfce/cobra:latest

  .. note:: Make to read the section :ref:`docs-DockerNotes` properly!

1b. **The "more complicated" way**:  Install the necessary tools (*Snakemake*, *samtools*, *bedtools*, *Homer*, and *R* along with various packages).

  .. note:: Note that most tools in CoBRA are avaiable via conda.

  We recommend installing all tools via conda, in which case the installation then becomes as easy as

  .. code-block:: Bash

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install snakemake bedtools samtools subread

  If conda is not yet installed, follow the `installation instructions <https://conda.io/docs/user-guide/install/index.html>`_. Installation is quick and easy. Make sure to restrat the terminal after installation, so that *conda* is available.

  .. note:: You do not need to uninstall other Python installations or packages in order to use conda. Even if you already have a system Python, another Python installation from a source such as the macOS Homebrew package manager and globally installed packages from pip such as pandas and NumPy, you do not need to uninstall, remove, or change any of them before using conda.

  If you want to install the tools manually and outside of the conda framework, see the following instructions for each of the tools: `snakemake  <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io/en/latest/content/installation.html>`_, `Subread <http://subread.sourceforge.net>`_.

  In addition, *R* is needed along with various packages (see below for details).

2. **Clone the Git repository:**

    .. code-block:: Bash

      git clone https://git.embl.de/grp-zaugg/CoBRA.git

    If you receive an error, *Git* may not be installed on your system. If you run Ubuntu, try the following command:

    .. code-block:: Bash

      sudo apt-get install git

    For macOS, there are multiple ways of installing it. If you already have *Homebrew* (http://brew.sh) installed, simply type:

    .. code-block:: Bash

      brew install git

    Otherwise, consult the internet on how to best install Git for your system.

3. **To run CoBRA with an example ATAC-Seq / RNA-seq dataset for 50 TF, simply perform the following steps (see section**  :ref:`exampleDataset` **for dataset details)**:

  * Change into the ``example/input`` directory within the Git repository

      .. code-block:: Bash

        cd CoBRA/example/input

  * Download the data via the download script

        .. code-block:: Bash

          sh downloadAllData.sh

  * To test if the setup is correct, start a dryrun via the first helper script

        .. code-block:: Bash

          sh startAnalysisDryRun.sh

  * Once the dryrun is successful, start the analysis via the second helper script.

    .. code-block:: Bash

      sh startAnalysis.sh

    If you want to include ``Docker`` (which we strongly recommend), simply edit the file and add the ``--use-Docker`` and ``--Docker-args`` command line arguments in addition to the other arguments (see the Snakemake documentation and the section :ref:`docs-DockerNotes` for more details).

    Thus, the command you execute should look like this:

        .. code-block:: Bash

          snakemake --snakefile ../../src/Snakefile --cores 2 --configfile config.json \
           --use-Docker --Docker-args "--bind /your/CoBRA/path"

    Read in section :ref:`docs-DockerNotes` about the ``--bind`` option and what ``/your/CoBRA/path`` means here , it is actually very easy!

    You can also run the example analysis with all TF instead of only 50. For this, simply modify the ``TF`` parameter and set it to the special word ``all`` that tells *CoBRA* to use all recognized TFs instead of a specific list only (see section :ref:`parameter_TFs` for details).

4. **To run your own analysis**, modify the files ``config.json`` and ``sampleData.tsv``. See the instructions in the section `Run your own analysis`_ for more details.
5. **If your analysis finished successfully**, take a look into the ``FINAL_OUTPUT`` folder within your specified output directory, which contains the summary tables and visualization of your analysis. If you received an error, take a look in Section :ref:`docs-errors` to troubleshoot.

.. _docs-prerequisites:

Prerequisites for the "easy" way
==================================

The only prerequisite here is that Snakemake and ``Docker`` must be installed on the system you want to run *CoBRA*. See above for details with respect to the supported versions etc. For details how to install Snakemake, see below.


Prerequisites for the "manual" way
=====================================

Note that most of this section is only relevant if you use Snakemake without ``Docker``. This section lists the required software and how to install them. As outlined in Section :ref:`docs-quickstart`, the easiest way is to install all of them via ``conda``. However, it is of course also possible to install the tools separately.

Snakemake
--------------------------

Please ensure that you have at least version 5.3 installed. Principally, there are `multiple ways to install Snakemake <http://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. We recommend installing it, along with all the other required software, via conda.

*samtools*, *bedtools*, *Subread*
----------------------------------

In addition, `samtools <http://www.htslib.org/download>`_, `bedtools <http://bedtools.readthedocs.io>`_ and `Subread <http://subread.sourceforge.net>`_ are needed to run *CoBRA*. We recommend installing them, along with all the other required software, via conda.


R and R packages
--------------------------

A working ``R`` installation is needed and a number of packages from either CRAN or Bioconductor have to be installed.  Type the following in ``R`` to install them:

.. code-block:: R

  install.packages(c("checkmate", "futile.logger", "tidyverse", "reshape2", "RColorBrewer", "ggrepel", "lsr", "modeest", "boot", "grDevices", "pheatmap", "matrixStats", "locfdr"))

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install(c("limma", "vsn", "csaw", "DESeq2", "DiffBind", "geneplotter", "Rsamtools", "preprocessCore", "apeglm"))


.. _docs-runOwnAnalysis:

Run your own analysis
============================================================

Running your own analysis is almost as easy as running the example analysis (see section :ref:`exampleDataset`). Carefully read and follow the following steps and notes:

1. Copy the files ``config.json`` and ``startAnalysis.sh`` to a directory of your choice.
2. Modify the file ``config.json`` accordingly. For example, we strongly recommend running the analysis for all TF instead of just 50 as for the example analysis. For this, simply change the parameter â€œTFsâ€ to â€œallâ€. See Section :ref:`configurationFile` for details about the meaning of the parameters. Do not delete or rename any parameters or sections.
3. Create a **tab-separated** file that defines the input data, in analogy to the file ``sampleData.tsv`` from the example analysis, and refer to that in the file ``config.json`` (parameter ``summaryFile``)
4. Adapt the file ``startAnalysis.sh`` if necessary (the exact command line call to Snakemake and the various Snakemake-related parameters). If you run with Docker, see the section below for modifications.
5. Since running the pipeline is often computationally demanding, read Section :ref:`timeMemoryRequirements` and decide on which machine to run the pipeline. In most cases, we recommend running *CoBRA* in a cluster environment (see Section :ref:`clusterEnvironment` for details). The pipeline is written in Snakemake, and we strongly suggest to also read Section :ref:`workingWithPipeline` to get a basic understanding of how the pipeline works.


.. _docs-DockerNotes:

Adaptations and notes when running with Docker
============================================================
 With ``Docker``, each rule will be executed in pre-configured isolated containers that contain all necessary tools.  To enable it, you only have to add the following arguments when you execute Snakemake:

1. ``--use-Docker``: Just type it like this!

2. ``--Docker-args``: You need to make all directories that contain files that are referenced in the *CoBRA* configuration file available within the container also. By default, only the directory and subdirectories from which you start the analysis are automatically mounted inside the container. Since the *CoBRA* source code is outside the ``input`` folder for the example analysis, however, at least the root directory of the Git repository has to be mounted. This is actually quite simple! Just use ``--Docker-args "--bind /your/CoBRA/path"`` and replace ``/your/CoBRA/path`` with the root path in which you cloned the *CoBRA* Git repository (the one that has the subfolders ``example``, ``src`` etc.). If you reference additional files, simply add one or multiple directories to the bind path (use the comma to separate them). For example, if you reference the files ``/g/group1/user1/mm10.fa`` and ``/g/group2/user1/files/bla.txt`` in the configuration file file, you may add ``/g/group1/user1,/g/group2/user1/files`` or even just ``/g`` to the bind path (as all files you reference are within ``/g``).

  .. note:: We note again that within a Docker container, you cannot access paths outside of the directory from where you started executing Snakemake. If you receive errors in the ``checkParameterValidity`` rule that a directory does not exist even though you can cd into it, you most likely forgot to include the path this folder or a parent path as part of the ``bind`` option.

3. ``--Docker-prefix /your/directory`` (optional): You do not have to, but you may want to add the ``--Docker-prefix`` argument to store all ``Docker`` containers in a central place (here: ``/your/directory``) instead of the local ``.snakemake`` directory. If you intend to run multiple *CoBRA* analyses in different folders, you can save space and time because the containers won't have to be downloaded each time and stored in multiple locations.

Please read the following additional notes and warnings related to ``Docker``:

- .. warning:: If you use ``Docker`` version 3, make sure you have at least version 3.0.3 installed, as there was an issue with Snakemake and particular ``Docker`` versions. For more details, see `here <https://bitbucket.org/snakemake/snakemake/issues/1017/snakemake-process-suspended-upon-execution>`_.