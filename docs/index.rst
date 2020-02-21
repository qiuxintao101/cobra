.. diffTF documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================================
Welcome to the documentation of *CoBRA*!
=========================================

CoBRA - Containerized Bioinformatics workflow for Reproducible ChIP/ATAC-seq Analysis

CoBRA is a comprehensive ChIP/ATAC‐seq analysis tool built using snakemake and Docker which allows for escalable, reproducible, portable and easy-to-use workflows.

CoBRA combines the use of several dozen ChIP/ATAC‐seq tools, suites, and packages to create a complete pipeline that takes ChIP/ATAC‐seq analysis from unsupervised analyses, differential peak calling, and downstream pathway analysis. In addition, CoBRA has been outfitted with several recently published tools that allow for better normalziation and CNV correction. The results are compiled in a simple and highly visual report containing the key figures to explain the analysis, and then compiles all of the relevant files, tables, and pictures into an easy to navigate folder.


This site is organized into the following three parts:


.. toctree::
   :maxdepth: 2
   :caption: Quick Start and Installation

   install.rst

.. toctree::
   :maxdepth: 2
   :caption: Workflow Details

   detail.rst

.. toctree::
  :maxdepth: 2
  :caption: Project Information

  projectInfo.rst
