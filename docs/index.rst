.. CoBRA documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
=========================================
Welcome to the documentation of *CoBRA*!
=========================================

===========================
 *CoBRA*: Containerized Bioinformatics workflow for Reproducible ChIP/ATAC-seq Analysis
===========================

CoBRA is a comprehensive ChIP/ATAC‐seq analysis workflow built using snakemake and Docker which allows for a scalable, reproducible, portable and easy-to-use pipeline.

CoBRA combines the use of several dozen ChIP/ATAC‐seq tools, suites, and packages to create a complete pipeline that takes ChIP/ATAC‐seq data and produces unsupervised analyses, differential peak calling, and downstream pathway analysis. In addition, CoBRA has been outfitted with several recently published tools that allow for better normalziation and CNV correction. The results are compiled in a simple and highly organized output format containing key figures and tables.

The documentation is organized into the following two sections:


.. toctree::
   :maxdepth: 2
   :caption: Installation and Quick Start 

   install.rst

.. toctree::
   :maxdepth: 2
   :caption: Workflow Details

   detail.rst

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials.rst

