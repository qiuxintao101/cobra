.. CoBRA documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===========================
 *CoBRA*: Containerized Bioinformatics workflow for Reproducible ChIP/ATAC-seq Analysis
===========================

CoBRA is a comprehensive ChIP/ATAC‐seq analysis workflow built using snakemake and Docker which allows for a scalable, reproducible, portable and easy-to-use pipeline.

CoBRA combines the use of several dozen ChIP/ATAC‐seq tools, suites, and packages to create a complete pipeline that takes ChIP/ATAC‐seq data and produces unsupervised analyses, differential peak calling, and downstream pathway analysis. In addition, CoBRA has been outfitted with several recently published tools that allow for better normalziation and CNV correction. The results are compiled in a simple and highly organized output format containing key figures and tables.


.. _citation:

Citation
============================

If you use this pipeline, please cite the following reference:

Xintao Qiu*, Avery S. Feit*, Ariel Feiglin*, Yingtian Xie, Nikolas Kesten, Len Taing, Joseph Perkins, Shengqing Gu, Yihao Li, Paloma Cejas, Ningxuan Zhou, Rinath Jeselsohn, Myles Brown, X. Shirley Liu, Henry W. Long. CoBRA: Containerized Bioinformatics Workflow for Reproducible ChIP/ATAC-seq Analysis. 2021. Genomics, Proteomics & Bioinformatics.

Open Access. DOI: `https://doi.org/10.1016/j.gpb.2020.11.007 <https://doi.org/10.1016/j.gpb.2020.11.007>`_


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

License
GNU General Public License v2.0
