
.. _docs-quickstart:

*****
Tutorials
*****

The three tutorials below gives you a taste of different capabilities of our CoBRA workflow. The data sets we used are, a GR ChIP-seq data set from the ENCODE project, a set of H3K27ac ChIP-seq data from colon cancer cell lines, and an ATAC-seq experiment on HL-60 promyelocytes differentiating into macrophages, each of which demonstrate a subfield of expertise of our CoBRA pipeline. 

GP_ChIP Tutorial
============================================================

####Background####

Our first tutorial includes an in depth analysis of the sample GR_ChIP dataset.

1. **Initiate Environment**: 
  
  Use the following command to activate the cobra environment:
  
  .. code-block:: Bash

    source activate cobra


2. **Retrieve the Latest Version of Cobra:**

  If installed using docker, run the following command to change the working directory. Otherwise, skip to next command:
   
  .. code-block:: Bash
   
     cd cobra
   
  Create a new directory, change to the new working directory and pull the latest version of CoBRA using ``git clone`` :

  .. code-block:: Bash

     mkdir gr_chip
     git clone https://bitbucket.org/cfce/cobra.git .

  If you receive an error, *Git* may not be installed on your system. Please consult the internet on how to best install Git for your system.


3. **Download the GR_ChIP sample dataset:**

  Please use the following command to download the GR_ChIP sample dataset. This dataset is of moderate size (3.9 G) and may take 5-10 minutes to download. 

  .. code-block:: Bash
   
     snakemake download_example_GR_ChIP
   
  The following command will pull the most recent version code that is required to run cobra. This command must be run in a empty directory:

  .. code-block:: Bash

     git clone https://bitbucket.org/cfce/cobra.git .

  If you receive an error, *Git* may not be installed on your system. Please consult the internet on how to best install Git for your system.

4. **To run CoBRA with an example ChIP-Seq / ATAC-seq dataset, execute the following steps (see section**  :ref:`exampleDataset` **for details)**:

  Download the example data for the GR-ChIP experiement. This command will download all bam, bigwig and bed files that are needed to run the example:

  .. code-block:: Bash

     snakemake download_example_GR_ChIP

  To check if the setup is correct, begin a dry run via the following command:

  .. code-block:: Bash

     snakemake all -np

  Once the dry run completes without errors, start the analysis using the command (using 6 cores).

  .. code-block:: Bash

     snakemake all --cores 6

5. **To run CoBRA on your experiment**, setup the files ``config.yaml`` and ``metasheet.csv`` according to your own experiment. Instructions can be found in the section `Run CoBRA on your experiment`_.
6. **If CoBRA runs succesfully**, explore the ``analysis`` folder which contains the data and visualization of your analysis. If you encountered an error, look in Section :ref:`docs-errors` to troubleshoot.

.. _docs-prerequisites:

Dependencies for quick installation
==================================

For quick instalation, the only dependency required is ``Docker``. Once ``Docker`` is installed, *CoBRA* dependencies can be pulled by the ``Docker`` container. Details on how to pull the ``Docker`` image can be found above.


Dependencies for manual installation
=====================================

This section is only relevant if you do not use ``Docker`` to run *CoBRA*. All dependencies are listed above in section 1b. 


.. _docs-runOwnAnalysis:

Run CoBRA on your experiment
============================================================

See section :ref:`exampleDataset`for a guide to running the example analysis. Running CoBRA on your own data is straitforward. To do so, please execute the following steps:

1. Modify the file ``config.yaml`` accordingly. See Section :ref:`configurationFile` for details about the meaning of the parameters. Do not delete or rename any parameters or sections.
2. Change the ``metasheet.csv`` file to match the input data. Just as in the example dataset, the metasheet contains the same sample names as the ``config.yaml`` file, the same must be done when running your own analysis.

3. Activate the cobra environment and start a dry run with the following command
   
   .. code-block:: Bash
      
      source activate cobra
      snakemake all -np
   
   If dry run is successfull, proceeding with the following to start (using 6 cores).
   
   .. code-block:: Bash
      
      snakemake all --cores 6 
      
4. Running *CoBRA* is computationally demanding (see section :ref:`timeMemoryRequirements`). As such, we suggest running *CoBRA* in a multicore machine capable of handling parallelization. 


.. _docs-DockerNotes:

Notes for running with Docker
============================================================
 With ``Docker``, the *CoBRA* workflow will be executed in pre-configured isolated container that contains all dependency tools. You only need to pay attention to the following arguments when running *CoBRA* in the ``Docker`` container.

.. code-block:: Bash

   docker run --rm -v $PWD:/cobra -it cfce/cobra:latest

1. ``--rm``: This option will help delete the container immediately after it exits. This helps to prevent having to clean up containers after the workflow has finished running.

2. ``-v``: The ``-v`` flag mounts the current directory ``$PWD`` into /cobra in the container. You need to make all directories that contain files that are referenced in the *CoBRA* config file available within the container. If you reference additional files, simply add multiple ``-v`` flags to the mount path (use the space to separate them). For example, if you reference the files ``/mnt/home/user1/AR_ChIP.bam`` and ``/mnt/home/user1/AR_ChIP.bed`` in the configuration file file, you may add ``-v /mnt/home/user1:/mnt/home/user1 `` or even just ``-v /mnt:/mnt`` to the bind path.
  
  .. note:: We emphasize that within a Docker container, files outside of the directory from where you started executing *CoBRA* are not accessible. You will receive errors if you forgot to include the path of this folder as part of the ``-v`` option.

3. ``-it``: The ``-it`` option allows you to interact with the container’s shell and run any command inside of it.

4. ``cfce/cobra:latest``: ``cfce/cobra`` is the name of the container that we created in the dockerhub. ``lastest`` is the version of the container.

Once you start runing the *CoBRA* containter, it's bash shell will be attached to the terminal, and the command prompt will change:

.. code-block:: Bash

   (base) root@5d8bf16cd2cb:/#

The above command prompt change means that you have suceefully started the cotainer of *CoBRA*, and you may proceed to run the example or your own data.

You do not have to, but you may go through the following tutorial related to ``Docker``. This will help you gain a better understanding of ``Docker``. For more details, see `here <https://docker-curriculum.com/>`_.
