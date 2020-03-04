.. _docs-quickstart:

Try it out now!
============================================================

*CoBRA* runs on Linux and macOS and is even independent on the operating system if combined with ``Docker``. The following quick start briefly summarizes the necessary steps to install and use it.

Principally, there are two ways of installing *CoBRA* and the proper tools:

1a. **The "easy" way**: Using ``Docker`` and our preconfigured *CoBRA* containers that contain all necessary tools, R, and R libraries

  You only need to install ``Docker``. *CoBRA* supports Docker in Versions >=1.8. You can check whether you already have ``Docker`` installed by simply typing

  .. code-block:: Bash

    docker --version

  If docker is not installed, please install `docker <https://docs.docker.com/install/>`_ on your platform first. Once docker is installed, you can simply type following command to pull and run the CoBRA image.
  

  .. code-block:: Bash

    docker run --rm -v $PWD:/cobra -it cfce/cobra:latest
  
  .. note:: Make to read the section :ref:`docs-DockerNotes` properly!

1b. **The "more complicated" way**:  Install the necessary tools (*Snakemake*, *samtools*, *bedtools*, *Homer*, and *R* along with various packages).

  .. note:: Note that most tools in CoBRA are avaiable via conda.

  We recommend installing all tools via conda, in which case the installation then becomes as easy as download the :download:`environment.yml <environment.yml>` file. Then run the following command:

  .. code-block:: Bash

     conda env create -f environment.yml -n cobra

  If conda is not yet installed, follow the `installation instructions <https://conda.io/docs/user-guide/install/index.html>`_. Installation is quick and easy. Make sure to restrat the terminal after installation, so that *conda* is available.

  In addition, Cobra is needed along with following packages that ourside the conda framework. See the following instructions for each of the tools: `giggle  <https://github.com/ryanlayer/giggle>`_, `liftover <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver>`_.

2. **Activate cobra conda environment:**

  Once you have installed ``Cobra`` by docker or conda, you now can use the following command activate the environment
  
  .. code-block:: Bash

     source activate cobra
  
  .. note:: If install the cobra environemnt by docker, please read the :ref:`docs-DockerNotes` properly! Especially folllowing the instructions to mount the current folder inside into the containter.


3. **Clone the Git repository:**

  If install by docker, please run the following command to change the working directory:
   

  .. code-block:: Bash
   
     cd cobra
   
  Following command will get most recent version of codes that require to run cobra, this command need to run in a empty directory:

  .. code-block:: Bash

     git clone https://bitbucket.org/cfce/cobra.git .

  If you receive an error, *Git* may not be installed on your system. Please consult the internet on how to best install Git for your system.

3. **To run CoBRA with an example ChIP-Seq / ATAC-seq dataset, simply perform the following steps (see section**  :ref:`exampleDataset` **for dataset details)**:

  Download the example data for the GR-ChIP experiement, this will download all bam, bigwig and bed files that needed to run the example:

  .. code-block:: Bash

     snakemake download_example_GR_ChIP

  To test if the setup is correct, start a dryrun via the following command:

  .. code-block:: Bash

     snakemake all -np

  Once the dryrun is successful, start the analysis via the command (using 6 cores).

  .. code-block:: Bash

     snakemake all --cores 6

4. **To run your own analysis**, modify the files ``config.yaml`` and ``metasheet.csv``. See the instructions in the section `Run your own analysis`_ for more details.
5. **If your analysis finished successfully**, take a look into the ``analysis`` folder within your specified output directory, which contains the data and visualization of your analysis. If you received an error, take a look in Section :ref:`docs-errors` to troubleshoot.

.. _docs-prerequisites:

Prerequisites for the "easy" way
==================================

The only prerequisite here is that ``Docker`` must be installed on the system you want to run *CoBRA*. See above for details with respect to the supported versions etc.


Prerequisites for the "manual" way
=====================================

Note that most of this section is only relevant if you use *CoBRA* without ``Docker``. This section lists the required software and how to install them. As outlined in Section :ref:`docs-quickstart`, the easiest way is to install all of them via ``conda``. However, it is of course also possible to install the tools separately.


.. _docs-runOwnAnalysis:

Run your own analysis
============================================================

Running your own analysis is almost as easy as running the example analysis (see section :ref:`exampleDataset`). Carefully read and follow the following steps and notes:

1. Modify the file ``config.yaml`` accordingly. See Section :ref:`configurationFile` for details about the meaning of the parameters. Do not delete or rename any parameters or sections.
2. Change the ``metasheet.csv`` file that match the input data, in analogy to the file ``metasheet.csv`` from the example analysis, and refer to that in the file ``config.yaml`` (parameter ``bam``, ``bigwig``, ``bed``)
3. Activate the cobra environment and start a dryrun with the following command
   
   .. code-block:: Bash
      
      source activate cobra
      snakemake all -np
   
   If dryrun is successfull, proceeding with the following to start (using 6 cores).
   
   .. code-block:: Bash
      
      snakemake all --cores 6 
      
4. Since running the pipeline is often computationally demanding, read Section :ref:`timeMemoryRequirements` and decide on which machine to run the pipeline. In most cases, we recommend running *CoBRA* in a cluster environment (see Section :ref:`clusterEnvironment` for details). The pipeline is written in Snakemake, and we strongly suggest to also read Section :ref:`workingWithPipeline` to get a basic understanding of how the pipeline works.


.. _docs-DockerNotes:

Adaptations and notes when running with Docker
============================================================
 With ``Docker``, the *CoBRA* workflow will be executed in pre-configured isolated container that contain all necessary tools.  To use it, you only have to add the following arguments when you initial the docker container:

.. code-block:: Bash

   docker run --rm -v $PWD:/cobra -it cfce/cobra:latest

1. ``--rm``: This option will help delete the container immediately after it exits. This helps to prevent having to clean up containers after finish runing the workflow.

2. ``-v``: The ``-v`` flag mounts the current directory ``$PWD`` into /cobra in the container. You need to make all directories that contain files that are referenced in the *CoBRA* config file available within the container. If you reference additional files, simply add multiple ``-v`` flags to the mount path (use the space to separate them). For example, if you reference the files ``/mnt/home/user1/AR_ChIP.bam`` and ``/mnt/home/user1/AR_ChIP.bed`` in the configuration file file, you may add ``-v /mnt/home/user1:/mnt/home/user1 `` or even just ``-v /mnt:/mnt`` to the bind path.

  .. note:: We note again that within a Docker container, you cannot access paths outside of the directory from where you started executing Snakemake. If you receive errors indecate that a directory does not exist even though you can cd into it, you most likely forgot to include the path this folder or a parent path as part of the ``-v`` option.

3. ``-it``: The ``-it`` options allows you to interact with the container’s shell and run any command inside of it.

4. ``cfce/cobra:latest``: ``cfce/cobra`` is the name of the container that we create in the dockerhub. ``lastest`` is the version of the container.

Once you start runing the cobra containter, it's bash shell will be attached to the terminal, and the command prompt will change:

.. code-block:: Bash

   (base) root@5d8bf16cd2cb:/#

Above command prompt change means you have suceefully start the cotainer of cobra, you may proceed to run the example or your own data.

You do not have to, but you may go through the following tutorial related to ``Docker``, this will help you understand the docker better. For more details, see `here <https://docker-curriculum.com/>`_.
