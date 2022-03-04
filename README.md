[title]: - "High-Throughput BWA Read Mapping"
[TOC]


## Introduction

This tutorial focuses on a subset of the [Data Carpentry Genomics workshop curriculum](https://datacarpentry.org/genomics-workshop/) - specifically, this page cover's how to run a BWA workflow on OSG resources. It will use the same general flow as the BWA segment of the Data Carpentry workshop with minor adjustments. The goal of this tutorial is to learn how to convert an existing BWA workflow to run on the OSPool.  

## Get Tutorial Files

Logged into the submit node, we will run the tutorial command, that will 
create a folder for our analysis, as well as some sample files. 

```
tutorial bwa
```

## Install and Prepare BWA
First, we need to install BWA, also called Burrows-Wheeler Aligner. To do this, we will create and navigate to a new folder in our /home directory called `software`. We will then follow the developer's instructions (https://github.com/lh3/bwa) for using `git clone` to clone the software and then build the tool using `make`. 

```
cd ~/tutorial-bwa
cd software
git clone https://github.com/lh3/bwa.git
cd bwa
make
```

Next, BWA needs to be added to our PATH variables, to test if the installation worked: 

```
export PATH=$PATH:/home/$USER/tutorial-bwa/software/bwa/
``` 

To check that BWA has been installed correctly, type `bwa`. You should receive output similar to the following: 

```
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1198-dirty
Contact: Heng Li <hli@ds.dfci.harvard.edu>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
...

```

Now that we have successfully installed `bwa`, we will create a portable compressed tarball of this software so that it is smaller and quicker to transport when we submit our jobs to the OSPool. 

```
cd ~/tutorial-bwa/software
tar -czvf bwa.tar.gz bwa
```

Checking the size of this compressed tarball using `ls -lh bwa.tar.gz` reveals the file is approximately 4MB. The tarball should stay in /home.


## Download Data to Analyze

Now that we have installed BWA, we need to download data to analyze. For this tutorial, we will be downloading data used in the Data Carpentry workshop. This data includes both the genome of Escherichia coli (E. coli) and paired-end RNA sequencing reads obtained from a study carried out by Blount et al. published in [PNAS](http://www.pnas.org/content/105/23/7899). Additional information about how the data was modified in preparation for this analysis can be found on the [Data Carpentry's workshop website](https://datacarpentry.org/wrangling-genomics/aio.html).

``` 
cd ~/tutorial-bwa
./download_data.sh
```
Investigating the size of the downloaded genome by typing:

```
ls -lh data/ref_genome/
```

reveals the file is 1.4 MB. Therefore, this file should remain in /home and does not need to be moved to /public. We should also check the trimmed fastq paired-end read files: 

```
ls -lh data/trimmed_fastq_small
```

Once everything is downloaded, make sure you're still in the `tutorial-bwa` directory. 
```
cd ~/tutorial-bwa
```

## Run a Single Test Job

Now that we have all items in our analysis ready, it is time to submit a single test job to map our RNA reads to the E. coli genome. For a single test job, we will choose a single sample to analyze. In the following example, we will align both the forward and reverse reads of SRR2584863 to the E. coli genome. Using a text editor such as `nano` or `vim`, we can create an example submit file for this test job called `bwa-test.sub` containing the following information:

```
universe    = vanilla
executable  = bwa-test.sh
# arguments = 

# need to transfer bwa.tar.gz file, the reference
# genome, and the trimmed fastq files
transfer_input_files = software/bwa.tar.gz, data/ref_genome/ecoli_rel606.fasta.gz, data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq, data/trimmed_fastq_small/SRR2584863_2.trim.sub.fastq

log         = TestJobOutput/bwa_test_job.log
output      = TestJobOutput/bwa_test_job.out
error       = TestJobOutput/bwa_test_job.error

+JobDurationCategory = "Medium"
request_cpus    = 1
request_memory  = 2GB
request_disk    = 1GB

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
requirements = (OSGVO_OS_STRING == "RHEL 7")

queue 1
```
You will notice that the .log, .out, and .error files will be saved to a folder called `TestJobOutput`. We need to create this folder using `mkdir TestJobOutput` before we submit our job. 

We will call the script for this analysis `bwa-test.sh` and it should contain the following information: 

```
#!/bin/bash
# Script name: bwa-test.sh

echo "Unpacking software"
tar -xzf bwa.tar.gz

echo "Setting PATH for bwa" 
export PATH=$_CONDOR_SCRATCH_DIR/bwa/:$PATH

echo "Indexing E. coli genome"
bwa index ecoli_rel606.fasta.gz

echo "Starting bwa alignment for SRR2584863"
bwa mem ecoli_rel606.fasta.gz SRR2584863_1.trim.sub.fastq SRR2584863_2.trim.sub.fastq > SRR2584863.aligned.sam

echo "Done with bwa alignment for SRR2584863!"

echo "Cleaning up files generated from genome indexing"
rm ecoli_rel606.fasta.gz.amb
rm ecoli_rel606.fasta.gz.ann
rm ecoli_rel606.fasta.gz.bwt
rm ecoli_rel606.fasta.gz.pac
rm ecoli_rel606.fasta.gz.sa
```

We can submit this single test job to HTCondor by typing: 

```
condor_submit bwa-test.sub
```

To check the status of the job, we can use `condor_q`. 

Upon the completion of the test job, we should investigate the output to ensure that it is what we expected and also review the `.log` file to help optimize future resource requests in preparation for scaling up. 

For example, when we investigate the `bwa_test_job.log` file created in this analysis, at the bottom of the file we see a resource table: 

```       
       Partitionable Resources :    Usage  Request Allocated 
           Cpus                 :                 1         1 
           Disk (KB)            :   253770  1048576  27945123 
           Memory (MB)          :      144     2048      2500
```

Here we see that we used less than half of both the disk space and memory we requested. In future jobs, we should request a smaller amount of each resource, such as 0.5 GB of disk space and 0.5 GB of memory. Prior to scaling up our analysis, we should run additional test jobs using these resource requests to ensure that they are sufficient to allow our job to complete successfully.


## Scaling Up to Analyze Multiple Samples

In preparation for scaling up, please review our [guide on how to scale up after a successful test job](https://support.opensciencegrid.org/support/solutions/articles/12000076552-scaling-up-after-success-with-test-jobs) and how to 
[easily submit multiple jobs with a single submit file](https://support.opensciencegrid.org/support/solutions/articles/12000073165-easily-submit-multiple-jobs).

After reviewing how to submit multiple jobs with a single submit file, it is possible to determine that the most appropriate way to submit multiple jobs for this analysis is to use `queue <var> from <list.txt>`. 

To use this option, we first need to create a file with just the sample names/IDs that we want to analyze. To do this, we want to cut all information after the "_" symbol to remove the forward/reverse read information and file extensions. For example, we want SRR2584863_1.trim.sub.fastq to become just SRR2584863. 

We will save the sample names in a file called `samples.txt`:

```
cd ~/tutorial-bwa
cd data/trimmed_fastq_small/
ls *.fastq | cut -f 1 -d '_' | uniq > samples.txt
cd ~/tutorial-bwa
```

Now, we can create a new submit file called `bwa-alignment.sub` to queue a new job for each sample. To make it simpler to start, you can copy the `bwa-test.sub` file (`cp bwa-test.sub bwa-alignment.sub`) and modify it. 

```
universe    = vanilla
executable  = bwa-alignment.sh
arguments   = $(sample)

transfer_input_files = software/bwa.tar.gz, data/ref_genome/ecoli_rel606.fasta.gz, data/trimmed_fastq_small/$(sample)_1.trim.sub.fastq, data/trimmed_fastq_small/$(sample)_2.trim.sub.fastq

transfer_output_remaps = "SRR2584863.aligned.sam=results/SRR2584863.aligned.sam; SRR2584866.aligned.sam=results/SRR2584866.aligned.sam; SRR2589044.aligned.sam=results/SRR2589044.aligned.sam"

log         = log/bwa_$(sample)_job.log
output      = output/bwa_$(sample)_job.out
error       = error/bwa_$(sample)_job.error

+JobDurationCategory = "Medium"
request_cpus    = 1 
request_memory  = 0.5GB
request_disk    = 0.5GB

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
requirements = (OSGVO_OS_STRING == "RHEL 7")

queue sample from data/trimmed_fastq_small/samples.txt
```

In addition to restructuring our submit file to queue a new job for each sample, it is also advantageous to have our standard output, log, and error files saved to dedicated folders called "log", "output", and "error" to help keep our output files organized.  Therefore, we need to make these folders in our /home directory prior to submitting our job. We will also create an additional folder to store our aligned sequencing files in a folder called `results`:

```
mkdir log
mkdir output
mkdir error
mkdir results
```

To store the aligned sequencing files in the `results` folder, we can add the `transfer_output_remaps` feature to our submit file. This feature allows us to specify a name and a path to save our output files in the format of "file1 = path/to/save/file2", where file1 is the origional name of the document and file2 is the name that we want to save the file using. In the example above, we do not change the name of the resulting output files. This feature also helps us keep an organized working space, rather than having all of our resulting sequencing files be saved to our /home directory. 

Once our submit file has been updated, we can update our script to look like and call it something like `bwa-alignment.sh`: 

```
#!/bin/bash
# Script name: bwa-alignment.sh

echo "Unpackage software"
tar -xzf bwa.tar.gz

echo "Set PATH for bwa" 
export PATH=$_CONDOR_SCRATCH_DIR/bwa/:$PATH

# Renaming first argument
sample=$1

echo "Index E.coli genome"
bwa index ecoli_rel606.fasta.gz

echo "Starting bwa alignment for ${sample}"
bwa mem ecoli_rel606.fasta.gz ${sample}_1.trim.sub.fastq ${sample}_2.trim.sub.fastq > ${sample}.aligned.sam

echo "Done with bwa alignment for ${sample}!"

echo "Cleaning up workspace"
rm ecoli_rel606.fasta.gz.amb
rm ecoli_rel606.fasta.gz.ann
rm ecoli_rel606.fasta.gz.bwt
rm ecoli_rel606.fasta.gz.pac
rm ecoli_rel606.fasta.gz.sa
```

Once ready, we can submit our job to HTCondor by using `condor_submit bwa-alignment.sub`. 

When we type `condor_q`, we see that three jobs have entered the queue (one for each of our three experimental samples).

When our jobs are completed, we can confirm that our alignment output results files were created by typing:

```
ls -lh results/*
``` 

We can also investigate our log, error, and output files in their respective folders to ensure we obtained the resulting output of these files that we expected. 



_For more information about running bioinformatics workflows on the OSG, we recommend our [BLAST tutorial](https://support.opensciencegrid.org/support/solutions/articles/12000062020-running-a-blast-workflow) as well as our [Samtools](https://support.opensciencegrid.org/support/solutions/articles/12000074984-example-software-compilation) instillation guide._
