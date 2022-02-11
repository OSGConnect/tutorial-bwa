[TOC]


# Introduction
This tutorial focuses on a subset of the [Data Carpentry Genomics workshop curriculum](https://datacarpentry.org/genomics-workshop/) - specifically, this page cover's how to run a BWA workflow on OSG resources. It will use the same general flow as the BWA segment of the Data Carpentry workshop with minor adjustments. The goal of this tutorial is to learn how to convert an existing BWA workflow to run on the OS Pool.  


# Install and Prepare BWA
First, we need to install BWA, also called Burrows-Wheeler Aligner. To do this, we will create and navigate to a new folder in our /home directory called `software`. We will then follow the creator's instructions (https://github.com/lh3/bwa) for using `git clone` to clone the software and then build the tool using `make`. 

```
mkdir software
cd software
git clone https://github.com/lh3/bwa.git
cd bwa; make
```

Next, BWA needs to be added to our PATH variables. To do so, use the following commands but replace "userID" with your own user ID. 

```
export PATH=$PATH:/home/userID/software/bwa/
source ~/.bashrc
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

Now that we have successfully installed `bwa`, we will create a portable compressed tarball of this software so that it is smaller and quicker to transport when we submit our jobs to the OS Pool. 

```
cd ~/software
tar -czvf bwa.tar.gz bwa
```

Checking the size of this compressed tarball using `ls -lh bwa.tar.gz` reveals the file is approximately 4MB, which means it should stay in /home. 


# Download Data to Analyze
Now that we have installed BWA, we need to download data to analyze. For this tutorial, we will be downloading data used in the Data Carpentry workshop. This data includes both the genome of Escherichia coli (E. coli) and paired-end RNA sequencing reads obtained from a study carried out by Zachary D. Blount, Christina Z. Borland, and Richard E. Lenski published in [PNAS](http://www.pnas.org/content/105/23/7899). Additional information about how the data was modified in preparation for this analysis can be found on the [Data Carpentry's workshop website](https://datacarpentry.org/wrangling-genomics/aio.html).

``` 
cd ~
mkdir -p data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
```
Investigating the size of the downloaded genome by typing:

```
ls -lh data/ref_genome/
```

reveals the file is 1.4 MB. Therefore, this file should remain in /home and does not need to be moved to /public. 

Next, we will download a small subset of quality-trimmed fastq paired-end read files. 

```
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mv sub/ ~/data/trimmed_fastq_small
rm sub.tar.gz
```


# Run a Single Test Job
Now that we have all items in our analysis ready, it is time to submit a single test job to map our RNA reads to the E. coli genome. For a single test job, we will choose a single sample to analyze. In the following example, we will analyze the forward and reverse reads of both the forward and reverse reads of SRR2584863. An example submit file for this test job may be called `bwa-test.sub` and may look like: 

```
universe    = vanilla
executable  = bwa-test.sh

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
requirements = (OSGVO_OS_STRING == "RHEL 7")
+JobDurationCategory = "Medium"

transfer_input_files = software/bwa.tar.gz, data/ref_genome/ecoli_rel606.fasta.gz, data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq, data/trimmed_fastq_small/SRR2584863_2.trim.sub.fastq

# arguments = 

log         = TestJobOutput/bwa_test_job.log
output      = TestJobOutput/bwa_test_job.out
error       = TestJobOutput/bwa_test_job.error

request_cpus    = 1
request_memory  = 2GB
request_disk    = 1GB

queue 1
```
You will notice that the .log, .out, and .error files will be saved to a folder called `TestJobOutput`. We need to create this folder using `mkdir TestJobOutput` before we submit our job. 

We will call the script for this analysis `bwa-test.sh` and should contain the following information: 

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

Here we see that we used less than half of both the disk space and memory we requested. In future jobs, we could request a smaller amount of each resource, such as 0.5 GB of disk space and 0.5 GB of memory. We should run additional test jobs using these resource requests to ensure that they will allow our job to complete successfully before scaling up our analysis. 


# Scaling Up to Analyze Multiple Samples
In preparation for scaling up, please review our [guide on how to scale up after a successful test job](https://support.opensciencegrid.org/support/solutions/articles/12000076552-scaling-up-after-success-with-test-jobs) and how to 
[easily submit multiple jobs with a single submit file](https://support.opensciencegrid.org/support/solutions/articles/12000073165-easily-submit-multiple-jobs)

After reviewing how to submit multiple jobs with a single submit file, we see that an option for scaling up our analysis is to use `queue <var> from <list.txt>`. 

To use this option, we first need to create a file with just the names of our samples that we want to analyze and to cut all information after the "_" symbol to remove the forward/reverse read information and file extensions. We will save the sample names in a file called `samples.txt`. 

```
cd data/trimmed_fastq_small/
ls *.fastq | cut -f 1 -d '_' | uniq > samples.txt
cd ~
```

Now, we can create a new submit file called `bwa-alignment.sub` to queue a new job for each sample. 

```
universe        = vanilla     
executable      = bwa-alignment.sh

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
requirements = (OSGVO_OS_STRING == "RHEL 7")
+JobDurationCategory = "Medium"

transfer_input_files = software/bwa.tar.gz, data/ref_genome/ecoli_rel606.fasta.gz, data/trimmed_fastq_small/$(sample)_1.trim.sub.fastq, data/trimmed_fastq_small/$(sample)_2.trim.sub.fastq

transfer_output_remaps = "SRR2584863.aligned.sam=results/SRR2584863.aligned.sam; SRR2584866.aligned.sam=results/SRR2584866.aligned.sam; SRR2589044.aligned.sam=results/SRR2589044.aligned.sam"

"output_small.txt=outputs/output_small.txt.$(Cluster).$(Process)"

arguments = $(sample)

log         = log/bwa_$(sample)_job.log
output      = output/bwa_$(sample)_job.out
error       = error/bwa_$(sample)_job.error

request_cpus    = 1 
request_memory  = 0.5GB
request_disk    = 0.5GB

queue sample from data/trimmed_fastq_small/samples.txt
```

In addition to restructuring our submit file to queue a new job for each sample listed in our submit file, it is also advantageous to have our standard output, log, and error files saved to dedicated folders called "log", "output", and "error".  Therefore, we need to make these folders in our /home directory prior to submitting our job: 

```
mkdir log
mkdir output
mkdir error
mkdir results
```

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

When our jobs have completed, we can type confirm that our alignment output results files were created by typing:

```
cd results
ls -lah *
``` 

We can also investigate our log, error, and output files in their respective folders. 



For more information about running bioinformatics workflows on the OSG, we recommend our [BLAST tutorial](https://support.opensciencegrid.org/support/solutions/articles/12000062020-running-a-blast-workflow) as well as our [Samtools](https://support.opensciencegrid.org/support/solutions/articles/12000074984-example-software-compilation) instillation guide. 


