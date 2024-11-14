#!/bin/bash

# Sample slurm submission script for the Chimera compute cluster
# Lines beginning with # are comments, and will be ignored by
# the interpreter.Lines beginning with #SBATCH are directives
# to the scheduler.These in turn can be commented out by
# adding a second # (e.g. ##SBATCH lines will not be processed
# by the scheduler).
#
#
# set name of job
#SBATCH --job-name=ohmer 16S

# set the number of processors/tasks needed
#SBATCH -n 24

#set an account to use
#if not used then default will be used
# for scavenger users, use this format:
#BATCH --account=patrick.kearns
# for contributing users, use this format:
##SBATCH --account=

# set max wallclock timeDD-HH:MM:SS

# the default time will be 1 hour if not set
#SBATCH --time=00-24:00:00

# set a memory request
#SBATCH --mem=64gb

# Set filenames for stdout and stderr.%j can be used for the jobid.
# see "filename patterns" section of the sbatch man page for
# additional options
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#

# set the partition where the job will run.Multiple partitions can
# be specified as a comma separated list
# Use command "sinfo" to get the list of partitions
#SBATCH --partition=DGXA100
##SBATCH --partition=Intel6240,Intel6248,DGXA100

#When submitting to the GPU node, these following three lines are needed:

##SBATCH --gres=gpu:1
##SBATCH --export=NONE
#source /etc/profile
 

#Optional
# mail alert at start, end and/or failure of execution
# see the sbatch man page for other options
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=patrick.kearns@umb.edu

# Put your job commands here, including loading any needed
# modules or diagnostic echos.
module load anaconda3-2021.05
source activate qiime2-2023.5

## Processing of microbiome samples
cd /hpcstor6/scratch01/p/patrick.kearns/10_26_24_Ohmer-lab

#load raw FASTQ reads into QIIME
qiime tools import --type EMPPairedEndSequences --input-path ./data --output-path ohmer_seqs.qza

#demultiplex reads
qiime demux emp-paired \
  --i-seqs ohmer_seqs.qza \
 --m-barcodes-file ohmer_lab_map.txt \
 --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences ohmer_demux.qza \
  --o-error-correction-details  ohmer_demux-details.qza \
  --p-no-golay-error-correction 

 
 #export fastq files
qiime tools export --input-path ohmer_demux.qza --output-path ohmer_demux
 
 
