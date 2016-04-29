# Lab Notebook for TARA project

## 4/21/16
    * downloaded all the predicted genes for each station into "PredictedGenes" 
    folder in lustre scr space with DownloadORFs.sh
        * run: bsub -q day sh DownloadORFs.sh
        

## 4/23/16

Going to map and count all the reads to the OM-RGC 

* Download OM-RGC reference
    curl -O ftp://ftp.sra.ebi.ac.uk/vol1/ERA412/ERA412970/tab/OM-RGC_seq.release.tsv.gz
    (size 7.8GB compressed, make sure you ask for enough memory to map on the cluster)

Need download all the raw runs. It looks like for each station, if they ran multiple runs 
for that sample, they pooled all the runs per sample. I think. 

split DownloadRunsProk.sh into 10 files that will download stuff in chunks; submit 10 jobs for each file
 ftp DownloadRunsProk1-10.sh into TARA folder
 
* run: `bsub -q day sh DownloadRunsProk1.sh` for each; 
    * submitted the jobs at 3pm on Sunday, all but the last one (which had 43 files) ran in that time
    * had to run DownloadRunsProkPatch.sh to download the last 9 files that hadn't downloaded in Prok6
    
## Week of 4/25/16

Quality and adapter trimming of all the reads

wrote a shell script TrimReads*.sh to do trimming with trimmomatic
Broke it into six + 1 files TrimReads1.sh, TrimReads2.sh... TrimReads6.sh, then TrimReadsPatch.sh
(TrimReadsPatch.sh is last 7 or so runs that hadn't finished downloading at time of submission, submitted job next day) 

Used Phred score of 5, dropped anything below 45bp, and used TruSeq3 adapters (which is supposed to be right for HiSeq Illumina)

In the scripts, used 8 threads to run trimmomatic, so make sure you reserve 8 CPUs 

e.g. `bsub -q day -n 8 -R "span[hosts=1]" sh TrimReads5.sh`

## 4/26/16

add the bowtie2 module to my killdevil session `module add bowtie2`

build an index from the OM-RGC fasta file :


## Indexing OM-RGC

This was kind of a nightmare - it is a fragmented 40-some-million-gene reference catalog, 
at a total of 25G. 

I tried it with increasing amounts of memory on the campus cluster, but just couldn't do it with available resources.

I spun up an AWS EC2 on an r3.8xlarge machine with 244GB of memory, and attached 500GB storage (100GB was not sufficient, for some reason). 
This worked (although bowtie2 had to auto-adjust parameters to get it to work even then!)
And it took a looong time to figure out - both installing bowtie2, and errors in actually running bowtie2-build. 
AND my ssh connection timed out mid-run a couple of times, necessitating I updated my ServerAliveInterval and ServerAliveCountMax in my /.ssh/config file to 60 and 10, respectively. (DON'T misspell those lines! Trust me, don't).

remember to mount your storage: (make sure the 'xvdf' is correct, it might look like xvda, xvdx, etc)
`mkdir /mnt`
`mkfs -t ext4 /dev/xvdf`
`mount /dev/xvdf /mnt`

I copied the OM-RGC fasta file to my amazon instance from my killdevil terminal using scp
see here: http://angus.readthedocs.io/en/2014/amazon/transfer-files-between-instance.html


### Downloading bowtie2:

Make sure you download the linux software

`curl -O -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip`

unzip it:
`unzip bowtie2-2.2.9-linux-x86_64.zip`

copy all the bowtie2 executables to your path:
`cd bowtie2-2.2.9/`
`sudo cp bowtie2* /usr/local/bin`

move into your directory with the OM-RGC fasta file, and run bowtie2-build.

### Running bowtie2-build

In the end, this worked: (these were parameters bowtie2 chose automatically, when I reran the command after timing out and restarting that I used when I restarted the build command)

`bowtie2-build -f -a --bmax 1977633598 --dcv 1024 --threads 14 OM-RGC_seq.release.fna OM-RGC_idx`

it took about 6 hours to index it total on the EC2 machine in the end; 3 hours for the forward index and 3 hours for the mirror index

    