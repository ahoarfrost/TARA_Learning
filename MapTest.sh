runs="ERR598950
ERR599095"

forward="_1P.fq"
reverse="_2P.fq"

#"RawRuns/Trimmed/$run$forward"
#ERR315856_1P.fq
#ERR315856_2P.fq
#trimmed runs are in TARA/RawRuns/Trimmed

#bowtie indexFile fastqFile outputFile
#bowtie2 -p 12 -x reference/OM-RGC_idx -1 /data/drosophila/RAL357_1.fastq -2 /data/drosophila/RAL357_2.fastq -S RAL357_bowtie.sam

for run in $runs
do
    bowtie2 -x reference/OM-RGC_idx -1 RawRuns/Trimmed/$run$forward -2 RawRuns/Trimmed/$run$reverse -S $run_bowtie2map.sam
done
