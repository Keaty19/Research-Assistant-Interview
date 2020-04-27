#!/usr/bin/env nextflow
 
/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads="/scratch/c.c1866917/NextFlow/input/ATACseq/*_{1,2}.fastq.gz"

/*
 * Set the locations of the directories as variables for use in
 * pathing to files
 */
baseDir="/scratch/c.c1866917/NextFlow/bin"
outputDir="/scratch/c.c1866917/NextFlow/output"
resourcesDir="/scratch/c.c1866917/NextFlow/resources"
inputDir="/scratch/c.c1866917/NextFlow/input"
TRIMMOMATIC="/apps/genomics/trimmomatic/0.39/trimmomatic-0.39.jar"
picard="/apps/genomics/picard/2.20.2/bin/picard.jar"
contam=file("/scratch/c.c1866917/NextFlow/resources/contaminants.fa")

/*
 * Set the reference genome file name and location
 */

genome_file=Channel.fromPath( "/scratch/c.c1866917/NextFlow/resources/GRCm38_sm.fa" )
gen_index="/scratch/c.c1866917/NextFlow/resources/GRCm38_sm.fa"

Channel
    .fromFilePairs( params.reads, size: 2, flat: true)                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .into { read_pairs; read_pairs2 }

/*
 * Step 0. Builds the genome index required by the mapping process
 */
process buildIndex {

    executor "slurm"
    cpus 8
    memory "64 GB"
    clusterOptions="-A scw1448"

    publishDir "${baseDir}/../resources", mode: "copy"

    input:
    file genome from genome_file
      
    output:
    file "${genome}" into genome_index

    script:    
    """
    module purge
    module load bwa

    bwa index ${genome}
    """
}

/*
* Step 1. Running the FastQC analysis on the paired-end fastq files
*/
process fastqcOriginal {

    executor "slurm" 
    cpus 8
    memory "64 GB"
    clusterOptions="-A scw1448"

    publishDir "$baseDir/../output/FastQC/aa047/Original", mode: "copy"

    input:
    set pair_id, file(fwd_reads), file(rev_reads) from read_pairs

    output:
    file "*_fastqc.{zip,html}" into fastqc_original

    script:
    """     
    module purge
    module load java
    module load FastQC

    fastqc -f fastq -q ${fwd_reads}
    fastqc -f fastq -q ${rev_reads}
    """
}

/*
* Step 2. Running Trimmomatic on the original fastq files
*/
process trimmomatic {
    executor "slurm"
    cpus 8
    errorStrategy {task.attempt < 4 ? 'retry' : 'ignore'}
    memory { 35.GB * task.attempt }
    runTime="24:00:00"
    clusterOptions="-A scw1448"

    publishDir "${baseDir}/../output/Trimmomatic/aa047", mode: "copy"

    input:
    set pair_id, file(fwd_reads), file(rev_reads) from read_pairs2

    output:
    set pair_id, "${fwd_reads.getSimpleName()}_trimmed.fastq.gz", "${rev_reads.getSimpleName()}_trimmed.fastq.gz" into trimmed_fastqc, bwa_align
    file "*_unpaired.fastq.gz" into trimmed_unpaired

    script:    
    """
    module purge
    module load java
    module load trimmomatic

    java -jar ${TRIMMOMATIC} PE -phred33 \
    ${fwd_reads} ${rev_reads} \
    ${fwd_reads.getSimpleName()}_trimmed.fastq.gz ${fwd_reads.getSimpleName()}_unpaired.fastq.gz \
    ${rev_reads.getSimpleName()}_trimmed.fastq.gz ${rev_reads.getSimpleName()}_unpaired.fastq.gz \
    ILLUMINACLIP:${contam}:2:30:10 HEADCROP:16 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:55
    """
}

/*
* Step 3. Running the FastQC analysis on the trimmed fastq files
*/
process fastqcTrimmed {

    executor "slurm" 
    cpus 8
    memory "64 GB"
    clusterOptions="-A scw1448"

    publishDir "$baseDir/../output/FastQC/aa047/Trimmed", mode: "copy"

    input:
    set trimmed_id, file(fwd_reads), file(rev_reads) from trimmed_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_trimmed

    script:
    """     
    module purge
    module load java
    module load FastQC

    fastqc -f fastq -q ${fwd_reads}
    fastqc -f fastq -q ${rev_reads}
    """
}

/*
* Step 4. Running bwa alignment on the trimmed read against the reference genome
*/
process bwaAlign {

    executor "slurm"
    cpus 10
    errorStrategy {task.attempt < 5 ? 'retry' : 'ignore'}
    memory { 45.GB * task.attempt }
    runTime="24:00:00"
    clusterOptions="-A scw1448"
     
    publishDir "${baseDir}/../output/bwa/aa047", mode: "copy"

    input:
    set trimmed_id, file(fwd_reads), file(rev_reads) from bwa_align

    output:
    file "*.bam" into alignment_file

    script:
    """
    module purge
    module load bwa
    module load samtools

    bwa mem -t 8 ${gen_index} ${fwd_reads} ${rev_reads} | samtools view -Shu - > ${fwd_reads.getSimpleName()}.bam
    """
}

/*
* Step 5. Sort the alignment bam file
*/
process sortAlign {

    executor "slurm"
    cpus 8
    errorStrategy {task.attempt < 5 ? 'retry' : 'ignore'}
    memory { 45.GB * task.attempt }
    runTime="24:00:00"
    clusterOptions="-A scw1448"
     
    publishDir "${baseDir}/../output/bwa/aa047", mode: "copy"

    input:
    file alignment from alignment_file

    output:
    file "*.sorted.bam" into sorted_file

    script:
    """
    module purge
    module load samtools

    samtools sort -m 80G ${alignment} > ${alignment.getSimpleName()}.sorted.bam
    """
}

/*
* Step 6. Marking and removing duplicates on the sorted bam file.
*/
process rmDuplicates {

    executor "slurm"
    cpus 8
    memory "64 GB"
    clusterOptions="-A scw1448"
     
    publishDir "${baseDir}/../output/bwa/aa047", mode: "copy"

    input:
    file sorted from sorted_file

    output:
    file "*.rmdup.bam" into index_file, macs2_bam
    file "*.rmduplicates.txt" into rmdup_log

    script:
    """
    module purge
    module load samtools
    module load java
    module load picard

    java -jar ${picard} MarkDuplicates INPUT=${sorted} OUTPUT=${sorted.getSimpleName()}.rmdup.bam METRICS_FILE=${sorted.getSimpleName()}.rmduplicates.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
    """
}

/*
* Step 7. Index the newly created rmdup file.
*/
process samIndex {

    executor "slurm"
    cpus 8
    memory "64 GB"
    clusterOptions="-A scw1448"
     
    publishDir "${baseDir}/../output/bwa/aa047", mode: "copy"

    input:
    file index from index_file

    output:
    file "*.bam.bai" into alignment_index

    script:
    """
    module purge
    module load samtools

    samtools index -b ${index} > ${index.baseName}.bam.bai
    """
}

/*
* Step 8. Carry out peak calling, with 3 different parameters, using macs2
*/
process macs2Peak {
    
    executor "slurm"
    cpus 8
    memory "64GB"
    clusterOptions="-A scw1448"
    
    publishDir "${baseDir}/../output/macs2/aa047", mode: "copy"

    input:
    file aligned from macs2_bam

    output:
    file "*.macs2_control_lambda.bdg" into macs2_conbdg
    file "*.macs2_model.r" into macs2_modelr
    file "*.macs2_peaks.narrowPeak" into macs2_peaksnar
    file "*.macs2_peaks.xls" into macs2_peaksxls
    file "*.macs2_summits.bed" into macs2_sumbed
    file "*.macs2_treat_pileup.bdg" into macs2_treatpilebdg

    script:
    organism_id="mm"
    """
    module purge
    module load macs2/2.1.2

    macs2 callpeak -t ${aligned} -f AUTO -g ${organism_id} -n ${aligned.getSimpleName()}.q0.01.macs2 -B -q 0.01 --shift -100 --extsize 200
    macs2 callpeak -t ${aligned} -f AUTO -g ${organism_id} -n ${aligned.getSimpleName()}.q0.05.macs2 -B -q 0.05 --shift -100 --extsize 200
    macs2 callpeak -t ${aligned} -f AUTO -g ${organism_id} -n ${aligned.getSimpleName()}.q0.1.macs2 -B -q 0.1 --shift -100 --extsize 200
    """
}