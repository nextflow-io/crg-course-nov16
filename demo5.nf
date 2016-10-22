/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/reads/ggal_gut_{1,2}.fq"
params.annot = "$baseDir/data/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

/* 
 * prints user convenience 
 */
println "R N A T O Y   P I P E L I N E    "
println "================================="
println "genome             : ${params.genome}"
println "annotat            : ${params.annot}"
println "reads              : ${params.reads}"


/*
 * get a file object for the given param string
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)
 
/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    
    input:
    file genome_file
     
    output:
    file 'genome.index*' into genome_index
       
    """
    bowtie2-build --threads ${task.cpus} ${genome_file} genome.index
    """
}

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
read_pairs = Channel.fromFilePairs(params.reads, flat: true)

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
     
    input:
    file 'genome.index.fa' from genome_file 
    file genome_index from genome_index.first()
    set pair_id, file(read1), file(read2) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam
 
    """
    tophat2 -p ${task.cpus} genome.index ${read1} ${read2}
    """
}

/*
 * Step 3. Assembles the transcript by using the "cufflinks" tool
 */
process makeTranscript {
       
    input:
    file annot from annotation_file 
    set pair_id, file(bam_file) from bam
     
    output:
    set pair_id, file('transcript_*.gtf') into transcripts
 
    """
    cufflinks --no-update-check -q -p ${task.cpus} -G $annot ${bam_file}
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}
