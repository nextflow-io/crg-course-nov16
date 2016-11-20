/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/reads/ggal_gut_{1,2}.fq"
params.annot = "$baseDir/data/ggal/annotation.gff"
params.genome = "$baseDir/data/ggal/genome.fa"

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
    file genome from genome_file
     
    output:
    file 'genome.index*' into genome_index
       
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
read_pairs = Channel.fromFilePairs(params.reads, flat: true)

