/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/reads/*_gut_{1,2}.fq"
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