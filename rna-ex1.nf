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

