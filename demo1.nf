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

