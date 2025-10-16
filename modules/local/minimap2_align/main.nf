process MINIMAP2_ALIGN {
    tag { "${meta.id}" }
    publishDir "results/aln", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    container 'biocontainers/minimap2:v2.24--h5bf99c6_1'
    cpus 4
    memory '4 GB'
    input:
    tuple val(meta), path(reads)     // reads = List[fastq.gz]
    path  reference

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: minimap2

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads.join(' ')} \
      | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
    samtools index ${meta.id}.bam
    """
}










