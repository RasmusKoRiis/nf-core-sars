process MINIMAP2_ALIGN {
    tag { "${meta.id}" }
    publishDir "results/aln", mode: 'copy', overwrite: true
    //errorStrategy 'ignore'

    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
            
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










