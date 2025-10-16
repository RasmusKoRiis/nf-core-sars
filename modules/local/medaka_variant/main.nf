process MEDAKA_VARIANT {
    tag { "${meta.id}" }
    publishDir "results/medaka", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    container 'ontresearch/medaka:1.11.3'
    cpus 4
    memory '6 GB'

    input:
    tuple val(meta), path(bam)
    path  reference
    val   medaka_model

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: medaka_var

    script:
    """
    # Medaka wants reads; use bam via bam2fq then realign inside medaka_consensus? 
    # Simpler path: medaka_variant supports BAM with --regions whole reference via pileup
    # We'll extract FASTQ to be model-safe (small cost).
    samtools fastq ${bam} > ${meta.id}.tmp.fq

    medaka_variant \
      -i ${meta.id}.tmp.fq \
      -f ${reference} \
      -o medaka_${meta.id} \
      -m ${medaka_model} \
      -t ${task.cpus}

    bgzip -c medaka_${meta.id}/round_*/variants.vcf > ${meta.id}.vcf.gz
    tabix -p vcf ${meta.id}.vcf.gz
    rm -f ${meta.id}.tmp.fq
    """
}