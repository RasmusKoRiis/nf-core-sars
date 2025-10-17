process BCFTOOLS_CONSENSUS {
    tag { "${meta.id}" }
    publishDir "results/consensus", mode: 'copy', overwrite: true
    //errorStrategy 'ignore'
    
    container 'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509'
    input:
    tuple val(meta), path(vcf)
    path  reference
    tuple val(meta2), path(mask_bed)
    output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), path("${meta.id}.consensus.report.txt")

    script:
    """
    # Apply variants onto reference
    bcftools consensus -f ${reference} ${vcf} > ${meta.id}.tmp.fa

    # Mask low-depth / primer regions with Ns
    if [ -s ${mask_bed} ]; then
      bcftools consensus -f ${meta.id}.tmp.fa -m ${mask_bed} ${vcf} > ${meta.id}.consensus.fasta
    else
      mv ${meta.id}.tmp.fa ${meta.id}.consensus.fasta
    fi

    # Tiny report
    printf "sample\\t${meta.id}\\nvcf\\t${vcf}\\nmask\\t${mask_bed}\\n" > ${meta.id}.consensus.report.txt
    """
}