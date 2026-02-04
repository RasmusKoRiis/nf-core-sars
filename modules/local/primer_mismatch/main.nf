process PRIMER_MISMATCH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir ?: 'results'}/primer_metrics", mode: 'copy', overwrite: true

    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    tuple val(meta), path(consensus)
    path primer_db
    val run_id

    output:
    tuple val(meta), path("${meta.id}_primer_mismatches.csv"), emit: primer_mismatches

    script:
    """
    python3 /project-bin/primer_metrics.py mismatch \
        --consensus ${consensus} \
        --primer-db ${primer_db} \
        --sample-id ${meta.id} \
        --run-id ${run_id} \
        --output ${meta.id}_primer_mismatches.csv
    """
}
