process BUILD_PRIMER_DB {
    tag "primer_db"
    label 'process_single'
    publishDir "results/primer_metrics", mode: 'copy', overwrite: true

    container 'docker.io/rasmuskriis/nextclade-python'
    containerOptions = "-v ${baseDir}/bin:/project-bin"

    input:
    path primer_bed
    path primer_fasta
    val primer_set_name

    output:
    path("${primer_set_name}.primer_db.json"), emit: primer_db

    script:
    """
    python3 /project-bin/primer_metrics.py build-db \
        --primer-bed ${primer_bed} \
        --primer-fasta ${primer_fasta} \
        --primer-set-name ${primer_set_name} \
        --output ${primer_set_name}.primer_db.json
    """
}
