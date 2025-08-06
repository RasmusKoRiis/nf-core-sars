process REPORTFASTA {
    label 'process_single'



    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    path(nextclade_stats)
    path(nextclade_mutations)
    path(resistance_mutations)
    val runid


    output:
    path("${runid}.csv"), emit: report




    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    # Generate date
    current_date=\$(date '+%Y-%m-%d')

    python /project-bin/report-fasta.py

    #Add constant parameters to the report
    # Add RunID column
    awk -v runid=${runid} -v OFS=',' '{ if (NR == 1) { print  \$0, "RunID" } else { print  \$0, runid } }' merged_report.csv > ${runid}_temp1.csv


    # Add Date column
    awk -v date="\$current_date" -v OFS=',' '{ if (NR == 1) { print \$0, "Date" } else { print \$0, date } }' ${runid}_temp1.csv > ${runid}.csv


    """

}
