process REPORT {
    label 'process_single'



    //conda "bioconda::blast=2.15.0"
    container 'docker.io/rasmuskriis/blast_python_pandas:amd64'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory

    input:
    path(nextclade_stats)
    path(nextclade_mutations)
    path(resistance_mutations)
    val runid
    val seq_instrument
    path(samplesheet)
    
    output:
    path("${runid}.csv"), emit: report



    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    python /project-bin/report.py ${samplesheet}

    #Add constant parameters to the report
    # Add RunID column
    awk -v runid=${runid} -v OFS=',' '{ if (NR == 1) { print  \$0, "RunID" } else { print  \$0, runid } }' merged_report.csv > ${runid}_temp1.csv

    # Add Instrument ID column
    awk -v seq_instrument=${seq_instrument} -v OFS=',' '{ if (NR == 1) { print  \$0, "Instrument ID" } else { print  \$0, seq_instrument } }' ${runid}_temp1.csv > ${runid}_temp2.csv




    # Rename the final file to runID
    mv ${runid}_temp2.csv ${runid}.csv

    
    """

}
