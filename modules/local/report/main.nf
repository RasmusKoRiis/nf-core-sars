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
    val release_version
    val seq_instrument
    path(samplesheet)
    val primer
    
    output:
    path("${runid}.csv"), emit: report



    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    # Generate date
    current_date=\$(date '+%Y-%m-%d')

    #turn csv int tsv
    sed 's/,/\t/g' ${samplesheet} > samplesheet.tsv

    python /project-bin/report.py samplesheet.tsv

    #Add constant parameters to the report
    # Add RunID column
    awk -v runid=${runid} -v OFS=',' '{ if (NR == 1) { print  \$0, "RunID" } else { print  \$0, runid } }' merged_report.csv > ${runid}_temp1.csv

    # Add Instrument ID column
    awk -v seq_instrument=${seq_instrument} -v OFS=',' '{ if (NR == 1) { print  \$0, "Instrument ID" } else { print  \$0, seq_instrument } }' ${runid}_temp1.csv > ${runid}_temp2.csv

    # Add primer scheme used
    awk -v primer=\$(basename ${primer}) -v OFS=',' '{
        if (NR == 1) { 
            print \$0, "Primer" 
        } else { 
            print \$0, primer 
        } 
    }' ${runid}_temp2.csv > ${runid}_temp3.csv

    # Add Date column
    awk -v date="\$current_date" -v OFS=',' '{ if (NR == 1) { print \$0, "Date" } else { print \$0, date } }' ${runid}_temp3.csv > ${runid}_temp4.csv

    # Add Release Version column
    awk -v version="${release_version}" -v OFS=',' '{ if (NR == 1) { print \$0, "Release Version" } else { print \$0, version } }' ${runid}_temp4.csv > ${runid}.csv



    
    """

}
