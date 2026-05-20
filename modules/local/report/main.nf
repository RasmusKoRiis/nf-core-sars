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
    path(fasta)
    path(tool_versions, stageAs: 'tool_versions/version_??.yml')
    val version_control_metadata
    
    output:
    path("${runid}.csv"), emit: report
    path("${runid}_with_software_versions.csv"), emit: report_with_versions
    path("${runid}.fasta"), emit: report_fasta
    path("versions.yml"), emit: versions



    when:
    task.ext.when == null || task.ext.when


    script:
    """ 
    set -euo pipefail

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
    awk -v version="${release_version}" -v OFS=',' '{ if (NR == 1) { print \$0, "Release Version" } else { print \$0, version } }' ${runid}_temp4.csv > ${runid}_temp5.csv

    cat > version_control_static.txt <<'END_STATIC_VERSION_CONTROL_METADATA'
    ${version_control_metadata}
    END_STATIC_VERSION_CONTROL_METADATA

    python3 - version_control_static.txt tool_versions/* > version_control_metadata.txt <<'PY'
    import sys

    static_path = sys.argv[1]
    version_files = sys.argv[2:]

    parts = []
    with open(static_path) as fh:
        static_metadata = " ".join(line.strip() for line in fh if line.strip())
    if static_metadata:
        parts.append(static_metadata)

    process = None
    for path in version_files:
        try:
            lines = open(path)
        except OSError:
            continue
        with lines:
            for raw_line in lines:
                line = raw_line.rstrip()
                stripped = line.strip()
                if stripped.endswith(":") and line == stripped:
                    process = stripped[:-1].strip('"').strip("'")
                    continue
                if process and line[:1].isspace() and ":" in stripped:
                    tool, version = stripped.split(":", 1)
                    tool = tool.strip().strip('"').strip("'")
                    version = version.strip().strip('"').strip("'")
                    if tool and version:
                        parts.append(f"{process}.{tool}={version}")

    report_python = sys.version.split()[0]
    parts.append(f"REPORT.python={report_python}")
    print("; ".join(parts))
    PY

    # QC calculations for the primary report without software-version metadata
    python /project-bin/report_QC_calculation.py ${runid}_temp5.csv -o ${runid}.csv

    # Add software-version metadata to a separate full report
    python3 - ${runid}.csv version_control_metadata.txt ${runid}_with_software_versions.csv <<'PY'
    import csv
    import sys

    in_csv, metadata_path, out_csv = sys.argv[1:]
    with open(metadata_path) as fh:
        metadata = fh.read().strip()

    with open(in_csv, newline="") as in_fh, open(out_csv, "w", newline="") as out_fh:
        reader = csv.reader(in_fh)
        writer = csv.writer(out_fh)
        header = next(reader)
        writer.writerow(header + ["NGS_Script_vers"])
        for row in reader:
            writer.writerow(row + [metadata])
    PY

    # Make Multiple FASTA file for all samples
    cat ${fasta} > ${runid}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > ${runid}.csv
    Sample,RunID,Release Version,NGS_QC_Sum,GISAID_Comment,GISAID_Kommentar
    sample1,${runid},${release_version},,,
    EOF

    cat <<EOF > ${runid}_with_software_versions.csv
    Sample,RunID,Release Version,NGS_QC_Sum,GISAID_Comment,GISAID_Kommentar,NGS_Script_vers
    sample1,${runid},${release_version},,,,${version_control_metadata}; REPORT.python=stub
    EOF

    cat <<EOF > ${runid}.fasta
    >sample1
    ACGT
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    
    """

}
