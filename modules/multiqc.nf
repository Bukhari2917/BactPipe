process MULTIQC {
    publishDir "${params.outdir}/report", mode: 'copy'
    
    input:
    path fastqc
    path trim
    path assembly
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc . --filename multiqc_report.html --force
    """
}