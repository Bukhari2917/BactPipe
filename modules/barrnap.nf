process BARRNAP {
    tag "${sample_id}"
    publishDir "${params.outdir}/genes", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "rrna.gff"
    
    script:
    """
    barrnap ${assembly} > rrna.gff
    """
}