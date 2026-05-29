process ABRICATE_VIR {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "virulence.tsv"
    
    script:
    """
    abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo "No virulence genes found" > virulence.tsv
    """
}