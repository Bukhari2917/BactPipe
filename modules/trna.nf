process TRNA {
    tag "${sample_id}"
    publishDir "${params.outdir}/genes", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "trna.out"
    
    script:
    """
    tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo "No tRNA found" > trna.out
    """
}