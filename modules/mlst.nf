process MLST {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "mlst.txt"
    
    script:
    """
    mlst ${assembly} > mlst.txt 2>/dev/null || echo "No MLST scheme found" > mlst.txt
    """
}