process QUAST {
    tag "${sample_id}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "quast_results"
    
    script:
    """
    quast.py ${assembly} -o quast_results --threads ${task.cpus}
    """
}