process BUSCO {
    tag "${sample_id}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "busco_results"
    
    script:
    """
    busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${task.cpus}
    """
}