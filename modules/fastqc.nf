process FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    path "fastqc_results"
    
    script:
    """
    mkdir -p fastqc_results
    fastqc ${r1} ${r2} -o fastqc_results/ -t ${task.cpus}
    """
}