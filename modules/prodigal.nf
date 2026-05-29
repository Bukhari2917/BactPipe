process PRODIGAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/genes", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    tuple val(sample_id), path("proteins.faa")
    path "genes.gbk"
    
    script:
    """
    prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single
    """
}