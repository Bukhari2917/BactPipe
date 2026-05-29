process ASSEMBLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("contigs.fasta")
    
    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${task.cpus}
    cp spades_out/contigs.fasta .
    """
}