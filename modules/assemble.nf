process ASSEMBLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("contigs.fasta")

    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${task.cpus} --rename
    
    # Fix contig names for Prokka
    awk '/^>/ {print ">contig_" ++i; next} {print}' spades_out/contigs.fasta > spades_out/contigs.fixed.fasta
    cp spades_out/contigs.fixed.fasta contigs.fasta
    """
}
