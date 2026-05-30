#!/usr/bin/env nextflow

params.reads = "data/*_R{1,2}.fastq"
params.outdir = "results"

Channel
    .fromFilePairs(params.reads)
    .ifEmpty { error "No reads found" }
    .set { read_pairs }

process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    path "*.html"
    script:
    """
    fastqc ${reads} -o .
    """
}

process TRIMMING {
    tag "Trimming"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("*_R1_trimmed.fastq"), path("*_R2_trimmed.fastq")
    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \
        ${sample_id}_R1_trimmed.fastq ${sample_id}_R1_unpaired.fastq \
        ${sample_id}_R2_trimmed.fastq ${sample_id}_R2_unpaired.fastq \
        ILLUMINACLIP:adapters.fa:2:30:10 \
        SLIDINGWINDOW:4:5 MINLEN:50
    """
}

process ASSEMBLY {
    tag "Assembly"
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)
    output:
    path "contigs.fasta", emit: contigs
    script:
    """
    spades.py -1 ${reads_R1} -2 ${reads_R2} -o spades_out --isolate -t ${task.cpus}
    cp spades_out/contigs.fasta .
    """
}

process PRODIGAL {
    tag "Prodigal"
    publishDir "${params.outdir}/05_prodigal", mode: 'copy'
    input:
    path contigs
    output:
    path "proteins.faa"
    script:
    """
    prodigal -i ${contigs} -a proteins.faa -o genes.gff -p meta
    """
}

process ABRICATE_AMR {
    tag "AMR"
    publishDir "${params.outdir}/08_amr", mode: 'copy'
    input:
    path contigs
    output:
    path "amr_results.txt"
    script:
    """
    abricate ${contigs} --db resfinder > amr_results.txt
    """
}

process ABRICATE_VIRULENCE {
    tag "Virulence"
    publishDir "${params.outdir}/09_virulence", mode: 'copy'
    input:
    path contigs
    output:
    path "virulence_results.txt"
    script:
    """
    abricate ${contigs} --db vfdb > virulence_results.txt
    """
}

process CRISPR {
    tag "CRISPR"
    publishDir "${params.outdir}/10_crispr", mode: 'copy'
    input:
    path contigs
    output:
    path "crispr_results.txt"
    script:
    """
    minced ${contigs} > crispr_results.txt 2>&1
    """
}

process MLST {
    tag "MLST"
    publishDir "${params.outdir}/11_mlst", mode: 'copy'
    input:
    path contigs
    output:
    path "mlst_results.txt"
    script:
    """
    mlst ${contigs} > mlst_results.txt
    """
}

workflow {
    FASTQC(read_pairs)
    TRIMMING(read_pairs)
    TRIMMING.out.map { tuple(sample, r1, r2) }.set { trimmed }
    ASSEMBLY(trimmed)
    PRODIGAL(ASSEMBLY.out.contigs)
    ABRICATE_AMR(ASSEMBLY.out.contigs)
    ABRICATE_VIRULENCE(ASSEMBLY.out.contigs)
    CRISPR(ASSEMBLY.out.contigs)
    MLST(ASSEMBLY.out.contigs)
}
