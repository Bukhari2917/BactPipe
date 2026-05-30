#!/usr/bin/env nextflow

// BactPipe - Bacterial Genome Analysis Pipeline
// Analyses: QC, Trimming, Assembly, QUAST, Prodigal, rRNA, tRNA, AMR, Virulence, CRISPR, MLST

params.reads = "data/*_R{1,2}.fastq"
params.outdir = "results"

Channel
    .fromFilePairs(params.reads)
    .ifEmpty { error "No reads found matching pattern: ${params.reads}" }
    .set { read_pairs }

// 1. FASTQC
process FASTQC {
    tag "FASTQC: ${sample_id}"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    path "*.html"
    path "*.zip"
    script:
    """
    fastqc ${reads} -o .
    """
}

// 2. Trimming
process TRIMMING {
    tag "Trimming: ${sample_id}"
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

// 3. Assembly
process ASSEMBLY {
    tag "Assembly: ${sample_id}"
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

// 4. QUAST
process QUAST {
    tag "QUAST"
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input:
    path contigs
    output:
    path "quast_results/"
    script:
    """
    quast.py ${contigs} -o quast_results
    """
}

// 5. Prodigal
process PRODIGAL {
    tag "Prodigal"
    publishDir "${params.outdir}/05_prodigal", mode: 'copy'
    input:
    path contigs
    output:
    path "proteins.faa"
    path "genes.gff"
    script:
    """
    prodigal -i ${contigs} -a proteins.faa -o genes.gff -p meta
    """
}

// 6. rRNA
process RRNA {
    tag "rRNA"
    publishDir "${params.outdir}/06_rrna", mode: 'copy'
    input:
    path contigs
    output:
    path "rrna_results.txt"
    script:
    """
    barrnap ${contigs} > rrna_results.txt
    """
}

// 7. tRNA
process TRNA {
    tag "tRNA"
    publishDir "${params.outdir}/07_trna", mode: 'copy'
    input:
    path contigs
    output:
    path "trna_results.txt"
    script:
    """
    tRNAscan-SE ${contigs} -o trna_results.txt
    """
}

// 8. AMR Detection
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

// 9. Virulence Detection
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

// 10. CRISPR Detection
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

// 11. MLST
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

// Main Workflow
workflow {
    FASTQC(read_pairs)
    TRIMMING(read_pairs)
    TRIMMING.out.map { tuple(sample, r1, r2) }.set { trimmed_reads }
    ASSEMBLY(trimmed_reads)
    QUAST(ASSEMBLY.out.contigs)
    PRODIGAL(ASSEMBLY.out.contigs)
    RRNA(ASSEMBLY.out.contigs)
    TRNA(ASSEMBLY.out.contigs)
    ABRICATE_AMR(ASSEMBLY.out.contigs)
    ABRICATE_VIRULENCE(ASSEMBLY.out.contigs)
    CRISPR(ASSEMBLY.out.contigs)
    MLST(ASSEMBLY.out.contigs)
    log.info "=========================================="
    log.info "BactPipe Pipeline Finished!"
    log.info "Results in: ${params.outdir}"
    log.info "=========================================="
}
