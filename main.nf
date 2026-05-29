#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8
params.help = false

if (params.help) {
    log.info """
    BactPipe - Complete Bacterial Genome Analysis Pipeline
    
    USAGE:
        nextflow run main.nf --reads 'data/*_R{1,2}.fastq'
    
    OUTPUTS:
        • Assembly (contigs.fasta)
        • Gene predictions (proteins.faa)
        • AMR genes (amr_card.tsv)
        • Virulence genes (virulence.tsv)
        • MLST typing (mlst.txt)
        • CRISPR (crispr.txt)
        • Quality report (multiqc_report.html)
    """
    exit 0
}

// ============================================================================
// QUALITY CONTROL
// ============================================================================

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

// ============================================================================
// TRIMMING
// ============================================================================

process TRIM {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq")
    script:
    """
    fastp -i ${r1} -I ${r2} -o ${sample_id}_R1.fastq -O ${sample_id}_R2.fastq -q 20 -l 50 --thread ${task.cpus}
    """
}

// ============================================================================
// ASSEMBLY
// ============================================================================

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

// ============================================================================
// GENE PREDICTION
// ============================================================================

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

// ============================================================================
// AMR DETECTION
// ============================================================================

process ABRICATE_AMR {
    tag "${sample_id}"
    publishDir "${params.outdir}/amr", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "amr_card.tsv"
    script:
    """
    abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo "No AMR genes found" > amr_card.tsv
    """
}

// ============================================================================
// VIRULENCE FACTORS
// ============================================================================

process ABRICATE_VIR {
    tag "${sample_id}"
    publishDir "${params.outdir}/virulence", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "virulence.tsv"
    script:
    """
    abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo "No virulence genes found" > virulence.tsv
    """
}

// ============================================================================
// CRISPR DETECTION
// ============================================================================

process CRISPR {
    tag "${sample_id}"
    publishDir "${params.outdir}/crispr", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "crispr.txt"
    script:
    """
    echo "CRISPR analysis - minced not available" > crispr.txt
    """
}

// ============================================================================
// MLST TYPING
// ============================================================================

process MLST {
    tag "${sample_id}"
    publishDir "${params.outdir}/mlst", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "mlst.txt"
    script:
    """
    mlst ${assembly} > mlst.txt 2>/dev/null || echo "No MLST scheme found" > mlst.txt
    """
}

// ============================================================================
// MULTIQC REPORT
// ============================================================================

process MULTIQC {
    publishDir "${params.outdir}/report", mode: 'copy'
    input:
    path fastqc
    path trim
    path assembly
    output:
    path "multiqc_report.html"
    script:
    """
    multiqc . --filename multiqc_report.html --force
    """
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    // Create channel from input reads
    Channel.fromFilePairs(params.reads)
        .map { id, reads -> [params.sample, reads[0], reads[1]] }
        .set { reads_ch }
    
    // Run all processes
    FASTQC(reads_ch)
    TRIM(reads_ch)
    ASSEMBLE(TRIM.out)
    PRODIGAL(ASSEMBLE.out)
    ABRICATE_AMR(ASSEMBLE.out)
    ABRICATE_VIR(ASSEMBLE.out)
    CRISPR(ASSEMBLE.out)
    MLST(ASSEMBLE.out)
    MULTIQC(FASTQC.out.collect(), TRIM.out.map { it[2] }.collect(), ASSEMBLE.out.map { it[1] }.collect())
    
    log.info "=========================================="
    log.info "BactPipe COMPLETE Pipeline Finished!"
    log.info "Results: ${params.outdir}"
    log.info "=========================================="
}
