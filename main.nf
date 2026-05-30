#!/usr/bin/env nextflow

// ============================================
// BactPipe - Complete Bacterial Analysis Pipeline
// 13 analyses: QC, Trimming, Assembly, QUAST, 
// Prodigal, rRNA, tRNA, AMR, Virulence, CRISPR, 
// MLST, antiSMASH, eggNOG
// ============================================

// Define parameters
params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"

// Create channel from read pairs
Channel
    .fromFilePairs(params.reads)
    .ifEmpty { error "No reads found matching pattern: ${params.reads}" }
    .set { read_pairs }

// ============================================
// 1. FASTQC - Quality control
// ============================================
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

// ============================================
// 2. Trimming - Remove adapters and low quality
// ============================================
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
        SLIDINGWINDOW:4:5 \
        MINLEN:20
    """
}

// ============================================
// 3. Assembly - SPAdes genome assembly
// ============================================
process ASSEMBLY {
    tag "Assembly: ${sample_id}"
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)
    
    output:
    path "${sample_id}/"
    path "contigs.fasta"
    
    script:
    """
    spades.py -1 ${reads_R1} -2 ${reads_R2} \
        -o ${sample_id} \
        --isolate \
        -t ${task.cpus}
    cp ${sample_id}/contigs.fasta .
    """
}

// ============================================
// 4. QUAST - Assembly evaluation
// ============================================
process QUAST {
    tag "QUAST: ${sample_id}"
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

// ============================================
// 5. Prodigal - Gene prediction
// ============================================
process PRODIGAL {
    tag "Prodigal: ${sample_id}"
    publishDir "${params.outdir}/05_prodigal", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "genes.gff"
    path "proteins.faa"
    path "nucleotides.fna"
    
    script:
    """
    prodigal -i ${contigs} -a proteins.faa -d nucleotides.fna -f gff -o genes.gff
    """
}

// ============================================
// 6. rRNA detection - Barrnap
// ============================================
process RRNA {
    tag "rRNA: ${sample_id}"
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

// ============================================
// 7. tRNA detection - tRNAscan-SE
// ============================================
process TRNA {
    tag "tRNA: ${sample_id}"
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

// ============================================
// 8. AMR detection - Abricate (ResFinder)
// ============================================
process ABRICATE_AMR {
    tag "AMR: ${sample_id}"
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

// ============================================
// 9. Virulence detection - Abricate (VFDB)
// ============================================
process ABRICATE_VIRULENCE {
    tag "Virulence: ${sample_id}"
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

// ============================================
// 10. CRISPR detection - MinCED
// ============================================
process CRISPR {
    tag "CRISPR: ${sample_id}"
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

// ============================================
// 11. MLST - Multi-locus sequence typing
// ============================================
process MLST {
    tag "MLST: ${sample_id}"
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

// ============================================
// 12. eggNOG - Functional annotation (KEGG/COG/GO)
// ============================================
process EGGNOG {
    tag "eggNOG: ${sample_id}"
    publishDir "${params.outdir}/12_eggnog", mode: 'copy'
    
    input:
    path proteins_faa
    
    output:
    path "eggnog_results/"
    
    script:
    """
    emapper.py -i ${proteins_faa} \
        --output eggnog_results \
        --cpu ${task.cpus} \
        --dmnd_db /data/databases/eggnog_proteins.dmnd \
        --data_dir /data/databases/eggnog_data || true
    """
}

// ============================================
// MAIN WORKFLOW - Connects all processes
// ============================================
workflow {
    // Step 1: Quality control
    FASTQC(read_pairs)
    
    // Step 2: Trim reads
    TRIMMING(read_pairs)
    TRIMMING.out.map { tuple(sample, r1, r2) }.set { trimmed_reads }
    
    // Step 3: Assemble genome
    ASSEMBLY(trimmed_reads)
    ASSEMBLY.out.contigs.set { assembly_contigs }
    
    // Step 4: Assembly evaluation
    QUAST(assembly_contigs)
    
    // Step 5: Gene prediction
    PRODIGAL(assembly_contigs)
    
    // Step 6-7: RNA detection
    RRNA(assembly_contigs)
    TRNA(assembly_contigs)
    
    // Step 8-9: AMR and Virulence
    ABRICATE_AMR(assembly_contigs)
    ABRICATE_VIRULENCE(assembly_contigs)
    
    // Step 10-11: CRISPR and MLST
    CRISPR(assembly_contigs)
    MLST(assembly_contigs)
    
    // Step 12: Functional annotation (optional, may fail if no database)
    EGGNOG(PRODIGAL.out.proteins_faa)
    
    // Completion message
    log.info "=========================================="
    log.info "BactPipe Pipeline Finished! (12 analyses)"
    log.info "Results in: ${params.outdir}"
    log.info "=========================================="
}
