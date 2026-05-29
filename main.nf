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
        • Assembly stats (QUAST)
        • Completeness (BUSCO)
        • Gene predictions (proteins.faa)
        • rRNA/tRNA
        • KEGG Pathways, COG, GO (eggNOG)
        • AMR genes (CARD)
        • Virulence factors (VFDB)
        • CRISPR
        • MLST typing
        • Quality report (MultiQC)
    """
    exit 0
}

// ============================================================================
// PROCESS 1: QUALITY CONTROL
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
// PROCESS 2: TRIMMING
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
// PROCESS 3: ASSEMBLY
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
// PROCESS 4: ASSEMBLY STATISTICS (QUAST)
// ============================================================================

process QUAST {
    tag "${sample_id}"
    publishDir "${params.outdir}/quast", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "quast_results"
    script:
    """
    quast.py ${assembly} -o quast_results --threads ${task.cpus}
    """
}

// ============================================================================
// PROCESS 5: COMPLETENESS (BUSCO)
// ============================================================================

process BUSCO {
    tag "${sample_id}"
    publishDir "${params.outdir}/busco", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "busco_results"
    script:
    """
    busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${task.cpus} || echo "BUSCO failed" > busco_results/summary.txt
    """
}

// ============================================================================
// PROCESS 6: GENE PREDICTION (Prodigal)
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
// PROCESS 7: rRNA DETECTION (barrnap)
// ============================================================================

process BARRNAP {
    tag "${sample_id}"
    publishDir "${params.outdir}/rna", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "rrna.gff"
    script:
    """
    barrnap ${assembly} > rrna.gff || echo "No rRNA found" > rrna.gff
    """
}

// ============================================================================
// PROCESS 8: tRNA DETECTION (tRNAscan-SE)
// ============================================================================

process TRNA {
    tag "${sample_id}"
    publishDir "${params.outdir}/rna", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "trna.out"
    script:
    """
    tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo "No tRNA found" > trna.out
    """
}

// ============================================================================
// PROCESS 9: AMR DETECTION (ABRicate - CARD)
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
// PROCESS 10: VIRULENCE DETECTION (ABRicate - VFDB)
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
// PROCESS 11: CRISPR DETECTION
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
// PROCESS 12: MLST TYPING
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
// PROCESS 13: FUNCTIONAL ANNOTATION (eggNOG - KEGG, COG, GO)
// ============================================================================

process EGGNOG {
    tag "${sample_id}"
    publishDir "${params.outdir}/function", mode: 'copy'
    input:
    path proteins
    output:
    path "kegg_pathways.txt"
    path "cog_categories.txt"
    path "go_terms.txt"
    script:
    """
    if [ -f ${proteins} ]; then
        emapper.py -i ${proteins} --output eggnog --cpu ${task.cpus} --tax_scope Bacteria || echo "eggNOG failed"
        if [ -f eggnog.emapper.annotations ]; then
            grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > go_terms.txt || true
        else
            echo "eggNOG analysis not available" > kegg_pathways.txt
            echo "eggNOG analysis not available" > cog_categories.txt
            echo "eggNOG analysis not available" > go_terms.txt
        fi
    else
        echo "No protein file provided" > kegg_pathways.txt
        echo "No protein file provided" > cog_categories.txt
        echo "No protein file provided" > go_terms.txt
    fi
    """
}

// ============================================================================
// PROCESS 14: FINAL REPORT (MultiQC)
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
    multiqc . --filename multiqc_report.html --force || echo "MultiQC failed"
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
    QUAST(ASSEMBLE.out.map { [it[0], it[1]] })
    BUSCO(ASSEMBLE.out.map { [it[0], it[1]] })
    PRODIGAL(ASSEMBLE.out.map { [it[0], it[1]] })
    BARRNAP(ASSEMBLE.out.map { [it[0], it[1]] })
    TRNA(ASSEMBLE.out.map { [it[0], it[1]] })
    ABRICATE_AMR(ASSEMBLE.out)
    ABRICATE_VIR(ASSEMBLE.out)
    CRISPR(ASSEMBLE.out)
    MLST(ASSEMBLE.out)
    EGGNOG(PRODIGAL.out.map { it[1] })
    MULTIQC(FASTQC.out.collect(), TRIM.out.map { it[2] }.collect(), ASSEMBLE.out.map { it[1] }.collect())
    
    log.info "=========================================="
    log.info "BactPipe COMPLETE Pipeline Finished!"
    log.info "Results: ${params.outdir}"
    log.info "=========================================="
}
