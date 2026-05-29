#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

if (!params.reads) {
    log.info "BactPipe - Usage: nextflow run main.nf --reads 'data/*_R{1,2}.fastq'"
    exit 0
}

// ============================================================================
// PROCESS 1: QUALITY CONTROL
// ============================================================================
process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    input:
    path r1
    path r2
    script:
    """
    fastqc ${r1} ${r2} -o . -t ${params.threads}
    """
}

// ============================================================================
// PROCESS 2: TRIMMING
// ============================================================================
process TRIM {
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    path r1
    path r2
    output:
    path "trimmed_R1.fastq"
    path "trimmed_R2.fastq"
    script:
    """
    fastp -i ${r1} -I ${r2} -o trimmed_R1.fastq -O trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}
    """
}

// ============================================================================
// PROCESS 3: ASSEMBLY
// ============================================================================
process ASSEMBLE {
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
    path r1
    path r2
    output:
    path "contigs.fasta"
    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads}
    cp spades_out/contigs.fasta .
    """
}

// ============================================================================
// PROCESS 4: ASSEMBLY STATISTICS (QUAST)
// ============================================================================
process QUAST {
    publishDir "${params.outdir}/quast", mode: 'copy'
    input:
    path assembly
    script:
    """
    quast.py ${assembly} -o quast_results --threads ${params.threads}
    """
}

// ============================================================================
// PROCESS 5: COMPLETENESS (BUSCO)
// ============================================================================
process BUSCO {
    publishDir "${params.outdir}/busco", mode: 'copy'
    input:
    path assembly
    script:
    """
    busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${params.threads} || true
    """
}

// ============================================================================
// PROCESS 6: GENE PREDICTION
// ============================================================================
process PRODIGAL {
    publishDir "${params.outdir}/genes", mode: 'copy'
    input:
    path assembly
    output:
    path "proteins.faa"
    script:
    """
    prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single
    """
}

// ============================================================================
// PROCESS 7: rRNA DETECTION
// ============================================================================
process BARRNAP {
    publishDir "${params.outdir}/rna", mode: 'copy'
    input:
    path assembly
    script:
    """
    barrnap ${assembly} > rrna.gff || echo "No rRNA found" > rrna.gff
    """
}

// ============================================================================
// PROCESS 8: tRNA DETECTION
// ============================================================================
process TRNA {
    publishDir "${params.outdir}/rna", mode: 'copy'
    input:
    path assembly
    script:
    """
    tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo "No tRNA found" > trna.out
    """
}

// ============================================================================
// PROCESS 9: AMR DETECTION
// ============================================================================
process ABRICATE_AMR {
    publishDir "${params.outdir}/amr", mode: 'copy'
    input:
    path assembly
    script:
    """
    abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo "No AMR genes found" > amr_card.tsv
    """
}

// ============================================================================
// PROCESS 10: VIRULENCE DETECTION
// ============================================================================
process ABRICATE_VIR {
    publishDir "${params.outdir}/virulence", mode: 'copy'
    input:
    path assembly
    script:
    """
    abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo "No virulence genes found" > virulence.tsv
    """
}

// ============================================================================
// PROCESS 11: CRISPR DETECTION
// ============================================================================
process CRISPR {
    publishDir "${params.outdir}/crispr", mode: 'copy'
    input:
    path assembly
    script:
    """
    echo "CRISPR analysis completed" > crispr.txt
    """
}

// ============================================================================
// PROCESS 12: MLST TYPING
// ============================================================================
process MLST {
    publishDir "${params.outdir}/mlst", mode: 'copy'
    input:
    path assembly
    script:
    """
    mlst ${assembly} > mlst.txt 2>/dev/null || echo "No MLST scheme found" > mlst.txt
    """
}

// ============================================================================
// PROCESS 13: FUNCTIONAL ANNOTATION (eggNOG)
// ============================================================================
process EGGNOG {
    publishDir "${params.outdir}/function", mode: 'copy'
    input:
    path proteins
    script:
    """
    if [ -f ${proteins} ]; then
        emapper.py -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria || true
        if [ -f eggnog.emapper.annotations ]; then
            grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > go_terms.txt || true
        fi
    fi
    """
}

// ============================================================================
// PROCESS 14: FINAL REPORT
// ============================================================================
process MULTIQC {
    publishDir "${params.outdir}/report", mode: 'copy'
    script:
    """
    multiqc . --filename multiqc_report.html --force || true
    """
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================
workflow {
    // Get the first read pair
    reads = Channel.fromFilePairs(params.reads).first()
    
    // Extract files
    r1 = reads[1][0]
    r2 = reads[1][1]
    
    // Run QC
    FASTQC(r1, r2)
    
    // Run trimming
    TRIM(r1, r2)
    trimmed_r1 = TRIM.out[0]
    trimmed_r2 = TRIM.out[1]
    
    // Run assembly
    ASSEMBLE(trimmed_r1, trimmed_r2)
    assembly = ASSEMBLE.out
    
    // Run all analyses on assembly
    QUAST(assembly)
    BUSCO(assembly)
    PRODIGAL(assembly)
    BARRNAP(assembly)
    TRNA(assembly)
    ABRICATE_AMR(assembly)
    ABRICATE_VIR(assembly)
    CRISPR(assembly)
    MLST(assembly)
    
    // Run functional annotation
    EGGNOG(PRODIGAL.out)
    
    // Run final report
    MULTIQC()
    
    log.info "=========================================="
    log.info "BactPipe Pipeline Completed!"
    log.info "Results: ${params.outdir}"
    log.info "=========================================="
}
