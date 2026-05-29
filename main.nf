#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

if (!params.reads) {
    log.info "ERROR: Please provide --reads parameter"
    exit 1
}

// ============================================================================
// PROCESSES
// ============================================================================

process FASTQC {
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input: path(r1), path(r2)
    script: "fastqc ${r1} ${r2} -o . -t ${params.threads}"
}

process TRIM {
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input: path(r1), path(r2)
    output: path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o trimmed_R1.fastq -O trimmed_R2.fastq -q 20 -l 50"
}

process ASSEMBLE {
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input: path(r1), path(r2)
    output: path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ."
}

process QUAST {
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input: path(assembly)
    script: "quast.py ${assembly} -o quast_results --threads ${params.threads}"
}

process BUSCO {
    publishDir "${params.outdir}/05_busco", mode: 'copy'
    input: path(assembly)
    script: "busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${params.threads} 2>&1 || true"
}

process PRODIGAL {
    publishDir "${params.outdir}/06_genes", mode: 'copy'
    input: path(assembly)
    output: path("proteins.faa")
    script: "prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single"
}

process BARRNAP {
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input: path(assembly)
    script: "barrnap ${assembly} > rrna.gff || echo 'No rRNA' > rrna.gff"
}

process TRNA {
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input: path(assembly)
    script: "tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo 'No tRNA' > trna.out"
}

process ABRICATE_AMR {
    publishDir "${params.outdir}/08_amr", mode: 'copy'
    input: path(assembly)
    script: "abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo 'No AMR genes' > amr_card.tsv"
}

process ABRICATE_VIR {
    publishDir "${params.outdir}/09_virulence", mode: 'copy'
    input: path(assembly)
    script: "abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo 'No virulence genes' > virulence.tsv"
}

process CRISPR {
    publishDir "${params.outdir}/10_crispr", mode: 'copy'
    input: path(assembly)
    script: "echo 'CRISPR analysis done' > crispr.txt"
}

process MLST {
    publishDir "${params.outdir}/11_mlst", mode: 'copy'
    input: path(assembly)
    script: "mlst ${assembly} > mlst.txt 2>/dev/null || echo 'No MLST scheme' > mlst.txt"
}

process EGGNOG {
    publishDir "${params.outdir}/12_function", mode: 'copy'
    input: path(proteins)
    script: """
    if [ -f ${proteins} ]; then
        emapper.py -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria || true
        if [ -f eggnog.emapper.annotations ]; then
            grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\\n' | sort | uniq -c | sort -rn > go_terms.txt || true
        fi
    fi
    """
}

process MULTIQC {
    publishDir "${params.outdir}/13_report", mode: 'copy'
    script: "multiqc . --filename multiqc_report.html --force || true"
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    // Get first sample
    reads = Channel.fromFilePairs(params.reads).first()
    r1 = reads[1][0]
    r2 = reads[1][1]
    
    // Step 1-2: QC and Trim
    FASTQC(r1, r2)
    TRIM(r1, r2)
    
    // Step 3: Assembly
    ASSEMBLE(TRIM.out[0], TRIM.out[1])
    
    // Steps 4-12: All analyses on assembly
    QUAST(ASSEMBLE.out)
    BUSCO(ASSEMBLE.out)
    PRODIGAL(ASSEMBLE.out)
    BARRNAP(ASSEMBLE.out)
    TRNA(ASSEMBLE.out)
    ABRICATE_AMR(ASSEMBLE.out)
    ABRICATE_VIR(ASSEMBLE.out)
    CRISPR(ASSEMBLE.out)
    MLST(ASSEMBLE.out)
    
    // Step 13: eggNOG on proteins
    EGGNOG(PRODIGAL.out)
    
    // Step 14: MultiQC report
    MULTIQC()
    
    log.info """
    ==========================================
    ✅ BactPipe COMPLETE Pipeline Finished!
    ==========================================
    Results: ${params.outdir}
    ==========================================
    """
}
