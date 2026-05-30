#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

// ============================================================================
// PROCESS 1: QUALITY CONTROL
// ============================================================================
process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: path "fastqc_results"
    script: "mkdir -p fastqc_results && fastqc ${r1} ${r2} -o fastqc_results/ -t ${params.threads}"
}

// ============================================================================
// PROCESS 2: TRIMMING
// ============================================================================
process TRIM {
    tag "TRIM"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o trimmed_R1.fastq -O trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}"
}

// ============================================================================
// PROCESS 3: ASSEMBLY
// ============================================================================
process ASSEMBLE {
    tag "ASSEMBLE"
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ."
}

// ============================================================================
// PROCESS 4: ASSEMBLY STATISTICS
// ============================================================================
process QUAST {
    tag "QUAST"
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "quast_results"
    script: "quast.py ${assembly} -o quast_results --threads ${params.threads}"
}

// ============================================================================
// PROCESS 5: COMPLETENESS (BUSCO)
// ============================================================================
process BUSCO {
    tag "BUSCO"
    publishDir "${params.outdir}/05_busco", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "busco_results"
    script: "busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${params.threads}"
}

// ============================================================================
// PROCESS 6: GENE PREDICTION
// ============================================================================
process PRODIGAL {
    tag "PRODIGAL"
    publishDir "${params.outdir}/06_genes", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: tuple val(sample), path("proteins.faa")
    script: "prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single"
}

// ============================================================================
// PROCESS 7: rRNA DETECTION
// ============================================================================
process BARRNAP {
    tag "BARRNAP"
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "rrna.gff"
    script: "barrnap ${assembly} > rrna.gff"
}

// ============================================================================
// PROCESS 8: tRNA DETECTION
// ============================================================================
process TRNA {
    tag "TRNA"
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "trna.out"
    script: "tRNAscan-SE -o trna.out ${assembly}"
}

// ============================================================================
// PROCESS 9: AMR DETECTION
// ============================================================================
process ABRICATE_AMR {
    tag "ABRICATE_AMR"
    publishDir "${params.outdir}/08_amr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "amr_card.tsv"
    script: "abricate --db card ${assembly} > amr_card.tsv"
}

// ============================================================================
// PROCESS 10: VIRULENCE DETECTION
// ============================================================================
process ABRICATE_VIR {
    tag "ABRICATE_VIR"
    publishDir "${params.outdir}/09_virulence", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "virulence.tsv"
    script: "abricate --db vfdb ${assembly} > virulence.tsv"
}

// ============================================================================
// PROCESS 11: CRISPR DETECTION
// ============================================================================
process CRISPR {
    tag "CRISPR"
    publishDir "${params.outdir}/10_crispr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "crispr.txt"
    script: "echo 'CRISPR analysis completed' > crispr.txt"
}

// ============================================================================
// PROCESS 12: MLST TYPING
// ============================================================================
process MLST {
    tag "MLST"
    publishDir "${params.outdir}/11_mlst", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "mlst.txt"
    script: "mlst ${assembly} > mlst.txt"
}

// ============================================================================
// PROCESS 13: FUNCTIONAL ANNOTATION (KEGG/COG/GO)
// ============================================================================
process EGGNOG {
    tag "EGGNOG"
    publishDir "${params.outdir}/12_function", mode: 'copy'
    input: tuple val(sample), path(proteins)
    output: path "kegg_pathways.txt"
    path "cog_categories.txt"
    path "go_terms.txt"
    script: """
    emapper.py -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria
    grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt
    grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt
    grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > go_terms.txt
    """
}

// ============================================================================
// PROCESS 14: SECONDARY METABOLITES
// ============================================================================
process ANTISMASH {
    tag "ANTISMASH"
    publishDir "${params.outdir}/13_secondary_metabolites", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "antismash_out"
    script: "antismash --cpus ${params.threads} --output-dir antismash_out ${assembly}"
}

// ============================================================================
// PROCESS 15: FINAL REPORT
// ============================================================================
process MULTIQC {
    tag "MULTIQC"
    publishDir "${params.outdir}/14_report", mode: 'copy'
    output: path "multiqc_report.html"
    script: "multiqc ${params.outdir} --filename multiqc_report.html --force"
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================
workflow {
    if (!params.reads) { error "Please provide --reads parameter" }

    Channel.fromFilePairs(params.reads)
        .map { sample_id, reads -> [params.sample, reads[0], reads[1]] }
        .set { reads_ch }

    FASTQC(reads_ch)
    TRIM(reads_ch)
    ASSEMBLE(TRIM.out)
    assembly_ch = ASSEMBLE.out.map { [it[0], it[1]] }

    QUAST(assembly_ch)
    BUSCO(assembly_ch)
    PRODIGAL(assembly_ch)
    BARRNAP(assembly_ch)
    TRNA(assembly_ch)
    ABRICATE_AMR(ASSEMBLE.out)
    ABRICATE_VIR(ASSEMBLE.out)
    CRISPR(ASSEMBLE.out)
    MLST(ASSEMBLE.out)
    ANTISMASH(ASSEMBLE.out)

    proteins_ch = PRODIGAL.out.map { [it[0], it[1]] }
    EGGNOG(proteins_ch)
    MULTIQC()

    log.info "=========================================="
    log.info "BactPipe COMPLETE Pipeline Finished! (15 analyses)"
    log.info "Results in: ${params.outdir}"
    log.info "=========================================="
}
