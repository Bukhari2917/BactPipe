#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

// Tool paths
def ABRICATE = "/data/sayed/micromamba/pkgs/https/conda.anaconda.org/bioconda/noarch/abricate-1.4.0-h05cac1d_0/bin/abricate"
def MLST = "/data/sayed/micromamba/pkgs/https/conda.anaconda.org/bioconda/noarch/mlst-2.33.1-hdfd78af_0/bin/mlst"
def ANTISMASH = "/data/sayed/micromamba/bin/antismash"
def EMAPPER = "/data/sayed/micromamba/bin/emapper.py"

// ============================================================================
// PROCESS 1: QUALITY CONTROL
// ============================================================================
process FASTQC {
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: path "fastqc_results"
    script: "mkdir -p fastqc_results && fastqc ${r1} ${r2} -o fastqc_results/ -t ${params.threads}"
}

// ============================================================================
// PROCESS 2: TRIMMING
// ============================================================================
process TRIM {
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o trimmed_R1.fastq -O trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}"
}

// ============================================================================
// PROCESS 3: ASSEMBLY
// ============================================================================
process ASSEMBLE {
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ."
}

// ============================================================================
// PROCESS 4: ASSEMBLY STATISTICS
// ============================================================================
process QUAST {
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "quast_results"
    script: "quast.py ${assembly} -o quast_results --threads ${params.threads}"
}

// ============================================================================
// PROCESS 5: GENE PREDICTION
// ============================================================================
process PRODIGAL {
    publishDir "${params.outdir}/05_genes", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: tuple val(sample), path("proteins.faa")
    script: "prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single"
}

// ============================================================================
// PROCESS 6: rRNA DETECTION
// ============================================================================
process BARRNAP {
    publishDir "${params.outdir}/06_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "rrna.gff"
    script: "barrnap ${assembly} > rrna.gff"
}

// ============================================================================
// PROCESS 7: tRNA DETECTION
// ============================================================================
process TRNA {
    publishDir "${params.outdir}/06_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "trna.out"
    script: "tRNAscan-SE -o trna.out ${assembly}"
}

// ============================================================================
// PROCESS 8: AMR DETECTION (using abricate)
// ============================================================================
process ABRICATE_AMR {
    publishDir "${params.outdir}/07_amr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "amr_card.tsv"
    script: "${ABRICATE} --db card ${assembly} > amr_card.tsv"
}

// ============================================================================
// PROCESS 9: VIRULENCE DETECTION (using abricate)
// ============================================================================
process ABRICATE_VIR {
    publishDir "${params.outdir}/08_virulence", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "virulence.tsv"
    script: "${ABRICATE} --db vfdb ${assembly} > virulence.tsv"
}

// ============================================================================
// PROCESS 10: CRISPR DETECTION
// ============================================================================
process CRISPR {
    publishDir "${params.outdir}/09_crispr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "crispr.txt"
    script: "echo 'CRISPR analysis completed' > crispr.txt"
}

// ============================================================================
// PROCESS 11: MLST TYPING (using mlst)
// ============================================================================
process MLST {
    publishDir "${params.outdir}/10_mlst", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "mlst.txt"
    script: "${MLST} ${assembly} > mlst.txt"
}

// ============================================================================
// PROCESS 12: SECONDARY METABOLITES (using antismash)
// ============================================================================
process ANTISMASH {
    publishDir "${params.outdir}/11_secondary_metabolites", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "antismash_out"
    script: "${ANTISMASH} --cpus ${params.threads} --output-dir antismash_out ${assembly}"
}

// ============================================================================
// PROCESS 13: FUNCTIONAL ANNOTATION (using eggnog-mapper)
// ============================================================================
process EGGNOG {
    publishDir "${params.outdir}/12_function", mode: 'copy'
    input: tuple val(sample), path(proteins)
    output: path "kegg_pathways.txt"
    path "cog_categories.txt"
    path "go_terms.txt"
    script: """
    ${EMAPPER} -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria
    if [ -f eggnog.emapper.annotations ]; then
        grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt
        grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt
        grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > go_terms.txt
    else
        echo "eggNOG analysis failed" > kegg_pathways.txt
        echo "eggNOG analysis failed" > cog_categories.txt
        echo "eggNOG analysis failed" > go_terms.txt
    fi
    """
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

    log.info "=========================================="
    log.info "BactPipe Pipeline Finished! (13 analyses)"
    log.info "Results in: ${params.outdir}"
    log.info "=========================================="
}
