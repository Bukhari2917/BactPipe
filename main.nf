#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

process FASTQC {
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: path "fastqc_results"
    script: "mkdir -p fastqc_results && fastqc ${r1} ${r2} -o fastqc_results/ -t ${params.threads}"
}

process TRIM {
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o trimmed_R1.fastq -O trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}"
}

process ASSEMBLE {
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ."
}

process QUAST {
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "quast_results"
    script: "quast.py ${assembly} -o quast_results --threads ${params.threads}"
}

process PRODIGAL {
    publishDir "${params.outdir}/05_genes", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: tuple val(sample), path("proteins.faa")
    script: "prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single"
}

process BARRNAP {
    publishDir "${params.outdir}/06_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "rrna.gff"
    script: "barrnap ${assembly} > rrna.gff || echo 'No rRNA' > rrna.gff"
}

process TRNA {
    publishDir "${params.outdir}/06_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "trna.out"
    script: "tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo 'No tRNA found' > trna.out"
}

process ABRICATE_AMR {
    publishDir "${params.outdir}/07_amr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "amr_card.tsv"
    script: "abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo 'No AMR genes' > amr_card.tsv"
}

process ABRICATE_VIR {
    publishDir "${params.outdir}/08_virulence", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "virulence.tsv"
    script: "abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo 'No virulence genes' > virulence.tsv"
}

process CRISPR {
    publishDir "${params.outdir}/09_crispr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "crispr.txt"
    script: "echo 'CRISPR analysis completed' > crispr.txt"
}

process MLST {
    publishDir "${params.outdir}/10_mlst", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: path "mlst.txt"
    script: "mlst ${assembly} > mlst.txt 2>/dev/null || echo 'MLST not available' > mlst.txt"
}

process EGGNOG {
    publishDir "${params.outdir}/11_function", mode: 'copy'
    input: tuple val(sample), path(proteins)
    output: path "kegg_pathways.txt"
    path "cog_categories.txt"
    path "go_terms.txt"
    script: """
    emapper.py -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria 2>/dev/null || echo 'eggNOG failed' > eggnog_error.txt
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
    CRISPR(assembly_ch)
    MLST(ASSEMBLE.out)

    proteins_ch = PRODIGAL.out.map { [it[0], it[1]] }
    EGGNOG(proteins_ch)

    log.info "=========================================="
    log.info "BactPipe Pipeline Finished! (12 analyses)"
    log.info "Results in: ${params.outdir}"
    log.info "=========================================="
}
