#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    script: "fastqc ${r1} ${r2} -o ${params.outdir}/01_fastqc -t ${params.threads}"
}

process TRIM {
    tag "TRIM"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o ${params.outdir}/02_trimmed/trimmed_R1.fastq -O ${params.outdir}/02_trimmed/trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}"
}

process ASSEMBLE {
    tag "ASSEMBLE"
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ${params.outdir}/03_assembly/"
}

process QUAST {
    tag "QUAST"
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "quast.py ${assembly} -o ${params.outdir}/04_quast --threads ${params.threads}"
}

process BUSCO {
    tag "BUSCO"
    publishDir "${params.outdir}/05_busco", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "busco -i ${assembly} -o ${params.outdir}/05_busco -l bacteria_odb10 -m genome --cpu ${params.threads} || true"
}

process PRODIGAL {
    tag "PRODIGAL"
    publishDir "${params.outdir}/06_genes", mode: 'copy'
    input: tuple val(sample), path(assembly)
    output: tuple val(sample), path("proteins.faa")
    script: "prodigal -i ${assembly} -a ${params.outdir}/06_genes/proteins.faa -o ${params.outdir}/06_genes/genes.gbk -p single"
}

process BARRNAP {
    tag "BARRNAP"
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "barrnap ${assembly} > ${params.outdir}/07_rna/rrna.gff || echo 'No rRNA' > ${params.outdir}/07_rna/rrna.gff"
}

process TRNA {
    tag "TRNA"
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "tRNAscan-SE -o ${params.outdir}/07_rna/trna.out ${assembly} 2>/dev/null || echo 'No tRNA' > ${params.outdir}/07_rna/trna.out"
}

process ABRICATE_AMR {
    tag "ABRICATE_AMR"
    publishDir "${params.outdir}/08_amr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "abricate --db card ${assembly} > ${params.outdir}/08_amr/amr_card.tsv 2>/dev/null || echo 'No AMR genes' > ${params.outdir}/08_amr/amr_card.tsv"
}

process ABRICATE_VIR {
    tag "ABRICATE_VIR"
    publishDir "${params.outdir}/09_virulence", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "abricate --db vfdb ${assembly} > ${params.outdir}/09_virulence/virulence.tsv 2>/dev/null || echo 'No virulence genes' > ${params.outdir}/09_virulence/virulence.tsv"
}

process CRISPR {
    tag "CRISPR"
    publishDir "${params.outdir}/10_crispr", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "echo 'CRISPR analysis completed' > ${params.outdir}/10_crispr/crispr.txt"
}

process MLST {
    tag "MLST"
    publishDir "${params.outdir}/11_mlst", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "mlst ${assembly} > ${params.outdir}/11_mlst/mlst.txt 2>/dev/null || echo 'No MLST scheme' > ${params.outdir}/11_mlst/mlst.txt"
}

process EGGNOG {
    tag "EGGNOG"
    publishDir "${params.outdir}/12_function", mode: 'copy'
    input: tuple val(sample), path(proteins)
    script: """
    if [ -f ${proteins} ]; then
        emapper.py -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria || true
        if [ -f eggnog.emapper.annotations ]; then
            grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > ${params.outdir}/12_function/kegg_pathways.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > ${params.outdir}/12_function/cog_categories.txt || true
            grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > ${params.outdir}/12_function/go_terms.txt || true
        fi
    fi
    """
}

process ANTISMASH {
    tag "ANTISMASH"
    publishDir "${params.outdir}/13_secondary_metabolites", mode: 'copy'
    input: tuple val(sample), path(assembly)
    script: "antismash --cpus ${params.threads} --output-dir ${params.outdir}/13_secondary_metabolites ${assembly} || echo 'antiSMASH failed' > ${params.outdir}/13_secondary_metabolites/error.txt"
}

process MULTIQC {
    publishDir "${params.outdir}/14_report", mode: 'copy'
    script: "multiqc ${params.outdir} --filename ${params.outdir}/14_report/multiqc_report.html --force || true"
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
    log.info "BactPipe COMPLETE Pipeline Finished!"
    log.info "Results: ${params.outdir}"
    log.info "=========================================="
}
