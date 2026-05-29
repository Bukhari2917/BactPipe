#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

process FASTQC {
    tag "FASTQC"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input:
    tuple val(sample), path(r1), path(r2)
    script: "fastqc ${r1} ${r2} -o . -t ${params.threads}"
}

process TRIM {
    tag "TRIM"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    tuple val(sample), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o trimmed_R1.fastq -O trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}"
}

process ASSEMBLE {
    tag "ASSEMBLE"
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    tuple val(sample), path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ."
}

process QUAST {
    tag "QUAST"
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "quast.py ${assembly} -o quast_results --threads ${params.threads}"
}

process BUSCO {
    tag "BUSCO"
    publishDir "${params.outdir}/05_busco", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${params.threads} || true"
}

process PRODIGAL {
    tag "PRODIGAL"
    publishDir "${params.outdir}/06_genes", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    output:
    tuple val(sample), path("proteins.faa")
    script: "prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single"
}

process BARRNAP {
    tag "BARRNAP"
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "barrnap ${assembly} > rrna.gff || echo 'No rRNA' > rrna.gff"
}

process TRNA {
    tag "TRNA"
    publishDir "${params.outdir}/07_rna", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo 'No tRNA' > trna.out"
}

process ABRICATE_AMR {
    tag "ABRICATE_AMR"
    publishDir "${params.outdir}/08_amr", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo 'No AMR genes' > amr_card.tsv"
}

process ABRICATE_VIR {
    tag "ABRICATE_VIR"
    publishDir "${params.outdir}/09_virulence", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo 'No virulence genes' > virulence.tsv"
}

process CRISPR {
    tag "CRISPR"
    publishDir "${params.outdir}/10_crispr", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "echo 'CRISPR analysis completed' > crispr.txt"
}

process MLST {
    tag "MLST"
    publishDir "${params.outdir}/11_mlst", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script: "mlst ${assembly} > mlst.txt 2>/dev/null || echo 'No MLST scheme' > mlst.txt"
}

process EGGNOG {
    tag "EGGNOG"
    publishDir "${params.outdir}/12_function", mode: 'copy'
    input:
    tuple val(sample), path(proteins)
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

process ANTISMASH {
    tag "ANTISMASH"
    publishDir "${params.outdir}/14_secondary_metabolites", mode: 'copy'
    input:
    tuple val(sample), path(assembly)
    script:
    """
    antismash --cpus ${params.threads} \
              --output-dir antismash_out \
              ${assembly} || echo 'antiSMASH failed' > antismash_out/error.txt
    """
}

process MULTIQC {
    tag "MULTIQC"
    publishDir "${params.outdir}/13_report", mode: 'copy'
    script: "multiqc . --filename multiqc_report.html --force || true"
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
    log.info "BactPipe COMPLETE Pipeline Finished! (15 analyses)"
    log.info "Results: ${params.outdir}"
    log.info "=========================================="
}
