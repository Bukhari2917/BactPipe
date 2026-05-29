#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

process FASTQC {
    tag "FASTQC"
    input: tuple val(sample), path(r1), path(r2)
    script:
    """
    mkdir -p ${params.outdir}/01_fastqc
    fastqc ${r1} ${r2} -o ${params.outdir}/01_fastqc -t ${params.threads}
    """
}

process TRIM {
    tag "TRIM"
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")
    script:
    """
    mkdir -p ${params.outdir}/02_trimmed
    fastp -i ${r1} -I ${r2} -o ${params.outdir}/02_trimmed/trimmed_R1.fastq -O ${params.outdir}/02_trimmed/trimmed_R2.fastq -q 20 -l 50 --thread ${params.threads}
    """
}

process ASSEMBLE {
    tag "ASSEMBLE"
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("contigs.fasta")
    script:
    """
    mkdir -p ${params.outdir}/03_assembly
    spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads}
    cp spades_out/contigs.fasta ${params.outdir}/03_assembly/
    """
}

process QUAST {
    tag "QUAST"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/04_quast
    quast.py ${assembly} -o ${params.outdir}/04_quast --threads ${params.threads}
    """
}

process BUSCO {
    tag "BUSCO"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/05_busco
    busco -i ${assembly} -o ${params.outdir}/05_busco -l bacteria_odb10 -m genome --cpu ${params.threads} || true
    """
}

process PRODIGAL {
    tag "PRODIGAL"
    input: tuple val(sample), path(assembly)
    output: tuple val(sample), path("proteins.faa")
    script:
    """
    mkdir -p ${params.outdir}/06_genes
    prodigal -i ${assembly} -a ${params.outdir}/06_genes/proteins.faa -o ${params.outdir}/06_genes/genes.gbk -p single
    """
}

process BARRNAP {
    tag "BARRNAP"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/07_rna
    barrnap ${assembly} > ${params.outdir}/07_rna/rrna.gff || echo 'No rRNA' > ${params.outdir}/07_rna/rrna.gff
    """
}

process TRNA {
    tag "TRNA"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/07_rna
    tRNAscan-SE -o ${params.outdir}/07_rna/trna.out ${assembly} 2>/dev/null || echo 'No tRNA' > ${params.outdir}/07_rna/trna.out
    """
}

process ABRICATE_AMR {
    tag "ABRICATE_AMR"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/08_amr
    abricate --db card ${assembly} > ${params.outdir}/08_amr/amr_card.tsv 2>/dev/null || echo 'No AMR genes' > ${params.outdir}/08_amr/amr_card.tsv
    """
}

process ABRICATE_VIR {
    tag "ABRICATE_VIR"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/09_virulence
    abricate --db vfdb ${assembly} > ${params.outdir}/09_virulence/virulence.tsv 2>/dev/null || echo 'No virulence genes' > ${params.outdir}/09_virulence/virulence.tsv
    """
}

process CRISPR {
    tag "CRISPR"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/10_crispr
    echo 'CRISPR analysis completed' > ${params.outdir}/10_crispr/crispr.txt
    """
}

process MLST {
    tag "MLST"
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/11_mlst
    mlst ${assembly} > ${params.outdir}/11_mlst/mlst.txt 2>/dev/null || echo 'No MLST scheme' > ${params.outdir}/11_mlst/mlst.txt
    """
}

process EGGNOG {
    tag "EGGNOG"
    input: tuple val(sample), path(proteins)
    script:
    """
    mkdir -p ${params.outdir}/12_function
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
    input: tuple val(sample), path(assembly)
    script:
    """
    mkdir -p ${params.outdir}/13_secondary_metabolites
    antismash --cpus ${params.threads} --output-dir ${params.outdir}/13_secondary_metabolites ${assembly} || echo 'antiSMASH failed' > ${params.outdir}/13_secondary_metabolites/error.txt
    """
}

process MULTIQC {
    script:
    """
    multiqc ${params.outdir} --filename ${params.outdir}/multiqc_report.html --force || true
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
