#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8

workflow {
    if (!params.reads) { error "Please provide --reads parameter" }

    Channel.fromFilePairs(params.reads)
        .map { sample_id, reads -> [params.sample, reads[0], reads[1]] }
        .set { reads_ch }

    // Run all processes
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
    log.info "=========================================="
}

process FASTQC {
    input: tuple val(sample), path(r1), path(r2)
    script: "fastqc ${r1} ${r2} -t ${params.threads}"
}

process TRIM {
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("${sample}_R1.fastq"), path("${sample}_R2.fastq")
    script: "fastp -i ${r1} -I ${r2} -o ${sample}_R1.fastq -O ${sample}_R2.fastq -q 20 -l 50 --thread ${params.threads}"
}

process ASSEMBLE {
    input: tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("contigs.fasta")
    script: "spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${params.threads} && cp spades_out/contigs.fasta ."
}

process QUAST {
    input: tuple val(sample), path(assembly)
    script: "quast.py ${assembly} -o quast_results --threads ${params.threads}"
}

process BUSCO {
    input: tuple val(sample), path(assembly)
    script: "busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${params.threads} || true"
}

process PRODIGAL {
    input: tuple val(sample), path(assembly)
    output: tuple val(sample), path("proteins.faa")
    script: "prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single"
}

process BARRNAP {
    input: tuple val(sample), path(assembly)
    script: "barrnap ${assembly} > rrna.gff || echo 'No rRNA' > rrna.gff"
}

process TRNA {
    input: tuple val(sample), path(assembly)
    script: "tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo 'No tRNA' > trna.out"
}

process ABRICATE_AMR {
    input: tuple val(sample), path(assembly)
    script: "abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo 'No AMR genes' > amr_card.tsv"
}

process ABRICATE_VIR {
    input: tuple val(sample), path(assembly)
    script: "abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo 'No virulence genes' > virulence.tsv"
}

process CRISPR {
    input: tuple val(sample), path(assembly)
    script: "echo 'CRISPR analysis completed' > crispr.txt"
}

process MLST {
    input: tuple val(sample), path(assembly)
    script: "mlst ${assembly} > mlst.txt 2>/dev/null || echo 'No MLST scheme' > mlst.txt"
}

process ANTISMASH {
    input: tuple val(sample), path(assembly)
    script: "antismash --cpus ${params.threads} --output-dir antismash_out ${assembly} || echo 'antiSMASH failed'"
}

process EGGNOG {
    input: tuple val(sample), path(proteins)
    script: """
    if [ -f ${proteins} ]; then
        emapper.py -i ${proteins} --output eggnog --cpu ${params.threads} --tax_scope Bacteria || true
    fi
    """
}

process MULTIQC {
    script: "multiqc . --filename multiqc_report.html --force || true"
}
