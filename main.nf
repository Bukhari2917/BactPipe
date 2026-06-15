#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "data/raw_fastq/*_{R1,R2}*.fastq.gz"
params.outdir = "results"
params.threads = 8
params.memory = "32.GB"

workflow {
    Channel.fromFilePairs(params.reads, size: 2).set { read_pairs }
    FASTQC(read_pairs)
    TRIMMING(read_pairs)
    ASSEMBLY(TRIMMING.out.trimmed_reads)
    QUAST(ASSEMBLY.out.assembly)
    PROKKA(ASSEMBLY.out.assembly)
    ABRICATE_AMR(ASSEMBLY.out.assembly)
    ABRICATE_VIR(ASSEMBLY.out.assembly)
    MLST(ASSEMBLY.out.assembly)
    MULTIQC(FASTQC.out.reports.collect())
}

process FASTQC {
    tag "FastQC: ${sample_id}"
    publishDir "${params.outdir}/01_fastqc", mode: 'copy'
    input: tuple val(sample_id), path(reads)
    output: path "*.html" into fastqc_reports
    script: "fastqc ${reads} -t ${params.threads}"
}

process TRIMMING {
    tag "Trimming: ${sample_id}"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    input: tuple val(sample_id), path(reads)
    output: tuple val(sample_id), path("*_trimmed.fastq.gz") into trimmed_reads
    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    fastp -i ${r1} -I ${r2} \\
          -o ${sample_id}_R1_trimmed.fastq.gz \\
          -O ${sample_id}_R2_trimmed.fastq.gz \\
          --detect_adapter_for_pe \\
          --cut_front --cut_tail \\
          --cut_window_size 4 --cut_mean_quality 20 \\
          --length_required 50 \\
          --thread ${params.threads}
    """
}

process ASSEMBLY {
    tag "Assembly: ${sample_id}"
    publishDir "${params.outdir}/03_assembly", mode: 'copy'
    input: tuple val(sample_id), path(trimmed_reads)
    output: path "contigs.fasta" into assembly_results
    script:
    def r1 = trimmed_reads[0]
    def r2 = trimmed_reads[1]
    """
    spades.py -1 ${r1} -2 ${r2} \\
              -o spades_output \\
              --isolate \\
              -t ${params.threads} \\
              -m ${params.memory.replace('GB','')} \\
              --rename
    
    awk '/^>/ {print ">contig_" ++i; next} {print}' spades_output/contigs.fasta > contigs.fasta
    """
}

process QUAST {
    tag "QUAST: ${sample_id}"
    publishDir "${params.outdir}/04_quast", mode: 'copy'
    input: path assembly
    output: path "quast_results"
    script: "quast.py ${assembly} -o quast_results -t ${params.threads}"
}

process PROKKA {
    tag "Annotation: ${sample_id}"
    publishDir "${params.outdir}/05_annotation", mode: 'copy'
    input: path assembly
    output: path "*.gff", path "*.gbk", path "*.faa"
    script:
    """
    prokka ${assembly} \\
           --outdir prokka_out \\
           --prefix sample \\
           --kingdom Bacteria \\
           --centre X \\
           --cpus ${params.threads} \\
           --force
    
    cp prokka_out/* ./
    """
}

process ABRICATE_AMR {
    tag "AMR: ${sample_id}"
    publishDir "${params.outdir}/06_amr", mode: 'copy'
    input: path assembly
    output: path "amr_card.tsv"
    script: "abricate ${assembly} --db card > amr_card.tsv"
}

process ABRICATE_VIR {
    tag "Virulence: ${sample_id}"
    publishDir "${params.outdir}/07_virulence", mode: 'copy'
    input: path assembly
    output: path "virulence.tsv"
    script: "abricate ${assembly} --db vfdb > virulence.tsv"
}

process MLST {
    tag "MLST: ${sample_id}"
    publishDir "${params.outdir}/08_mlst", mode: 'copy'
    input: path assembly
    output: path "mlst.txt"
    script: "mlst ${assembly} > mlst.txt"
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    input: path fastqc_reports
    output: path "multiqc_report.html"
    script: "multiqc . -o multiqc_report"
}
