#!/usr/bin/env nextflow

// Define input
params.reads = "data/*_R{1,2}.fastq"
params.outdir = "results"

// Set up channels
Channel
    .fromFilePairs(params.reads)
    .set { read_pairs }

process FASTQC {
    tag "FASTQC: ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.html"
    path "*.zip"
    
    script:
    """
    fastqc ${reads} -o .
    """
}

process TRIMMING {
    tag "Trimming: ${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_R1_trimmed.fastq"), path("*_R2_trimmed.fastq")
    
    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \
        ${sample_id}_R1_trimmed.fastq ${sample_id}_R1_unpaired.fastq \
        ${sample_id}_R2_trimmed.fastq ${sample_id}_R2_unpaired.fastq \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    """
}

process ASSEMBLY {
    tag "Assembly: ${sample_id}"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)
    
    output:
    path "${sample_id}/"
    path "contigs.fasta" into assembly_contigs
    
    script:
    """
    spades.py -1 ${reads_R1} -2 ${reads_R2} -o ${sample_id} --isolate
    cp ${sample_id}/contigs.fasta contigs.fasta
    """
}

process QUAST {
    tag "QUAST: ${sample_id}"
    publishDir "${params.outdir}/quast", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "quast_results/"
    
    script:
    """
    quast.py ${contigs} -o quast_results
    """
}

process PRODIGAL {
    tag "Prodigal: ${sample_id}"
    publishDir "${params.outdir}/prodigal", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "genes.gff"
    path "proteins.faa"
    path "nucleotides.fna"
    
    script:
    """
    prodigal -i ${contigs} -a proteins.faa -d nucleotides.fna -f gff -o genes.gff
    """
}

process TRNA {
    tag "tRNA: ${sample_id}"
    publishDir "${params.outdir}/trna", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "trna_results.txt"
    
    script:
    """
    tRNAscan-SE ${contigs} -o trna_results.txt
    """
}

process RRNA {
    tag "rRNA: ${sample_id}"
    publishDir "${params.outdir}/rrna", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "rrna_results.txt"
    
    script:
    """
    barrnap ${contigs} > rrna_results.txt
    """
}

process ABRICATE_AMR {
    tag "AMR: ${sample_id}"
    publishDir "${params.outdir}/amr", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "amr_results.txt"
    
    script:
    """
    abricate ${contigs} --db resfinder > amr_results.txt
    """
}

process ABRICATE_VIRULENCE {
    tag "Virulence: ${sample_id}"
    publishDir "${params.outdir}/virulence", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "virulence_results.txt"
    
    script:
    """
    abricate ${contigs} --db vfdb > virulence_results.txt
    """
}

process CRISPR {
    tag "CRISPR: ${sample_id}"
    publishDir "${params.outdir}/crispr", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "crispr_results.txt"
    
    script:
    """
    minced ${contigs} > crispr_results.txt 2>&1
    """
}

process MLST {
    tag "MLST: ${sample_id}"
    publishDir "${params.outdir}/mlst", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "mlst_results.txt"
    
    script:
    """
    mlst ${contigs} > mlst_results.txt
    """
}

process ANTISMASH {
    tag "antiSMASH: ${sample_id}"
    publishDir "${params.outdir}/antismash", mode: 'copy'
    
    input:
    path contigs
    
    output:
    path "${sample_id}_antismash/"
    
    script:
    """
    antismash ${contigs} --output-dir ${sample_id}_antismash --genefinding-tool prodigal
    """
}

process EGGNOG {
    tag "eggNOG: ${sample_id}"
    publishDir "${params.outdir}/eggnog", mode: 'copy'
    
    input:
    path proteins_faa from PRODIGAL.out.proteins_faa
    
    output:
    path "eggnog_results/"
    
    script:
    """
    emapper.py -i ${proteins_faa} --output eggnog_results --cpu 4 --dmnd_db eggnog_proteins --data_dir eggnog_data
    """
}

// Workflow connecting everything
workflow {
    // Read QC
    FASTQC(read_pairs)
    
    // Trimming
    TRIMMING(read_pairs)
    TRIMMING.out.map { tuple(sample, r1, r2) }.set { trimmed_reads }
    
    // Assembly
    ASSEMBLY(trimmed_reads)
    
    // Assembly evaluation
    QUAST(ASSEMBLY.out.contigs)
    
    // Gene prediction
    PRODIGAL(ASSEMBLY.out.contigs)
    
    // RNA genes
    TRNA(ASSEMBLY.out.contigs)
    RRNA(ASSEMBLY.out.contigs)
    
    // Functional analyses
    ABRICATE_AMR(ASSEMBLY.out.contigs)
    ABRICATE_VIRULENCE(ASSEMBLY.out.contigs)
    CRISPR(ASSEMBLY.out.contigs)
    MLST(ASSEMBLY.out.contigs)
    ANTISMASH(ASSEMBLY.out.contigs)
    
    // Functional annotation (needs Prodigal proteins)
    EGGNOG(PRODIGAL.out.proteins_faa)
}
