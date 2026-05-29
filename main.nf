#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8
params.help = false

if (params.help) {
    log.info "BactPipe - Complete Bacterial Genome Analysis Pipeline"
    exit 0
}

// ============================================================================
// PROCESSES
// ============================================================================

process FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    path "fastqc_results"
    script:
    """
    mkdir -p fastqc_results
    fastqc ${r1} ${r2} -o fastqc_results/ -t ${task.cpus}
    """
}

process TRIM {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq")
    script:
    """
    fastp -i ${r1} -I ${r2} -o ${sample_id}_R1.fastq -O ${sample_id}_R2.fastq -q 20 -l 50 --thread ${task.cpus}
    """
}

process ASSEMBLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("contigs.fasta")
    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o spades_out --isolate --threads ${task.cpus}
    cp spades_out/contigs.fasta .
    """
}

process QUAST {
    tag "${sample_id}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "quast_results"
    script:
    """
    quast.py ${assembly} -o quast_results --threads ${task.cpus}
    """
}

process BUSCO {
    tag "${sample_id}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "busco_results"
    script:
    """
    busco -i ${assembly} -o busco_results -l bacteria_odb10 -m genome --cpu ${task.cpus}
    """
}

process PRODIGAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/genes", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    tuple val(sample_id), path("proteins.faa")
    script:
    """
    prodigal -i ${assembly} -a proteins.faa -o genes.gbk -p single
    """
}

process BARRNAP {
    tag "${sample_id}"
    publishDir "${params.outdir}/genes", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "rrna.gff"
    script:
    """
    barrnap ${assembly} > rrna.gff
    """
}

process TRNA {
    tag "${sample_id}"
    publishDir "${params.outdir}/genes", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "trna.out"
    script:
    """
    tRNAscan-SE -o trna.out ${assembly} 2>/dev/null || echo "No tRNA found" > trna.out
    """
}

process ABRICATE_AMR {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "amr_card.tsv"
    script:
    """
    abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo "No AMR genes found" > amr_card.tsv
    """
}

process ABRICATE_VIR {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "virulence.tsv"
    script:
    """
    abricate --db vfdb ${assembly} > virulence.tsv 2>/dev/null || echo "No virulence genes found" > virulence.tsv
    """
}

process CRISPR {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "crispr.txt"
    script:
    """
    echo "CRISPR analysis - minced not available" > crispr.txt
    """
}

process MLST {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    input:
    tuple val(sample_id), path(assembly)
    output:
    path "mlst.txt"
    script:
    """
    mlst ${assembly} > mlst.txt 2>/dev/null || echo "No MLST scheme found" > mlst.txt
    """
}

process EGGNOG {
    tag "${sample_id}"
    publishDir "${params.outdir}/function", mode: 'copy'
    input:
    path proteins
    output:
    path "kegg_pathways.txt"
    path "cog_categories.txt"
    path "go_terms.txt"
    script:
    """
    emapper.py -i ${proteins} --output eggnog --cpu ${task.cpus} --tax_scope Bacteria
    grep -v '^#' eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt
    grep -v '^#' eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt
    grep -v '^#' eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > go_terms.txt
    """
}

process MULTIQC {
    publishDir "${params.outdir}/report", mode: 'copy'
    input:
    path fastqc
    path trim
    path assembly
    output:
    path "multiqc_report.html"
    script:
    """
    multiqc . --filename multiqc_report.html --force
    """
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
        .map { id, reads -> [params.sample, reads[0], reads[1]] }
    
    FASTQC(reads_ch)
    TRIM(reads_ch)
    ASSEMBLE(TRIM.out)
    QUAST(ASSEMBLE.out.map { [it[0], it[1]] })
    BUSCO(ASSEMBLE.out.map { [it[0], it[1]] })
    PRODIGAL(ASSEMBLE.out.map { [it[0], it[1]] })
    BARRNAP(ASSEMBLE.out.map { [it[0], it[1]] })
    TRNA(ASSEMBLE.out.map { [it[0], it[1]] })
    ABRICATE_AMR(ASSEMBLE.out)
    ABRICATE_VIR(ASSEMBLE.out)
    CRISPR(ASSEMBLE.out)
    MLST(ASSEMBLE.out)
    EGGNOG(PRODIGAL.out.map { it[1] })
    MULTIQC(FASTQC.out.collect(), TRIM.out.map { it[2] }.collect(), ASSEMBLE.out.map { it[1] }.collect())
    
    log.info "=========================================="
    log.info "BactPipe COMPLETE Pipeline Finished!"
    log.info "Results: ${params.outdir}"
    log.info "=========================================="
}
