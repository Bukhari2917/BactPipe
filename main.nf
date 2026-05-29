#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = null
params.outdir = "./results"
params.sample = "BACTERIA"
params.threads = 8
params.help = false

if (params.help) {
    log.info "BactPipe - Bacterial Genome Analysis"
    exit 0
}

include { FASTQC } from './modules/fastqc.nf'
include { TRIM } from './modules/trim.nf'
include { ASSEMBLE } from './modules/assemble.nf'
include { QUAST } from './modules/quast.nf'
include { BUSCO } from './modules/busco.nf'
include { PRODIGAL } from './modules/prodigal.nf'
include { BARRNAP } from './modules/barrnap.nf'
include { TRNA } from './modules/trna.nf'
include { EGGNOG } from './modules/eggnog.nf'
include { ABRICATE_AMR } from './modules/abricate_amr.nf'
include { ABRICATE_VIR } from './modules/abricate_vir.nf'
include { CRISPR } from './modules/crispr.nf'
include { MLST } from './modules/mlst.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow {
    // Create channel from input reads
    Channel.fromFilePairs(params.reads)
        .map { id, reads -> 
            def sample_name = params.sample == "BACTERIA" ? id : params.sample
            [sample_name, reads[0], reads[1]]
        }
        .set { reads_ch }
    
    // Run QC
    FASTQC(reads_ch)
    TRIM(reads_ch)
    
    // Run assembly - TRIM.out is already a tuple [sample, r1, r2]
    ASSEMBLE(TRIM.out)
    
    // Run analysis on assembly output
    QUAST(ASSEMBLE.out.map { [it[0], it[1]] })
    BUSCO(ASSEMBLE.out.map { [it[0], it[1]] })
    PRODIGAL(ASSEMBLE.out.map { [it[0], it[1]] })
    BARRNAP(ASSEMBLE.out.map { [it[0], it[1]] })
    TRNA(ASSEMBLE.out.map { [it[0], it[1]] })
    ABRICATE_AMR(ASSEMBLE.out)
    ABRICATE_VIR(ASSEMBLE.out)
    CRISPR(ASSEMBLE.out)
    MLST(ASSEMBLE.out)
    
    // eggNOG needs just the protein file
    EGGNOG(PRODIGAL.out.map { it[1] })
    
    // MultiQC report
    MULTIQC(
        FASTQC.out.collect(), 
        TRIM.out.map { it[2] }.collect(), 
        ASSEMBLE.out.map { it[1] }.collect()
    )
    
    log.info "✅ Pipeline completed! Results: ${params.outdir}"
}
