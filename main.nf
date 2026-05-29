workflow {
    // Read input files
    Channel.fromFilePairs(params.reads)
        .map { id, reads -> [params.sample, reads[0], reads[1]] }
        .set { reads_ch }
    
    // Quality control
    FASTQC(reads_ch)
    
    // Trimming
    TRIM(reads_ch)
    
    // Assembly - TRIM.out is a channel, we need to use it directly
    TRIM.out
        .map { sample, r1, r2 -> [sample, r1, r2] }
        .set { assembly_input }
    
    ASSEMBLE(assembly_input)
    
    // Remaining processes
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
    
    log.info "Pipeline completed! Results: ${params.outdir}"
}
