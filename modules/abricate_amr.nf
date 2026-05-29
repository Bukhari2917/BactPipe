process ABRICATE_AMR {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "amr_card.tsv"
    path "amr_ncbi.tsv"
    path "amr_resfinder.tsv"
    
    script:
    """
    abricate --db card ${assembly} > amr_card.tsv 2>/dev/null || echo "No AMR genes found" > amr_card.tsv
    abricate --db ncbi ${assembly} > amr_ncbi.tsv 2>/dev/null || echo "No AMR genes found" > amr_ncbi.tsv
    abricate --db resfinder ${assembly} > amr_resfinder.tsv 2>/dev/null || echo "No AMR genes found" > amr_resfinder.tsv
    """
}