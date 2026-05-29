process CRISPR {
    tag "${sample_id}"
    publishDir "${params.outdir}/specialized", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "crispr.txt"
    
    script:
    """
    minced -i ${assembly} -o minced_out 2>/dev/null || echo "No CRISPR found" > crispr.txt
    if [ -f minced_out.gff ]; then
        echo "CRISPR arrays found" > crispr.txt
    else
        echo "No CRISPR arrays found" > crispr.txt
    fi
    """
}