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
    
    grep -v "^#" eggnog.emapper.annotations | cut -f12 | sort | uniq -c | sort -rn > kegg_pathways.txt
    grep -v "^#" eggnog.emapper.annotations | cut -f7 | sort | uniq -c | sort -rn > cog_categories.txt
    grep -v "^#" eggnog.emapper.annotations | cut -f9 | tr ',' '\n' | sort | uniq -c | sort -rn > go_terms.txt
    """
}