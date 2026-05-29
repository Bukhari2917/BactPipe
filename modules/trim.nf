process TRIM {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq")
    path "fastp.html"
    
    script:
    """
    fastp -i ${r1} -I ${r2} -o ${sample_id}_R1.fastq -O ${sample_id}_R2.fastq \\
          --html fastp.html -q 20 -l 50 --thread ${task.cpus}
    """
}