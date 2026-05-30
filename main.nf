cd /data/sayed/bacteria_test/BactPipe

cat > main.nf << 'EOF'
#!/usr/bin/env nextflow

params.reads = "data/*_R{1,2}.fastq"
params.outdir = "/data/sayed/bacteria_test/Out_Results"

Channel.fromFilePairs(params.reads).set { read_pairs }

process FASTQC {
    publishDir "${params.outdir}/01_fastqc", mode: 'copy', overwrite: true
    input: tuple val(sample_id), path(reads)
    output: path "*.html"
    script: "fastqc ${reads} -o ."
}

process ASSEMBLY {
    publishDir "${params.outdir}/03_assembly", mode: 'copy', overwrite: true
    input: tuple val(sample_id), path(reads)
    output: path "contigs.fasta", emit: contigs
    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out --isolate -t ${task.cpus}
    cp spades_out/contigs.fasta .
    """
}

process PRODIGAL {
    publishDir "${params.outdir}/05_prodigal", mode: 'copy', overwrite: true
    input: path contigs
    output: path "proteins.faa"
    script: "prodigal -i ${contigs} -a proteins.faa -o genes.gff -p meta"
}

process ABRICATE_AMR {
    publishDir "${params.outdir}/08_amr", mode: 'copy', overwrite: true
    input: path contigs
    output: path "amr_results.txt"
    script: "abricate ${contigs} --db resfinder > amr_results.txt"
}

process ABRICATE_VIRULENCE {
    publishDir "${params.outdir}/09_virulence", mode: 'copy', overwrite: true
    input: path contigs
    output: path "virulence_results.txt"
    script: "abricate ${contigs} --db vfdb > virulence_results.txt"
}

process CRISPR {
    publishDir "${params.outdir}/10_crispr", mode: 'copy', overwrite: true
    input: path contigs
    output: path "crispr_results.txt"
    script: "minced ${contigs} > crispr_results.txt 2>&1"
}

process MLST {
    publishDir "${params.outdir}/11_mlst", mode: 'copy', overwrite: true
    input: path contigs
    output: path "mlst_results.txt"
    script: "mlst ${contigs} > mlst_results.txt"
}

workflow {
    FASTQC(read_pairs)
    ASSEMBLY(read_pairs)
    PRODIGAL(ASSEMBLY.out.contigs)
    ABRICATE_AMR(ASSEMBLY.out.contigs)
    ABRICATE_VIRULENCE(ASSEMBLY.out.contigs)
    CRISPR(ASSEMBLY.out.contigs)
    MLST(ASSEMBLY.out.contigs)
}
EOF

# Clean everything
rm -rf work .nextflow*
rm -rf /data/sayed/bacteria_test/Out_Results/*

# Run pipeline fresh
nextflow run main.nf --reads 'data/*_R{1,2}.fastq'
