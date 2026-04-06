#!/usr/bin/env nextflow

/*
 * classify.nf – Kraken2 taxonomic classification of extracted reads
 *
 * Usage:
 *   nextflow run nextflow/classify.nf \
 *       --input_csv samples.csv \
 *       --kraken2_db /data/kraken2/PlusPF \
 *       --outdir results/
 *
 * The CSV must have columns: sample_id, fastq
 * For paired-end data, add a fastq2 column.
 *
 * This pipeline can be chained after extract_unmapped.nf by supplying
 * the extracted FASTQ paths in the input CSV.
 *
 * AI assistance acknowledgment: developed with AI assistance. Best practices
 * in bioinformatics should always supersede implementation specifics here.
 */

nextflow.enable.dsl = 2

params.input_csv    = null
params.outdir       = 'results'
params.kraken2_db   = null
params.confidence   = 0.0
params.threads      = 4
params.memory_mapping = false

// ── Processes ────────────────────────────────────────────────────────────────

process CLASSIFY_READS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sample_id), path(fastq1), val(fastq2)

    output:
    tuple val(sample_id), path("*.kraken2.report.txt"), path("*.kraken2.output.txt"), emit: results

    script:
    def paired_arg  = fastq2 ? "--paired ${fastq2}" : ''
    def mem_arg     = params.memory_mapping ? '--memory-mapping' : ''
    """
    csc-classify \\
        ${fastq1} ${paired_arg} \\
        --db ${params.kraken2_db} \\
        -o . \\
        --sample-id ${sample_id} \\
        --confidence ${params.confidence} \\
        --threads ${task.cpus} \\
        ${mem_arg}
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {
    if (!params.input_csv) {
        error "Please supply --input_csv with columns: sample_id, fastq [, fastq2]"
    }
    if (!params.kraken2_db) {
        error "Please supply --kraken2_db with the path to a Kraken2 database"
    }

    Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            def fq2 = row.fastq2 ?: ''
            tuple(row.sample_id, file(row.fastq), fq2 ? file(fq2) : '')
        }
        .set { samples_ch }

    CLASSIFY_READS(samples_ch)
}
