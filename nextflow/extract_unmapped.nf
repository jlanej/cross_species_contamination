#!/usr/bin/env nextflow

/*
 * extract_unmapped.nf – Batch extraction of unmapped reads from BAM/CRAM files
 *
 * Usage:
 *   nextflow run nextflow/extract_unmapped.nf \
 *       --input_csv samples.csv \
 *       --outdir results/
 *
 * The CSV must have columns: sample_id, file  (and optionally: reference)
 *
 * AI assistance acknowledgment: developed with AI assistance. Best practices
 * in bioinformatics should always supersede implementation specifics here.
 */

nextflow.enable.dsl = 2

params.input_csv  = null
params.outdir     = 'results'
params.mapq       = null       // set to an integer to include low-MAPQ reads
params.threads    = 4
params.reference  = null       // global reference; per-sample overrides in CSV

// ── Processes ────────────────────────────────────────────────────────────────

process EXTRACT_UNMAPPED {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sample_id), path(input_file), val(reference)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: fastqs

    script:
    def ref_arg  = reference ? "--reference ${reference}" : ''
    def mapq_arg = params.mapq ? "--mapq ${params.mapq}" : ''
    """
    extract-unmapped \\
        ${input_file} \\
        -o . \\
        --sample-id ${sample_id} \\
        --threads ${task.cpus} \\
        ${ref_arg} \\
        ${mapq_arg}
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {
    if (!params.input_csv) {
        error "Please supply --input_csv with columns: sample_id, file [, reference]"
    }

    Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            def ref = row.reference ?: params.reference ?: ''
            tuple(row.sample_id, file(row.file), ref)
        }
        .set { samples_ch }

    EXTRACT_UNMAPPED(samples_ch)
}
