/*
 * aggregate.nf – Process definition for report aggregation
 *
 * Wraps the ``csc-aggregate`` CLI to build a sample-by-taxon matrix
 * from the collection of Kraken2 reports produced by CLASSIFY_READS.
 *
 * AI assistance acknowledgment: developed with AI assistance.
 */

process AGGREGATE_REPORTS {
    tag "aggregate"
    publishDir "${params.outdir}/aggregate", mode: 'copy'

    cpus   params.aggregate_cpus   ?: 2
    memory params.aggregate_memory ?: '4 GB'
    time   params.aggregate_time   ?: '1h'

    input:
    path(reports)

    output:
    path("taxa_matrix.tsv"),           emit: matrix
    path("aggregation_metadata.json"), emit: metadata

    script:
    def norm_arg      = params.no_normalize ? '--no-normalize' : ''
    def min_reads_arg = params.min_reads ? "--min-reads ${params.min_reads}" : ''
    """
    csc-aggregate \\
        ${reports} \\
        -o . \\
        ${min_reads_arg} \\
        ${norm_arg}
    """
}
