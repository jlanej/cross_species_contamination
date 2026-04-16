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
    path("taxa_matrix_cpm.tsv"),          emit: matrix_cpm
    path("taxa_matrix_raw.tsv"),          emit: matrix_raw
    path("taxa_matrix_*_?.tsv"),          emit: rank_matrices, optional: true
    path("aggregation_metadata.json"),    emit: metadata
    path("rank_filter_metadata.json"),    emit: rank_metadata

    script:
    def min_reads_arg = params.min_reads ? "--min-reads ${params.min_reads}" : ''
    def rank_arg      = params.rank_filter ? "--rank-filter ${params.rank_filter}" : ''
    """
    csc-aggregate \\
        ${reports} \\
        -o . \\
        ${min_reads_arg} \\
        ${rank_arg}
    """
}
