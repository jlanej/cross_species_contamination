/*
 * detect.nf – Process definition for outlier detection
 *
 * Wraps the ``csc-detect`` CLI to perform statistical outlier
 * analysis on the aggregated sample-by-taxon matrix.
 *
 * AI assistance acknowledgment: developed with AI assistance.
 */

process DETECT_OUTLIERS {
    tag "detect"
    publishDir "${params.outdir}/detect", mode: 'copy'

    cpus   params.detect_cpus   ?: 2
    memory params.detect_memory ?: '4 GB'
    time   params.detect_time   ?: '1h'

    input:
    path(matrix)

    output:
    path("flagged_samples.tsv"),  emit: flagged
    path("qc_summary.json"),      emit: qc_summary
    path("quarantine_list.txt"),   emit: quarantine

    script:
    def method_arg    = params.detect_method ? "--method ${params.detect_method}" : ''
    def mad_arg       = params.mad_threshold ? "--mad-threshold ${params.mad_threshold}" : ''
    def iqr_arg       = params.iqr_multiplier ? "--iqr-multiplier ${params.iqr_multiplier}" : ''
    def kitome_arg    = params.kitome_taxa ? "--kitome-taxa ${params.kitome_taxa}" : ''
    def bg_arg        = params.no_subtract_background ? '--no-subtract-background' : ''
    def rank_arg      = params.rank_filter ? "--rank-filter ${params.rank_filter}" : ''
    """
    csc-detect \\
        ${matrix} \\
        -o . \\
        ${method_arg} \\
        ${mad_arg} \\
        ${iqr_arg} \\
        ${kitome_arg} \\
        ${bg_arg} \\
        ${rank_arg}
    """
}
