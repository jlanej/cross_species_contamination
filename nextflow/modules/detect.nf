/*
 * detect.nf – Process definition for outlier detection
 *
 * Wraps the ``csc-detect`` CLI to perform statistical outlier
 * analysis on the aggregated sample-by-taxon matrix.
 *
 * rank_matrices contains the per-rank typed matrices produced by
 * AGGREGATE_REPORTS (e.g. taxa_matrix_cpm_S.tsv).  Staging them
 * in the same work directory lets csc-detect find sibling rank
 * files automatically when --rank-filter is active.
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
    path(rank_matrices)
    path(tier_matrices)

    output:
    path("flagged_samples.tsv"),                  emit: flagged
    path("qc_summary.json"),                      emit: qc_summary
    path("quarantine_list.txt"),                  emit: quarantine
    // csc-detect writes per-rank results into subdirectories named after the
    // rank code (e.g. S/, G/, F/) when --rank-filter is active.
    path("*/flagged_samples.tsv"),                emit: rank_flagged,    optional: true
    path("*/qc_summary.json"),                    emit: rank_qc,         optional: true
    path("*/quarantine_list.txt"),                emit: rank_quarantine, optional: true
    // High-confidence tier outputs (e.g. conf0p50/) and their rank
    // sub-subdirectories (e.g. conf0p50/S/).
    path("conf*/flagged_samples.tsv"),            emit: tier_flagged,    optional: true
    path("conf*/qc_summary.json"),                emit: tier_qc,         optional: true
    path("conf*/quarantine_list.txt"),            emit: tier_quarantine, optional: true
    path("conf*/*/flagged_samples.tsv"),          emit: tier_rank_flagged, optional: true

    script:
    def method_arg    = params.detect_method ? "--method ${params.detect_method}" : ''
    def mad_arg       = params.mad_threshold ? "--mad-threshold ${params.mad_threshold}" : ''
    def iqr_arg       = params.iqr_multiplier ? "--iqr-multiplier ${params.iqr_multiplier}" : ''
    def gmm_arg       = params.gmm_threshold ? "--gmm-threshold ${params.gmm_threshold}" : ''
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
        ${gmm_arg} \\
        ${kitome_arg} \\
        ${bg_arg} \\
        ${rank_arg}
    """
}
