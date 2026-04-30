#!/usr/bin/env nextflow

/*
 * main.nf – End-to-end cross-species contamination detection pipeline
 *
 * Orchestrates four stages:
 *   1. EXTRACT   – pull unmapped/low-MAPQ reads from BAM/CRAM
 *   2. CLASSIFY  – Kraken2 taxonomic classification of extracted reads
 *   3. AGGREGATE – build a sample-by-taxon matrix from classification reports
 *   4. DETECT    – statistical outlier detection for contamination flagging
 *
 * A final SUMMARY step produces an HTML report and a MultiQC-compatible
 * YAML file that can be consumed by ``multiqc .``.
 *
 * Usage:
 *   nextflow run nextflow/main.nf \
 *       --input_csv  samples.csv \
 *       --kraken2_db /data/kraken2/PlusPF \
 *       --outdir     results/
 *
 * The input CSV must contain columns: sample_id, file
 * Optional columns: reference, fastq2
 *
 * AI assistance acknowledgment: developed with AI assistance. Best practices
 * in bioinformatics should always supersede implementation specifics here.
 */

nextflow.enable.dsl = 2

// ── Include modular processes ────────────────────────────────────────────────
include { EXTRACT_UNMAPPED  } from './modules/extract'
include { CLASSIFY_READS    } from './modules/classify'
include { AGGREGATE_REPORTS } from './modules/aggregate'
include { DETECT_OUTLIERS   } from './modules/detect'
include { PIPELINE_SUMMARY  } from './modules/summary'

// ── Pipeline parameters (defaults; override via CLI or nextflow.config) ──────
params.input_csv  = null
params.outdir     = 'results'

// Extract
params.mapq       = null       // null = unmapped-only; integer = include low-MAPQ
params.extract_cpus   = 4
params.extract_memory = '4 GB'
params.extract_time   = '2h'

// Classify
params.kraken2_db     = null
params.confidence     = 0.0
params.memory_mapping = false
params.classify_cpus   = 4
params.classify_memory = '8 GB'
params.classify_time   = '4h'

// Aggregate
params.min_reads       = 0
params.rank_filter     = 'S G F'
params.db_path         = null       // path to Kraken2 DB for domain annotation (optional)
params.detect_matrix   = 'cpm'    // 'cpm' or 'raw'
// High-confidence tier (per-read confidence recomputation; requires --db-path).
// Multiple values can be supplied space-separated (e.g. '0.1 0.5'); each
// produces a parallel matrix set with suffix _conf{T}.  An empty string
// disables high-confidence aggregation (only the sensitive tier is built).
params.confidence_thresholds = ''
params.aggregate_cpus   = 2
params.aggregate_memory = '4 GB'
params.aggregate_time   = '1h'

// Detect
params.detect_method          = 'all'
params.mad_threshold          = 3.5
params.iqr_multiplier         = 1.5
params.gmm_threshold          = 0.5
params.kitome_taxa            = null
params.no_subtract_background = false
params.no_abs_detection       = false
params.detect_cpus   = 2
params.detect_memory = '4 GB'
params.detect_time   = '1h'

// ── Workflow ─────────────────────────────────────────────────────────────────
workflow {

    // --- Parameter validation ------------------------------------------------
    if (!params.input_csv) {
        error "ERROR: --input_csv is required. Supply a CSV with columns: sample_id, file"
    }
    if (!params.kraken2_db) {
        error "ERROR: --kraken2_db is required. Supply the path to a Kraken2 database directory."
    }
    if (!(params.detect_matrix in ['cpm', 'raw'])) {
        error "ERROR: --detect_matrix must be either 'cpm' or 'raw'."
    }

    // --- 1. Read sample sheet ------------------------------------------------
    Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def ref = row.reference ?: ''
            tuple(row.sample_id, file(row.file), ref)
        }
        .set { samples_ch }

    // --- 2. Extract unmapped reads -------------------------------------------
    EXTRACT_UNMAPPED(samples_ch)

    // --- 3. Classify with Kraken2 --------------------------------------------
    CLASSIFY_READS(EXTRACT_UNMAPPED.out.fastqs)

    // --- 4. Aggregate all reports into a single matrix -----------------------
    // Collect all report files into a single list for the aggregate process
    CLASSIFY_READS.out.reports
        .map { sample_id, report -> report }
        .collect()
        .set { all_reports_ch }

    // Collect per-read kraken2 outputs alongside reports.  These are
    // required when the high-confidence tier is enabled
    // (params.confidence_thresholds non-empty); they are otherwise
    // ignored by csc-aggregate.
    CLASSIFY_READS.out.outputs
        .map { sample_id, output -> output }
        .collect()
        .ifEmpty([])
        .set { all_outputs_ch }

    AGGREGATE_REPORTS(all_reports_ch, all_outputs_ch)

    // --- 5. Detect outliers --------------------------------------------------
    def detect_input_ch = params.detect_matrix == 'raw'
        ? AGGREGATE_REPORTS.out.matrix_raw
        : AGGREGATE_REPORTS.out.matrix_cpm

    // Stage per-rank matrices alongside the primary matrix so that
    // csc-detect --rank-filter can locate sibling files in the work dir.
    def rank_matrices_ch = AGGREGATE_REPORTS.out.rank_matrices
        .flatten()
        .collect()
        .ifEmpty([])

    // Stage confidence-tier matrices (and their per-rank siblings) so
    // that csc-detect can auto-discover them and run detection per tier.
    def tier_matrices_ch = AGGREGATE_REPORTS.out.tier_matrices
        .flatten()
        .collect()
        .ifEmpty([])

    // Stage the absolute-burden matrix and its rank/tier siblings so
    // that csc-detect can run its absolute-burden side pass when an
    // abs sibling is present.  The abs matrix is optional (only
    // emitted when csc-aggregate received --idxstats sidecars); when
    // absent, csc-detect simply skips the side pass.
    def abs_matrix_ch = AGGREGATE_REPORTS.out.matrix_abs
        .ifEmpty([])
    // The rank- and tier-channels above already include any abs
    // siblings (their globs match taxa_matrix_*_?.tsv and
    // taxa_matrix_*_conf*.tsv).  We pass empty placeholder channels
    // here so the DETECT_OUTLIERS process declaration matches.
    def abs_rank_matrices_ch = Channel.value([])
    def abs_tier_matrices_ch = Channel.value([])

    DETECT_OUTLIERS(
        detect_input_ch,
        rank_matrices_ch,
        tier_matrices_ch,
        abs_matrix_ch,
        abs_rank_matrices_ch,
        abs_tier_matrices_ch,
    )

    // --- 6. Pipeline summary report ------------------------------------------
    PIPELINE_SUMMARY(
        DETECT_OUTLIERS.out.flagged,
        DETECT_OUTLIERS.out.qc_summary,
        DETECT_OUTLIERS.out.quarantine,
        AGGREGATE_REPORTS.out.metadata,
    )
}

// ── Completion handler ───────────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? 'SUCCEEDED' : 'FAILED'
    log.info """
    ============================================================
      CSC Pipeline ${status}
    ============================================================
      Completed at : ${workflow.complete}
      Duration     : ${workflow.duration}
      Exit status  : ${workflow.exitStatus}
      Work dir     : ${workflow.workDir}
      Output dir   : ${params.outdir}
    ============================================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    ============================================================
      CSC Pipeline ERROR
    ============================================================
      Error message: ${workflow.errorMessage}
      Work dir     : ${workflow.workDir}
    ============================================================
    """.stripIndent()
}
