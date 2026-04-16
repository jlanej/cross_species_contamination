/*
 * summary.nf – Pipeline summary report
 *
 * Collects output from all stages and writes a MultiQC-compatible
 * summary report (YAML header + TSV data) so downstream MultiQC
 * runs can pick up the results.
 *
 * AI assistance acknowledgment: developed with AI assistance.
 */

process PIPELINE_SUMMARY {
    tag "summary"
    publishDir "${params.outdir}", mode: 'copy'

    cpus   1
    memory '1 GB'
    time   '30m'

    input:
    path(flagged)
    path(qc_summary)
    path(quarantine)
    path(agg_metadata)

    output:
    path("pipeline_report.html"),    emit: report
    path("csc_pipeline_mqc.yaml"),   emit: multiqc_data

    script:
    """
    #!/usr/bin/env python3
    import json
    import html
    from pathlib import Path
    from datetime import datetime, timezone

    # ── Load inputs ──────────────────────────────────────────────────
    qc = json.loads(Path("${qc_summary}").read_text())
    agg = json.loads(Path("${agg_metadata}").read_text())
    flagged_lines = Path("${flagged}").read_text().strip().splitlines()
    quarantine_ids = [
        ln.strip() for ln in Path("${quarantine}").read_text().strip().splitlines()
        if ln.strip()
    ]

    # ── Derive stats ─────────────────────────────────────────────────
    n_flagged_rows = max(len(flagged_lines) - 1, 0)  # subtract header
    n_quarantined  = len(quarantine_ids)
    n_samples      = agg.get("sample_count", "N/A")
    n_taxa         = agg.get("taxon_count", "N/A")
    method         = qc.get("method", "N/A")
    total_analysed = qc.get("total_samples", "N/A")
    total_taxa     = qc.get("total_taxa_analysed", "N/A")
    flagged_count  = qc.get("flagged_count", n_flagged_rows)

    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # ── MultiQC custom-content YAML ──────────────────────────────────
    mqc = {
        "id": "csc_pipeline",
        "section_name": "Cross-Species Contamination",
        "description": "Summary of the CSC detection pipeline run.",
        "plot_type": "table",
        "pconfig": {"id": "csc_summary_table", "title": "CSC Pipeline Summary"},
        "data": {
            "Pipeline Run": {
                "Samples": n_samples,
                "Taxa in Matrix": n_taxa,
                "Detection Method": method,
                "Flagged (sample, taxon) Pairs": flagged_count,
                "Quarantined Samples": n_quarantined,
                "Timestamp (UTC)": timestamp,
            }
        },
    }

    Path("csc_pipeline_mqc.yaml").write_text(
        json.dumps({"custom_data": {"csc_pipeline": mqc}}, indent=2)
    )

    # ── HTML report ──────────────────────────────────────────────────
    quarantine_items = "".join(
        f"<li>{html.escape(s)}</li>" for s in quarantine_ids
    ) or "<li>None</li>"

    report = f\"\"\"<!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="utf-8"/>
      <title>CSC Pipeline Report</title>
      <style>
        body {{ font-family: Arial, sans-serif; margin: 2em; }}
        table {{ border-collapse: collapse; margin: 1em 0; }}
        th, td {{ border: 1px solid #ccc; padding: 0.5em 1em; text-align: left; }}
        th {{ background: #f4f4f4; }}
        h1 {{ color: #333; }}
        .section {{ margin: 1.5em 0; }}
      </style>
    </head>
    <body>
      <h1>Cross-Species Contamination — Pipeline Report</h1>
      <p>Generated: {timestamp}</p>

      <div class="section">
        <h2>Summary</h2>
        <table>
          <tr><th>Metric</th><th>Value</th></tr>
          <tr><td>Samples processed</td><td>{n_samples}</td></tr>
          <tr><td>Taxa in matrix</td><td>{n_taxa}</td></tr>
          <tr><td>Detection method</td><td>{method}</td></tr>
          <tr><td>Flagged (sample, taxon) pairs</td><td>{flagged_count}</td></tr>
          <tr><td>Quarantined samples</td><td>{n_quarantined}</td></tr>
        </table>
      </div>

      <div class="section">
        <h2>Quarantined Samples</h2>
        <ul>{quarantine_items}</ul>
      </div>

      <div class="section">
        <h2>Output Files</h2>
        <table>
          <tr><th>Stage</th><th>Path</th></tr>
          <tr><td>Aggregated matrix (CPM)</td><td>aggregate/taxa_matrix_cpm.tsv</td></tr>
          <tr><td>Aggregated matrix (raw)</td><td>aggregate/taxa_matrix_raw.tsv</td></tr>
          <tr><td>Aggregation metadata</td><td>aggregate/aggregation_metadata.json</td></tr>
          <tr><td>Flagged samples</td><td>detect/flagged_samples.tsv</td></tr>
          <tr><td>QC summary</td><td>detect/qc_summary.json</td></tr>
          <tr><td>Quarantine list</td><td>detect/quarantine_list.txt</td></tr>
          <tr><td>MultiQC data</td><td>csc_pipeline_mqc.yaml</td></tr>
        </table>
      </div>
    </body>
    </html>\"\"\"

    Path("pipeline_report.html").write_text(report)
    """
}
