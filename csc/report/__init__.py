"""Static HTML contamination summary report module.

Consumes outputs of :mod:`csc.aggregate` (and optionally :mod:`csc.detect`)
and renders a single, self-contained HTML file summarising non-human
content in a cohort of human WGS samples.  See :mod:`csc.report.report`
for the rendering logic and :mod:`csc.report.cli` for the ``csc-report``
command-line interface.

The report deliberately uses only the Python standard library so that it
can be embedded in any pipeline without adding heavyweight dependencies
(Jinja2, pandas, matplotlib, etc.).  All figures are emitted as inline
SVG so the HTML is a single self-contained artefact.

AI assistance acknowledgment: This module was developed with AI
assistance.  Best practices in the bioinformatics field should always
take precedence over specific implementation details.
"""

from __future__ import annotations

from csc.report.report import (
    REPORT_SCHEMA_VERSION,
    ReportInputs,
    generate_html_report,
    load_inputs,
)

__all__ = [
    "REPORT_SCHEMA_VERSION",
    "ReportInputs",
    "generate_html_report",
    "load_inputs",
]
