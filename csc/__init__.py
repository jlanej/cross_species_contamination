"""Cross-Species Contamination Detection Pipeline.

A multi-module pipeline for detecting cross-species contamination in
whole-genome sequencing data.

Modules
-------
extract
    Streaming extraction of unmapped/poorly-mapped reads from BAM/CRAM files.
classify
    Taxonomic classification of extracted reads (stub).
aggregate
    Aggregation of classification results across samples (stub).
detect
    Statistical detection of cross-species contamination (stub).
"""

__version__ = "0.2.0"
