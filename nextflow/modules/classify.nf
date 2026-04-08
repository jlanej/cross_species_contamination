/*
 * classify.nf – Process definition for taxonomic classification
 *
 * Wraps the ``csc-classify`` CLI to run Kraken2 on extracted FASTQ
 * files and emit per-sample classification reports.
 *
 * AI assistance acknowledgment: developed with AI assistance.
 */

process CLASSIFY_READS {
    tag "${sample_id}"
    publishDir "${params.outdir}/classify/${sample_id}", mode: 'copy'

    cpus   params.classify_cpus   ?: 4
    memory params.classify_memory ?: '8 GB'
    time   params.classify_time   ?: '4h'

    input:
    tuple val(sample_id), path(fastqs)

    output:
    tuple val(sample_id), path("*.kraken2.report.txt"), emit: reports
    tuple val(sample_id), path("*.kraken2.output.txt"), emit: outputs

    script:
    def fq_list     = fastqs instanceof List ? fastqs : [fastqs]
    def is_paired   = fq_list.size() == 2
    def paired_arg  = is_paired ? '--paired' : ''
    def input_files = fq_list.join(' ')
    def mem_arg     = params.memory_mapping ? '--memory-mapping' : ''
    """
    csc-classify \\
        ${input_files} \\
        --db ${params.kraken2_db} \\
        -o . \\
        --sample-id ${sample_id} \\
        --confidence ${params.confidence} \\
        --threads ${task.cpus} \\
        ${paired_arg} \\
        ${mem_arg}
    """
}
