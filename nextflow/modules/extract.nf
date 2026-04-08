/*
 * extract.nf – Process definition for unmapped-read extraction
 *
 * Wraps the ``csc-extract`` CLI to pull unmapped (and optionally
 * low-MAPQ) reads from BAM/CRAM inputs and emit FASTQ files.
 *
 * AI assistance acknowledgment: developed with AI assistance.
 */

process EXTRACT_UNMAPPED {
    tag "${sample_id}"
    publishDir "${params.outdir}/extract/${sample_id}", mode: 'copy'

    cpus   params.extract_cpus   ?: 4
    memory params.extract_memory ?: '4 GB'
    time   params.extract_time   ?: '2h'

    input:
    tuple val(sample_id), path(input_file), val(reference)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: fastqs

    script:
    def ref_arg  = reference ? "--reference ${reference}" : ''
    def mapq_arg = params.mapq ? "--mapq ${params.mapq}" : ''
    """
    csc-extract \\
        ${input_file} \\
        -o . \\
        --sample-id ${sample_id} \\
        --threads ${task.cpus} \\
        ${ref_arg} \\
        ${mapq_arg}
    """
}
