#!/usr/bin/env nextflow

// This script is part of the Bifrost pipeline. Please see
// the accompanying LICENSE document for licensing issues,
// and the WIKI for this repo for instructions.

// Which version do we have?
//if (workflow.commitId) {
  version = "v0.2.0 ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]"
//}
//else {
//  version = "v0.2.0 local"
//}

log.info ''
log.info "================================================="
log.info " Bifrost quality assessment module, version ${version}"
log.info "================================================="
log.info "Reads                   : ${params.reads}"
log.info "#files in read set      : ${params.setsize}"
log.info "Results can be found in : ${params.out_dir}"
log.info "================================================="
log.info ""

// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set{fastqc_reads}

// Second is to send all of through fastqc

process run_fastqc {
    publishDir "${params.out_dir}/${params.fastqc}", mode: "${params.savemode}"
    tag {pair_id}
    label 'one'

    input:
    set pair_id, file(reads) from fastqc_reads

    output:
    file "$pair_id" into fastqc_results

    """
    mkdir ${pair_id}
    fastqc -q ${reads} -o ${pair_id} -t $task.cpus
    """
}

process run_multiqc {
    publishDir "${params.out_dir}/multiqc", mode: "${params.savemode}"
    tag {"multiqc"}
    label 'one'

    input:
    file "fastqc_output/*" from fastqc_results.toSortedList()

    output:
    file "multiqc_report.html" into multiqc_report

    """
    multiqc fastqc_output
    """
}
