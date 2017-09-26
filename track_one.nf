#!/usr/bin/env nextflow

def rev = workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)

log.info ''
log.info "================================================="
log.info " Bifrost assembly module version ${rev}"
log.info "================================================="
log.info "Reads                   : ${params.reads}"
log.info "#files in read set      : ${params.setsize}"
log.info "Results can be found in : ${params.out_dir}"
log.info "================================================="
log.info ""

preCmd = """
if [ -f /cluster/bin/jobsetup ];
then set +u; source /cluster/bin/jobsetup; set -u; fi
"""

// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ fastqc_reads; read_pairs }

// Second is to send all of through fastqc

process run_fastqc {
    publishDir "${params.out_dir}/${params.fastqc}", mode: 'copy'

    tag { pair_id }

    input:
    set pair_id, file(reads) from fastqc_reads

    output:
    file "$pair_id" into fastqc_results

    """
    mkdir ${pair_id}
    $task.fastqc -q ${reads} -o ${pair_id} -t $task.threads
    """
}

process run_fastqc_eval {
    publishDir "${params.out_dir}/${params.fastqc_eval}", mode: 'copy'

    tag { pair_id }

    input:
    file "fastqc_output/*" from fastqc_results.toSortedList()

    output:
    file "fastqc_eval.results"

    """
    fastqc_eval.py -d fastqc_output -o fastqc_eval.results
    """
}