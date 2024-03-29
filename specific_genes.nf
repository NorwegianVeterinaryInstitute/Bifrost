#!/usr/bin/env nextflow

// This script is part of the Bifrost pipeline. Please see
// the accompanying LICENSE document for licensing issues,
// and the WIKI for this repo for instructions.

// Which version do we have?
if (workflow.commitId) {
  version = "v0.2.0 $workflow.revision"
}
else {
  version = "v0.2.0 local"
}

log.info "================================================="
log.info " Bifrost specific gene finding with Ariba v${version}"
log.info "================================================="
log.info "Reads                   : ${params.reads}"
log.info "#files in read set      : ${params.setsize}"
log.info "MLST Scheme used        : ${params.mlst_db}"
log.info "AMR database            : ${params.amr_db}"
log.info "Virulence db            : ${params.vir_db}"
log.info "Results can be found in : ${params.out_dir}"
log.info "================================================="
log.info ""

// First, define the input data that go into input channels
// The databases are input as value channels to enable reuse
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set{read_pairs}

mlst_db = Channel
    .value(params.mlst_db)

amr_db = Channel
    .value(params.amr_db)

vir_db = Channel
    .value(params.vir_db)


// if there are more than two data files, we need to cat them together
// because ariba does not permit us to run with more than two files

process collate_data {
    // Note, not publishing these because that would mean
    // triple copies of the files on the system

    tag {pair_id}
    label 'one'

    input:
    set pair_id, file(reads) from read_pairs

    output:
    //set pair_id, file("${pair_id}_R{1,2}${params.file_ending}") into read_pairs_mlst, read_pairs_amr, read_pairs_vir
    set pair_id, file("${pair_id}*_concat.fq.gz") into \
      (read_pairs_mlst, read_pairs_amr, read_pairs_vir)


    """
    shopt -s extglob
    cat ${pair_id}*_?(R)1[_.]*.gz > ${pair_id}_R1_concat.fq.gz
    cat ${pair_id}*_?(R)2[_.]*.gz > ${pair_id}_R2_concat.fq.gz
    """
}


// The following two processes are for MLST finding

// Run ariba on each dataset

process run_ariba_mlst_pred {
    publishDir "${params.out_dir}" + "/" + "${params.mlst_results}", mode: "${params.savemode}"
    tag {pair_id}

    input:
    set pair_id, file(reads) from read_pairs_mlst
    path mlst_db from mlst_db

    output:
    file "${pair_id}_mlst_report.tsv" into pair_id_mlst_tsv
    file "${pair_id}_ariba" into pair_id_mlst_aribadir

    when:
    params.do_mlst == "yes"

    """
    ariba run --threads $task.cpus ${mlst_db}/ref_db ${pair_id}_R*_concat.fq.gz ${pair_id}_ariba &> ariba.out
    echo -e "header\t" \$(head -1 ${pair_id}_ariba/mlst_report.tsv) > ${pair_id}_mlst_report.tsv
    echo -e "${pair_id}\t" \$(tail -1 ${pair_id}_ariba/mlst_report.tsv) >> ${pair_id}_mlst_report.tsv
    """
}

// Summarize MLST results
process run_ariba_mlst_summarize {
    publishDir "${params.out_dir}" + "/" + "${params.mlst_results}", mode: "${params.savemode}"
    tag {'Summarizing mlst'}
    label 'one'

    input:
    file pair_id_mlst_tsv from pair_id_mlst_tsv.collect()

    output:
    file "mlst_summarized_results.tsv" into mlst_summarized

    when:
    params.do_mlst == "yes"

    """
    cat ${pair_id_mlst_tsv} >> mlst_summarized_results_tmp.tsv
    head -1 mlst_summarized_results_tmp.tsv > mlst_summarized_results.tsv
    cat mlst_summarized_results_tmp.tsv | grep -v "ST" >> mlst_summarized_results.tsv
    """
}



//  These two processes are for AMR prediction
process run_ariba_amr_pred {
    publishDir "${params.out_dir}" + "/" + "${params.amr_results}", mode: "${params.savemode}"
    tag{pair_id}

    input:
    set pair_id, file(reads) from read_pairs_amr
    path db_amr_prepareref from amr_db

    output:
    file "${pair_id}_amr_report.tsv" into pair_id_amr_tsv
    file "${pair_id}_ariba" into pair_id_amr_aribadir

    when:
    params.do_amr == "yes"


    """
    ariba run --threads $task.cpus ${db_amr_prepareref} ${pair_id}_R*_concat.fq.gz ${pair_id}_ariba &> ariba.out
    cp ${pair_id}_ariba/report.tsv ${pair_id}_amr_report.tsv

    """
}

// Summarize AMR results
process run_ariba_amr_summarize {
    publishDir "${params.out_dir}" + "/" + "${params.amr_results}", mode: "${params.savemode}"
    tag{'Summarizing AMR'}
    label 'one'

    input:
    file pair_id_amr_tsv from pair_id_amr_tsv.collect()

    output:
    file "amr_summarized*" into amr_summarized

    when:
    params.do_amr == "yes"

    """
    ariba summary amr_summarized ${pair_id_amr_tsv}
    """
}

//  These two processes are for virulence prediction
process run_ariba_vir_pred {
    publishDir "${params.out_dir}" + "/" + "${params.vir_results}", mode: "${params.savemode}"
    tag{pair_id}

    input:
    set pair_id, file(reads) from read_pairs_vir
    path db_vir_prepareref from vir_db

    output:
    file "${pair_id}_vir_report.tsv" into pair_id_vir_tsv
    file "${pair_id}_ariba" into pair_id_vir_aribadir

    when:
    params.do_vir == "yes"

    """
    ariba run --threads $task.cpus ${db_vir_prepareref} ${pair_id}_R*_concat.fq.gz \
      ${pair_id}_ariba &> ariba.out
    cp ${pair_id}_ariba/report.tsv ${pair_id}_vir_report.tsv

    """
}

// Summarize virulence results
process run_ariba_vir_summarize {
    publishDir "${params.out_dir}" + "/" + "${params.vir_results}", mode: "${params.savemode}"
    tag{'Summarizing virulence'}
    label 'one'

    input:
    file pair_id_vir_tsv from pair_id_vir_tsv.collect()

    output:
    file "vir_summarized*" into vir_summarized

    when:
    params.do_vir == "yes"

    """
    ariba summary vir_summarized ${pair_id_vir_tsv}
    """
}
