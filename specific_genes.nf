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
log.info "MLST Scheme used        : ${params.mlst_scheme}"
log.info "AMR database            : ${params.amr_db}"
log.info "Virulence db            : ${params.vir_db}"
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
    .set{read_pairs}


// if there are more than two data files, we need to cat them together
// because ariba does not permit us to run with more than two files

process collate_data {
    // Note, not publishing these because that would mean
    // triple copies of the files on the system

    input:
    set pair_id, file(reads) from read_pairs

    output:
    //set pair_id, file("${pair_id}_R{1,2}${params.file_ending}") into read_pairs_mlst, read_pairs_amr, read_pairs_vir
    set pair_id, file("${pair_id}*_concat.fq.gz") into \
      (read_pairs_mlst, read_pairs_amr, read_pairs_vir)

    """
    ${preCmd}
    cat ${pair_id}*R1* > ${pair_id}_R1_concat.fq.gz
    cat ${pair_id}*R2* > ${pair_id}_R2_concat.fq.gz
    """
}


// The following three processes are for MLST finding

// Download database database
process run_ariba_mlst_prep {
    publishDir params.out_dir + "/" + params.mlst_results, mode: 'copy'

    output:
    file "mlst_db" into mlst_db

    when:
    params.do_mlst == "yes"

    """
    ${preCmd}
    ariba pubmlstget "${params.mlst_scheme}" mlst_db
    """
}

// Run ariba on each dataset

process run_ariba_mlst_pred {
    //tag {$pair_id}
    publishDir params.out_dir + "/" + params.mlst_results, mode: 'copy'

    input:
    set pair_id, file(reads) from read_pairs_mlst
    file "mlst_db" from mlst_db

    output:
    file "${pair_id}_mlst_report.tsv" into pair_id_mlst_tsv
    file "${pair_id}_ariba" into pair_id_mlst_aribadir

    when:
    params.do_mlst == "yes"

    """
    ${preCmd}
    ariba run --threads $task.cpus mlst_db/ref_db ${pair_id}_R*_concat.fq.gz ${pair_id}_ariba &> ariba.out
    echo -e "header\t" \$(head -1 ${pair_id}_ariba/mlst_report.tsv) > ${pair_id}_mlst_report.tsv
    echo -e "${pair_id}\t" \$(tail -1 ${pair_id}_ariba/mlst_report.tsv) >> ${pair_id}_mlst_report.tsv
    """
}

// Summarize MLST results
process run_ariba_mlst_summarize {
    publishDir params.out_dir + "/" + params.mlst_results, mode: 'copy'

    input:
    file pair_id_mlst_tsv from pair_id_mlst_tsv.collect()

    output:
    file "mlst_summarized_results.tsv" into mlst_summarized

    when:
    params.do_mlst == "yes"

    """
    ${preCmd}
    cat ${pair_id_mlst_tsv} >> mlst_summarized_results_tmp.tsv
    head -1 mlst_summarized_results_tmp.tsv > mlst_summarized_results.tsv
    cat mlst_summarized_results_tmp.tsv | grep -v "ST" >> mlst_summarized_results.tsv
    """
}



//  These three processes are for AMR prediction
process run_ariba_amr_prep {
    publishDir params.out_dir + "/" + params.amr_results, mode: 'copy'

    output:
    file "db_amr_prepareref" into db_amr_prepareref

    when:
    params.do_amr == "yes"

    """
    ${preCmd}
    ariba getref ${params.amr_db} amr_db
    ariba prepareref -f amr_db.fa -m amr_db.tsv db_amr_prepareref
    """
}

process run_ariba_amr_pred {
    publishDir params.out_dir + "/" + params.amr_results, mode: 'copy'

    input:
    set pair_id, file(reads) from read_pairs_amr
    file "db_amr_prepareref" from db_amr_prepareref

    output:
    file "${pair_id}_amr_report.tsv" into pair_id_amr_tsv
    file "${pair_id}_ariba" into pair_id_amr_aribadir

    when:
    params.do_amr == "yes"


    """
    ${preCmd}
    ariba run --threads $task.cpus db_amr_prepareref ${pair_id}_R*_concat.fq.gz ${pair_id}_ariba &> ariba.out
    cp ${pair_id}_ariba/report.tsv ${pair_id}_amr_report.tsv

    """
}

// Summarize AMR results
process run_ariba_amr_summarize {
    publishDir params.out_dir + "/" + params.amr_results, mode: 'copy'

    input:
    file pair_id_amr_tsv from pair_id_amr_tsv.collect()

    output:
    file "amr_summarized*" into amr_summarized

    when:
    params.do_amr == "yes"

    """
    ${preCmd}
    ariba summary amr_summarized ${pair_id_amr_tsv}
    """
}

//  These three processes are for virulence prediction
process run_ariba_vir_prep {
    publishDir params.out_dir + "/" + params.vir_results, mode: 'copy'

    output:
    file "db_vir_prepareref" into db_vir_prepareref

    when:
    params.do_vir == "yes"

    """
    ${preCmd}
    ariba getref ${params.vir_db} vir_db
    ariba prepareref -f vir_db.fa -m vir_db.tsv db_vir_prepareref
    """
}

process run_ariba_vir_pred {
    publishDir params.out_dir + "/" + params.vir_results, mode: 'copy'

    input:
    set pair_id, file(reads) from read_pairs_vir
    file "db_vir_prepareref" from db_vir_prepareref

    output:
    file "${pair_id}_vir_report.tsv" into pair_id_vir_tsv
    file "${pair_id}_ariba" into pair_id_vir_aribadir

    when:
    params.do_vir == "yes"

    """
    ${preCmd}
    ariba run --threads $task.cpus db_vir_prepareref ${pair_id}_R*_concat.fq.gz \
      ${pair_id}_ariba &> ariba.out
    cp ${pair_id}_ariba/report.tsv ${pair_id}_vir_report.tsv

    """
}

// Summarize virulence results
process run_ariba_vir_summarize {
    publishDir params.out_dir + "/" + params.vir_results, mode: 'copy'

    input:
    file pair_id_vir_tsv from pair_id_vir_tsv.collect()

    output:
    file "vir_summarized*" into vir_summarized

    when:
    params.do_vir == "yes"

    """
    ${preCmd}
    ariba summary vir_summarized ${pair_id_vir_tsv}
    """
}

// Display information about the completed run
// See https://www.nextflow.io/docs/latest/metadata.html for more
// information about available onComplete options
workflow.onComplete {
	log.info "Nextflow Version:	$workflow.nextflow.version"
  log.info "Command Line:		$workflow.commandLine"
	log.info "Container:		$workflow.container"
	log.info "Duration:		    $workflow.duration"
	log.info "Output Directory:	$params.out_dir"
}
