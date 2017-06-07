#!/usr/bin/env nextflow

// Adapted from https://github.com/cdeanj/nextflow-tychus
// should put all of these on top into a config file.

version = "something"

log.info ''
log.info "================================================="
log.info " Bifrost assembly moduoe v${version}"
log.info "================================================="
log.info "Reads                   : ${params.reads}"
log.info "#files in read set      : ${params.setsize}"
log.info "Quast reference         : ${params.quast_ref}"
log.info "Quast gene set          : ${params.quast_genes}"
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
    $task.fastqc -q ${reads} -o ${pair_id}
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
    $task.python bin/fastqc_eval.py -d fastqc_output -o fastqc_eval.results
    """
}

// if there are more than two data files, we need to cat them together
// because spades becomes complicated with more than two files
process collate_data {
    publishDir "${params.out_dir}/${params.raw_data}", mode: 'copy'

    tag { pair_id }

    input:
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, file("${pair_id}_raw") into reads

    """
    ${preCmd}
    mkdir ${pair_id}_raw
    cat ${pair_id}*R1*${params.file_ending} > ${pair_id}_raw/R1${params.file_ending}
    cat ${pair_id}*R2*${params.file_ending} > ${pair_id}_raw/R2${params.file_ending}
    """
}

/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process run_trim {
	publishDir "${params.out_dir}/${params.trimmed}", mode: "copy"

	tag { pair_id }

    input:
    set pair_id, file("${pair_id}_raw") from reads

    output:
    set pair_id, file("${pair_id}_trimmed") into spades_read_pairs

    """
    ${preCmd}
    mkdir ${pair_id}_trimmed
    $task.trimmomatic PE -threads 1 ${pair_id}_raw/*${params.file_ending} \
        -baseout ${pair_id}_trimmed ILLUMINACLIP:$task.adapter_dir/${params.adapters}:2:30:10:3:TRUE \
        LEADING:${params.leading} TRAILING:${params.trailing} \
        SLIDINGWINDOW:${params.slidingwindow} MINLEN:${params.minlen}
    mv ${pair_id}_trimmed_1P ${pair_id}_trimmed/R1_trimmed${params.file_ending}
    mv ${pair_id}_trimmed_2P ${pair_id}_trimmed/R2_trimmed${params.file_ending}
    cat ${pair_id}_trimmed_1U ${pair_id}_trimmed_2U > ${pair_id}_trimmed/single${params.file_ending}
    """
}


/*
 * Build assembly with SPAdes
 */
process spades_assembly {
	publishDir "${params.out_dir}/${params.assembly}", mode: "copy"

	tag { pair_id }

	input:
	set pair_id, file("${pair_id}_trimmed") from spades_read_pairs

	output:
	set pair_id, file("${pair_id}_spades") into spades_assembly_results
	file("${pair_id}_spades_scaffolds.fasta") into asms_for_quast

	"""
	${preCmd}
	$task.spades -1 ${pair_id}_trimmed/R1_trimmed${params.file_ending} \
	    -2 ${pair_id}_trimmed/R2_trimmed${params.file_ending} \
	    -s ${pair_id}_trimmed/single${params.file_ending} -t $task.cpus -o ${pair_id}_spades
	cp ${pair_id}_spades/scaffolds.fasta ${pair_id}_spades_scaffolds.fasta
	"""
}


//process AnnotateContigs {
//	publishDir "${params.out_dir}/AnnotatedContigs", mode: "copy"
//
//	tag { pair_id }
//
//	input:
//	set pair_id, file(cisa_contigs) from cisa_integrated_contigs
//
//	output:
//	file("${pair_id}.*") into prokka_annotations
//
//	shell:
//	'''
//	#!/bin/sh
//	if [ !{species} && !{genus} ]
//	then
//		prokka !{cisa_contigs} --genus !{genus} --species !{species} --centre tychus --prefix !{pair_id} --cpus !{threads} --outdir annotations
//	else
//		prokka !{cisa_contigs} --prefix !{pair_id} --cpus !{threads} --outdir annotations
//	fi
//	mv annotations/* .
//	'''
//}

//abyss_assembly_quast_contigs.concat(
//                velvet_assembly_quast_contigs,
//               spades_assembly_quast_contigs,
//                idba_assembly_quast_contigs,
//		cisa_integrated_quast_contigs
//       )
//        .groupTuple(sort: true, size: 5)
//        .into { grouped_assembly_quast_contigs }


/*
 * Evaluate ALL assemblies with QUAST
 */
process quast_eval {
	publishDir "${params.out_dir}/${params.quast}", mode: "copy"

	tag { pair_id }

	input:
	file asm_list from asms_for_quast.toSortedList()

	output:
	file "quast_evaluation_all"  into quast_evaluation_all

	"""
	${preCmd}
	$task.quast --threads $task.threads -o quast_evaluation_all \
		-G ${params.quast_genes} -R ${params.quast_ref} \
	    ${asm_list} \
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