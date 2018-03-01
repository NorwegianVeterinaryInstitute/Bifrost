#!/usr/bin/env nextflow

// Note on coding: with Nextflow, we can use wildcards when specifying
// output, but not input. Thus I'm using a dummy variable (dummyvar)
// in some places to get things shipped from one process to another.
// This eans that I'm controlling what files are shipped from one
// process to the other via the output statement, not the input
// statement.


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
    set pair_id, file(reads) from fastqc_reads//.view()

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    """
    mkdir ${pair_id}
    $task.fastqc -q ${reads} -o . -t $task.threads
    """
}

// TODO: need to fix the python script so that it evals a set of files
// or something

//process run_fastqc_eval {
//    publishDir "${params.out_dir}/${params.fastqc_eval}", mode: 'copy'
//
//    tag { pair_id }
//
//    input:
//    file "*" from fastqc_results.toSortedList()//.view()
//
//    output:
//    file "fastqc_eval.results"
//
//    """
//    fastqc_eval.py -d fastqc_output -o fastqc_eval.results
//    """
//}

// if there are more than two data files, we need to cat them together
// because spades becomes complicated with more than two files
process collate_data {
    publishDir "${params.out_dir}/${params.raw_data}", mode: 'copy'

    tag { pair_id }

    input:
    set pair_id, file(reads) from read_pairs//.view()

    output:
    set pair_id, file("${pair_id}*_concat.fq.gz") into reads

    """
    ${preCmd}
    cat ${pair_id}*R1*${params.file_ending} > ${pair_id}_R1_concat.fq.gz
    cat ${pair_id}*R2*${params.file_ending} > ${pair_id}_R2_concat.fq.gz
    """
}


/*
 * Strip PhiX with bbmap
 */
process run_strip {

	publishDir "${params.out_dir}/${params.bbmap}", mode: "copy"

	tag { pair_id }

	input:
	set pair_id, file(dummyvar) from reads.view()

    output:
    set pair_id, file("${pair_id}*_concat_stripped.fq.gz") into reads_stripped
    file "${pair_id}_concat_mapped.sam" into mapped_reads

    """
    ${preCmd}
    $task.bbmap threads=$task.threads ref=${params.stripgenome} path=${params.stripdir} \
     in=${pair_id}_R1_concat.fq.gz \
     in2=${pair_id}_R2_concat.fq.gz \
     out=${pair_id}_concat_mapped.sam \
     outu=${pair_id}_R1_concat_stripped.fq.gz \
     outu2=${pair_id}_R2_concat_stripped.fq.gz
    """
}


/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process run_trim {
	publishDir "${params.out_dir}/${params.trimmed}", mode: "copy"

	tag { pair_id }

    input:
    set pair_id, file(dummyvar) from reads_stripped.view()

    output:
    set pair_id, file("${pair_id}*_concat_stripped_trimmed.fq.gz") into reads_trimmed
    file "${pair_id}_concat_stripped_trimmed.log"

    """
    ${preCmd}
    $task.trimmomatic PE -threads $task.threads -trimlog ${pair_id}_concat_stripped_trimmed.log ${pair_id}*_concat_stripped.fq.gz \
        -baseout ${pair_id}_trimmed ILLUMINACLIP:$task.adapter_dir/${params.adapters}:${params.illuminaClipOptions} \
        SLIDINGWINDOW:${params.slidingwindow} \
        LEADING:${params.leading} TRAILING:${params.trailing} \
        MINLEN:${params.minlen} &> ${pair_id}_run.log
    mv ${pair_id}_trimmed_1P ${pair_id}_R1_concat_stripped_trimmed.fq.gz
    mv ${pair_id}_trimmed_2P ${pair_id}_R2_concat_stripped_trimmed.fq.gz
    cat ${pair_id}_trimmed_1U ${pair_id}_trimmed_2U > ${pair_id}_S_concat_stripped_trimmed.fq.gz
    """
}


/*
 * Build assembly with SPAdes
 */
process spades_assembly {
	publishDir "${params.out_dir}/${params.assembly}", mode: "copy"

	tag { pair_id }

	input:
	set pair_id, file(dummyvar) from reads_trimmed.view()

	output:
	set pair_id, file("${pair_id}_spades_scaffolds.fasta") into assembly_results
	file "${pair_id}_spades_scaffolds.fasta" into asms_for_quast
	file "${pair_id}_spades.log"

	"""
	${preCmd}
	$task.spades ${params.careful} --cov-cutoff=${params.cov_cutoff} \
	    -1 ${pair_id}_R1_concat_stripped_trimmed.fq.gz \
	    -2 ${pair_id}_R2_concat_stripped_trimmed.fq.gz \
	    -s ${pair_id}_S_concat_stripped_trimmed.fq.gz -t $task.cpus -o ${pair_id}_spades
	cp ${pair_id}_spades/scaffolds.fasta ${pair_id}_spades_scaffolds.fasta
	cp ${pair_id}_spades/spades.log ${pair_id}_spades.log
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
	publishDir "${params.out_dir}/", mode: "copy"

	tag { pair_id }

	input:
	file asm_list from asms_for_quast.toSortedList()

	output:
	file "${params.quast}"  into quast_evaluation_all

	"""
	${preCmd}
	$task.quast --threads $task.threads -o ${params.quast} \
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
