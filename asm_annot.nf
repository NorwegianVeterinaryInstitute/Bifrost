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

// TODO I need to incorporate some options for prokka here
// Also strip genome

log.info ''
log.info "================================================="
log.info " Bifrost assembly module version ${version}"
log.info "================================================="
log.info "Reads                   : ${params.reads}"
log.info "#files in read set      : ${params.setsize}"
log.info "Results can be found in : ${params.out_dir}"
log.info "================================================="
log.info ""

// Needed to run on the Abel cluster
preCmd = """
if [ -f /cluster/bin/jobsetup ];
then set +u; source /cluster/bin/jobsetup; set -u; fi
"""

// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ fastqc_reads; read_pairs }

// run_fastq and run_multiqc are exactly the same as qc_track
process run_fastqc {
    publishDir "${params.out_dir}/fastqc", mode: 'copy'

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

process run_multiqc {
    publishDir "${params.out_dir}/multiqc", mode: 'copy'
    tag {"multiqc"}

    input:
    file "fastqc_output/*" from fastqc_results.toSortedList()

    output:
    file "multiqc_report.html" into multiqc_report

    """
    multiqc fastqc_output
    """
}

// if there are more than two data files, we need to cat them together
// because spades becomes complicated with more than two files
process collate_data {
    // Note, not publishing these because that would mean
    // triple copies of the files on the system

    tag {pair_id}

    input:
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, file("${pair_id}*_concat.fq.gz") into reads

    """
    ${preCmd}
    cat ${pair_id}*R1* > ${pair_id}_R1_concat.fq.gz
    cat ${pair_id}*R2* > ${pair_id}_R2_concat.fq.gz
    """
}


/*
 * Strip PhiX with bbmap
 */
process run_strip {

	publishDir "${params.out_dir}/bbmap", mode: "copy"

	tag { pair_id }

	input:
	set pair_id, file(reads) from reads

    output:
    set pair_id, file("${pair_id}*_concat_stripped.fq.gz") into reads_stripped
    file "${pair_id}_concat_mapped.sam" into mapped_reads
    file "bbmap_output.log"

    """
    ${preCmd}
    $task.bbmap threads=$task.threads ref=${params.stripgenome} path=${params.stripdir} \
    covstats=constats.txt covhist=covhist.txt basecov=basecov.txt bincov=bincov.txt \
     in=${pair_id}_R1_concat.fq.gz \
     in2=${pair_id}_R2_concat.fq.gz \
     out=${pair_id}_concat_mapped.sam \
     outu=${pair_id}_R1_concat_stripped.fq.gz \
     outu2=${pair_id}_R2_concat_stripped.fq.gz > bbmap_output.log
    """
}


/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process run_trim {
	publishDir "${params.out_dir}/bbmap_trimmed", mode: "copy"

	tag { pair_id }

    input:
    set pair_id, file(reads) from reads_stripped

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
	publishDir "${params.out_dir}/spades", mode: "copy"

	tag { pair_id }

	input:
	set pair_id, file(reads) from reads_trimmed

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
  // The output here is a directory in and of itself
  // thus not creating a new one
	publishDir "${params.out_dir}/", mode: "copy"

	tag { pair_id }

	input:
	file asm_list from asms_for_quast.toSortedList()

  //TODO: fix this, is why output is not going anywhere
	output:
	file quast_evaluation_all into quast_evaluation_all

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
