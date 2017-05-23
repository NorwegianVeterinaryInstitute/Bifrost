#!/usr/bin/env nextflow

// Adapted from https://github.com/cdeanj/nextflow-tychus
// should put all of these on top into a config file.


// Display help menu
if(params.help) {
	log.info ''
	log.info 'Bifrost - Assembly Pipeline'
	log.info ''
	log.info 'Usage: '
	log.info '    nextflow assembly.nf -profile assembly [options]'
	log.info ''
	log.info 'General Options: '
	log.info '    --read_pairs      DIR		Directory of paired FASTQ files'
	log.info '    --threads         INT		Number of threads to use for each process'
	log.info '    --output          DIR		Directory to write output files to'
	log.info ''
	log.info 'Trimmomatic Options: '
	log.info '    --leading         INT		Remove leading low quality or N bases'
	log.info '    --trailing        INT		Remove trailing low quality or N bases'
	log.info '    --slidingwindow   INT		Scan read with a sliding window'
	log.info '    --minlen          INT		Drop reads below INT bases long'
	log.info '    --adapters        STR		FASTA formatted adapter file'
	log.info ''
	log.info 'Prokka Options: '
	log.info '    --genus           STR		Target genus'
	log.info '    --species         STR		Target species'
	log.info ''
	return
}

preCmd = """
if [ -f /cluster/bin/jobsetup ];
then set +u; source /cluster/bin/jobsetup; set -u; fi
"""

// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ read_pairs }


// if there are more than two data files, we need to cat them together
// because spades becomes complicated with more than two files

process collate_data {
    publishDir "${params.out_dir}/${params.raw_data}", mode: 'copy'

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
    java -jar /home/karinlag/bin/trimmomatic.jar PE -threads 1 ${pair_id}_raw/*${params.file_ending} -baseout ${pair_id}_trimmed ILLUMINACLIP:/home/karinlag/bin/Trimmomatic-0.36/adapters/${params.adapters}:2:30:10:3:TRUE LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.slidingwindow} MINLEN:${params.minlen}
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
	set pair_id, file("${pair_id}_trimmed")  from spades_read_pairs

	output:
	set pair_id, file("${pair_id}_spades") into (spades_assembly_results, spades_assembly_quast_contigs)

	"""
	${preCmd}
	spades.py -1 ${pair_id}_trimmed/R1_trimmed${params.file_ending} -2 ${pair_id}_trimmed/R2_trimmed${params.file_ending} -s ${pair_id}_trimmed/single${params.file_ending} -t 1 -o ${pair_id}_spades
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
//process EvaluateAssemblies {
//	publishDir "${params.out_dir}/AssemblyReport", mode: "move"
//
//	tag { pair_id }
//
//	input:
//	set pair_id, file(quast_contigs) from grouped_assembly_quast_contigs
//
//	output:
//	file("${pair_id}_*") into quast_evaluation
//
//	shell:
//	'''
//	#!/bin/sh
//	quast.py !{quast_contigs[0]} !{quast_contigs[1]} !{quast_contigs[2]} !{quast_contigs[3]} !{quast_contigs[4]} --space-efficient --threads !{threads} -o output
//        mkdir quast_output
//        find output/ -maxdepth 2 -type f | xargs mv -t quast_output
//        cd quast_output
//        ls * | xargs -I {} mv {} !{pair_id}_{}
//        mv * ../
//	'''
//}

// Display information about the completed run
// See https://www.nextflow.io/docs/latest/metadata.html for more
// information about available onComplete options
workflow.onComplete {
	log.info "Nextflow Version:	$workflow.nextflow.version"
  	log.info "Command Line:		$workflow.commandLine"
	log.info "Container:		$workflow.container"
	log.info "Duration:		$workflow.duration"
	log.info "Output Directory:	$params.out_dir"
}