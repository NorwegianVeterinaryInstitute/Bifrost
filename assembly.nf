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

// Returns a tuple of read pairs in the form
// [dataset_id, forward.fq, reverse.fq] where
// the dataset_id is the shared prefix from
// the two paired FASTQ files.
Channel
	.fromFilePairs(params.read_pairs, flat: true)
	.ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
	.into { fastq_read_pairs }

/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 * At present only seems to care about paired reads, not singleton output
 */
process RunQC {
	publishDir "${params.out_dir}/PreProcessing", mode: "copy"

	tag { dataset_id }

        input:
        set dataset_id, file(forward), file(reverse) from fastq_read_pairs

        output:
        set dataset_id, file("${dataset_id}_1P.fastq"), file("${dataset_id}_2P.fastq") into (velvet_read_pairs, spades_read_pairs)

        """
        java -jar ${HOME}/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads ${threads} $forward $reverse -baseout ${dataset_id} ILLUMINACLIP:${HOME}/bin/Trimmomatic-0.36/adapters/${adapters}:2:30:10:3:TRUE LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}
        mv ${dataset_id}_1P ${dataset_id}_1P.fastq
        mv ${dataset_id}_2P ${dataset_id}_2P.fastq
        """
}


/*
 * Build assembly with Velvet
 */
//process BuildVelvetAssembly {
//	publishDir "${params.out_dir}/VelvetContigs", mode: "copy"
//
//	tag { dataset_id }
//
//	input:
//	set dataset_id, file(forward), file(reverse) from velvet_kg_pairs
//	val best from best_velvet_kmer_results
//
//	output:
//	set dataset_id, file("${dataset_id}_velvet-contigs.fa") into (velvet_assembly_results, velvet_assembly_quast_contigs)
//
//
//	shell:
//	'''
//	#!/bin/sh
//	best_kmer=`cat !{best}`
//	velveth auto $best_kmer -separate -fastq -shortPaired !{forward} !{reverse}
//	velvetg auto -exp_cov auto -cov_cutoff auto
//	mv auto/contigs.fa !{dataset_id}_velvet-contigs.fa
//	'''
//}

/*
 * Build assembly with SPAdes
 */
process BuildSpadesAssembly {
	publishDir "${params.out_dir}/SPadesContigs", mode: "copy"

	tag { dataset_id }

	input:
	set dataset_id, file(forward), file(reverse) from spades_read_pairs

	output:
	set dataset_id, file("${dataset_id}_spades-contigs.fa") into (spades_assembly_results, spades_assembly_quast_contigs)

	"""
	${HOME}/bin/SPAdes-3.10.1-Linux/bin/spades.py --pe1-1 ${forward} --pe1-2 ${reverse} -t ${threads} -o spades_output
	mv spades_output/contigs.fasta ${dataset_id}_spades-contigs.fa
	"""
}




//process AnnotateContigs {
//	publishDir "${params.out_dir}/AnnotatedContigs", mode: "copy"
//
//	tag { dataset_id }
//
//	input:
//	set dataset_id, file(cisa_contigs) from cisa_integrated_contigs
//
//	output:
//	file("${dataset_id}.*") into prokka_annotations
//
//	shell:
//	'''
//	#!/bin/sh
//	if [ !{species} && !{genus} ]
//	then
//		prokka !{cisa_contigs} --genus !{genus} --species !{species} --centre tychus --prefix !{dataset_id} --cpus !{threads} --outdir annotations
//	else
//		prokka !{cisa_contigs} --prefix !{dataset_id} --cpus !{threads} --outdir annotations
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
//	tag { dataset_id }
//
//	input:
//	set dataset_id, file(quast_contigs) from grouped_assembly_quast_contigs
//
//	output:
//	file("${dataset_id}_*") into quast_evaluation
//
//	shell:
//	'''
//	#!/bin/sh
//	quast.py !{quast_contigs[0]} !{quast_contigs[1]} !{quast_contigs[2]} !{quast_contigs[3]} !{quast_contigs[4]} --space-efficient --threads !{threads} -o output
//        mkdir quast_output
//        find output/ -maxdepth 2 -type f | xargs mv -t quast_output
//        cd quast_output
//        ls * | xargs -I {} mv {} !{dataset_id}_{}
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