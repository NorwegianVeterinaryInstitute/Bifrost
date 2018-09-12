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
    .into{fastqc_reads; read_pairs}

// run_fastq and run_multiqc are exactly the same as qc_track
process run_fastqc {
    publishDir "${params.out_dir}/fastqc", mode: 'copy'
    tag { pair_id }
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
    publishDir "${params.out_dir}/multiqc", mode: 'copy'
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

// if there are more than two data files, we need to cat them together
// because spades becomes complicated with more than two files
process collate_data {
    // Note, not publishing these because that would mean
    // triple copies of the files on the system
    tag {pair_id}
    label 'one'

    input:
    set pair_id, file(reads) from read_pairs

    output:
    set pair_id, file("${pair_id}*_concat.fq.gz") into (reads, pilon_reads)

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

	publishDir "${params.out_dir}/bbduk", mode: "copy"
	tag { pair_id }

	input:
	set pair_id, file(reads) from reads

    output:
    set pair_id, file("${pair_id}*_concat_stripped.fq.gz") into reads_stripped
    file "${pair_id}_bbduk_output.log"

    """
    ${preCmd}
    bbduk.sh threads=$task.cpus ref=${params.stripgenome} \
     in1=${pair_id}_R1_concat.fq.gz \
     in2=${pair_id}_R2_concat.fq.gz \
     outm=${pair_id}_matched.fq.gz \
     out1=${pair_id}_R1_concat_stripped.fq.gz \
     out2=${pair_id}_R2_concat_stripped.fq.gz \
     k=31 hdist=1 stats=stats.txt &> ${pair_id}_bbduk_output.log
    """
}


/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process run_trim {
	publishDir "${params.out_dir}/bbduk_trimmed", mode: "copy"
	tag { pair_id }

    input:
    set pair_id, file(reads) from reads_stripped

    output:
    set pair_id, file("${pair_id}*_concat_stripped_trimmed.fq.gz") into reads_trimmed
    file "${pair_id}_concat_stripped_trimmed.log"

    """
    ${preCmd}
    trimmomatic PE -threads $task.cpus -trimlog ${pair_id}_concat_stripped_trimmed.log ${pair_id}*_concat_stripped.fq.gz \
        -baseout ${pair_id}_trimmed ILLUMINACLIP:${params.adapter_dir}/${params.adapters}:${params.illuminaClipOptions} \
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
process run_spadesasm {
	publishDir "${params.out_dir}/spades", mode: "copy"
	tag { pair_id }
  label 'longtime'

	input:
	set pair_id, file(reads) from reads_trimmed

	output:
	set pair_id, file("${pair_id}_spades_scaffolds.fasta") \
    into (assembly_results, tobwa_results)
	file "${pair_id}_spades.log"

	"""
	${preCmd}
	spades.py ${params.careful} --cov-cutoff=${params.cov_cutoff} \
	    -1 ${pair_id}_R1_concat_stripped_trimmed.fq.gz \
	    -2 ${pair_id}_R2_concat_stripped_trimmed.fq.gz \
	    -s ${pair_id}_S_concat_stripped_trimmed.fq.gz -t $task.cpus -o ${pair_id}_spades
	cp ${pair_id}_spades/scaffolds.fasta ${pair_id}_spades_scaffolds.fasta
	cp ${pair_id}_spades/spades.log ${pair_id}_spades.log
	"""
}

// integrate pilon. I need to have a mapping step, followed by a pilon
// step.

/*
 * Map reads to the spades assembly
 */
process run_bwamem {
	publishDir "${params.out_dir}/bwamem", mode: "copy"
	tag { pair_id }
  label 'longtime'

	input:
	set pair_id, file("${pair_id}_spades_scaffolds.fasta") from tobwa_results
  set pair_id, file(reads) from pilon_reads

  output:
	set pair_id, file("${pair_id}_mapped_sorted.bam"), \
    file("${pair_id}_mapped_sorted.bam.bai") into bwamem_results

  """
  bwa index ${pair_id}_spades_scaffolds.fasta
  bwa mem -t $task.cpus  ${pair_id}_spades_scaffolds.fasta \
  *.fq.gz | samtools sort -o ${pair_id}_mapped_sorted.bam -
  samtools index ${pair_id}_mapped_sorted.bam
  """
}

/*
* Incorporating pilon_reads
*/

process run_pilon {
  publishDir "${params.out_dir}/pilon", mode: "copy"
	tag { pair_id }

	input:
	set pair_id, file("${pair_id}_mapped_sorted.bam"), \
    file("${pair_id}_mapped_sorted.bam.bai") from bwamem_results
  set pair_id, file("${pair_id}_spades_scaffolds.fasta") from assembly_results

  output:
	set pair_id, file("${pair_id}_pilon_spades.*") into pilon_results
  set pair_id, file("${pair_id}_pilon_spades.fasta") into to_prokka
  file "${pair_id}_pilon_spades.fasta" into asms_for_quast

  """
  pilon --threads $task.cpus --genome ${pair_id}_spades_scaffolds.fasta \
  --bam ${pair_id}_mapped_sorted.bam --output ${pair_id}_pilon_spades \
  --changes --vcfqe &> ${pair_id}_pilon_spades.log
  """

}

/*
* Annotation using PROKKA
*/
process run_prokka {
	publishDir "${params.out_dir}/prokka", mode: "copy"
	tag { pair_id }

	input:
	set pair_id, file("${pair_id}_pilon_spades.fasta") from to_prokka

	output:
	set pair_id, file("${pair_id}.*") into annotation_results

  """
	${preCmd}
  prokka --compliant --force --usegenus --cpus $task.cpus \
  --centre ${params.centre} --prefix ${pair_id} --locustag ${params.locustag} \
  --genus ${params.genus} --species ${params.species} \
  --kingdom ${params.kingdom} ${params.prokka_additional} \
  --outdir . ${pair_id}_pilon_spades.fasta
  """

}

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
	quast --threads $task.cpus -o quast_evaluation_all \
		-G ${params.quast_genes} -R ${params.quast_ref} \
	    ${asm_list} \
	"""
}
