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

// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{fastqc_reads; read_pairs}

// run_fastq and run_multiqc are exactly the same as qc_track
process run_fastqc {
    publishDir "${params.out_dir}/fastqc", mode: "${params.savemode}"
    tag { pair_id }
    label 'one'

    input:
    set pair_id, file(reads) from fastqc_reads

    output:
    file "$pair_id" into fastqc_multiqc

    """
    mkdir ${pair_id}
    fastqc -q ${reads} -o ${pair_id} -t $task.cpus
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
    shopt -s extglob
    cat ${pair_id}*_?(R)1[_.]*.gz > ${pair_id}_R1_concat.fq.gz
    cat ${pair_id}*_?(R)2[_.]*.gz > ${pair_id}_R2_concat.fq.gz
    """
}


/*
 * Strip PhiX with bbmap
 */
process run_strip {
    publishDir "${params.out_dir}/bbduk", 
                saveAs: {filename -> filename.endsWith('.gz') ? null:filename}, 
                mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from reads

    output:
    set pair_id, file("${pair_id}*_concat_stripped.fq.gz") into reads_stripped
    file "${pair_id}_bbduk_output.log"
    file "${pair_id}_stats.txt"
    file "${pair_id}*_stats.txt" into bbduk_stats_stripped_multiqc


    """
    bbduk.sh threads=$task.cpus ref=${params.stripgenome} \
    in1=${pair_id}_R1_concat.fq.gz \
    in2=${pair_id}_R2_concat.fq.gz \
    outm=${pair_id}_matched.fq.gz \
    out1=${pair_id}_R1_concat_stripped.fq.gz \
    out2=${pair_id}_R2_concat_stripped.fq.gz \
    k=31 hdist=1 stats=${pair_id}_stats.txt &> ${pair_id}_bbduk_output.log
    """
}


/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process run_trim {
    publishDir "${params.out_dir}/bbduk_trimmed", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from reads_stripped

    output:
    set pair_id, file("${pair_id}*_concat_stripped_trimmed.fq.gz") into (reads_trimmed, fq2)
    file "${pair_id}_stripped_trimmed_stderr.log"
    file "${pair_id}_stripped_trimmed_stdout.log"
    file "${pair_id}_stripped_trimmed_stderr.log" into bbduk_trimmed_multiqc 
    

    """
    trimmomatic PE -threads $task.cpus -trimlog ${pair_id}_concat_stripped_trimmed.log ${pair_id}*_concat_stripped.fq.gz \
    -baseout ${pair_id}_trimmed.fq.gz ILLUMINACLIP:${params.adapter_dir}/${params.adapters}:${params.illuminaClipOptions} \
    LEADING:${params.leading} TRAILING:${params.trailing} \
    SLIDINGWINDOW:${params.slidingwindow} \
    MINLEN:${params.minlen}  \
    2> ${pair_id}_stripped_trimmed_stderr.log 1> ${pair_id}_stripped_trimmed_stdout.log
    mv ${pair_id}_trimmed_1P.fq.gz ${pair_id}_R1_concat_stripped_trimmed.fq.gz
    mv ${pair_id}_trimmed_2P.fq.gz ${pair_id}_R2_concat_stripped_trimmed.fq.gz
    cat ${pair_id}_trimmed_1U.fq.gz ${pair_id}_trimmed_2U.fq.gz > ${pair_id}_S_concat_stripped_trimmed.fq.gz
    """
}

process run_fastqc_trimmed {
    publishDir "${params.out_dir}/fastqc_bbduk_trimmed", mode: "${params.savemode}"
    tag { pair_id }
    label 'one'

    input:
    set pair_id, file(reads) from fq2

    output:
    file "$pair_id" 
    file "${pair_id}" into fastqc_bbduk_trimmed

    """
    mkdir ${pair_id}
    fastqc -q ${reads} -o ${pair_id} -t $task.cpus
    """
}



/*
 * Build assembly with SPAdes
 */
process run_spadesasm {
    publishDir "${params.out_dir}/spades", mode: "${params.savemode}"
    tag { pair_id }
    label 'longtime'

    input:
    set pair_id, file(reads) from reads_trimmed

    output:
    set pair_id, file("${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta") \
    into (assembly_results, tobwa_results)
    file "${pair_id}_spades_scaffolds.fasta"
    file "${pair_id}_spades.log"


    // For 2022 version, params.careful was removed to do --isolate and --only-assembler
    """
    spades.py --cov-cutoff=${params.cov_cutoff} \
    -1 ${pair_id}_R1_concat_stripped_trimmed.fq.gz \
    -2 ${pair_id}_R2_concat_stripped_trimmed.fq.gz \
    -s ${pair_id}_S_concat_stripped_trimmed.fq.gz \
    -t $task.cpus --isolate --only-assembler \
    -o ${pair_id}_spades
    filter_fasta_length.py -i ${pair_id}_spades/scaffolds.fasta \
       -o ${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta \
       -m ${params.min_contig_len}
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
    publishDir "${params.out_dir}/bwamem", mode: "${params.savemode}"
    tag { pair_id }
    label 'longtime'

    input:
    set pair_id, file("${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta"), \
    file(reads) from tobwa_results.join(pilon_reads)

    output:
    set pair_id, file("${pair_id}_mapped_sorted.bam"), \
    file("${pair_id}_mapped_sorted.bam.bai") into bwamem_results
    file "${pair_id}_bwa.log"

    """
    bwa index ${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta
    bwa mem -t $task.cpus  ${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta \
    *.fq.gz 2> ${pair_id}_bwa.log| samtools sort -o ${pair_id}_mapped_sorted.bam -
    samtools index ${pair_id}_mapped_sorted.bam
    """
}

/*
* Incorporating pilon_reads
*/

process run_pilon {
    publishDir "${params.out_dir}/pilon", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file("${pair_id}_mapped_sorted.bam"), \
    file("${pair_id}_mapped_sorted.bam.bai"), \
    file("${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta") \
      from bwamem_results.join(assembly_results)

    output:
    set pair_id, file("${pair_id}_pilon_spades.*") into pilon_results
    set pair_id, file("${pair_id}_pilon_spades.fasta") into to_prokka
    file "${pair_id}_pilon_spades.fasta" into asms_for_quast

    """
    export _JAVA_OPTIONS=$task.javaopts
    pilon --threads $task.cpus --genome ${pair_id}_spades_scaffolds_min${params.min_contig_len}.fasta \
    --bam ${pair_id}_mapped_sorted.bam --output ${pair_id}_pilon_spades \
    --changes --vcfqe &> ${pair_id}_pilon_spades.log
    """
}

/*
* Annotation using PROKKA
*/
process run_prokka {
    publishDir "${params.out_dir}/prokka", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file("${pair_id}_pilon_spades.fasta") from to_prokka

    output:
    set pair_id, file("${pair_id}.*") into annotation_results
    file "${pair_id}.*" into annotation_multiqc

    """
    prokka --compliant --force --usegenus --cpus $task.cpus \
    --centre ${params.centre} --prefix ${pair_id} --locustag ${params.locustag} \
    --genus ${params.genus} --species ${params.species} \
    --kingdom ${params.kingdom} --strain ${pair_id}_prokka_info ${params.prokka_additional} \
    --outdir . ${pair_id}_pilon_spades.fasta
    """
}

/*
 * Evaluate ALL assemblies with QUAST
 */
process run_quast {
    // The output here is a directory in and of itself
    // thus not creating a new one
    publishDir "${params.out_dir}/", mode: "${params.savemode}"
    tag { pair_id }

    input:
    file asm_list from asms_for_quast.toSortedList()

    output:
    file quast_evaluation_all into quast_evaluation_all
    file quast_evaluation_all into quast_multiqc

    """
    quast --threads $task.cpus -o quast_evaluation_all \
    -g ${params.quast_genes} -R ${params.quast_ref} ${asm_list}
    """
}

process run_multiqc_final {
    publishDir "${params.out_dir}/multiqc_final", mode: "${params.savemode}"
    tag {"multiqc"}
    label 'one'

    input:
    file "fastqc_output/*" from fastqc_multiqc.collect()
    file "bbduk/*" from bbduk_stats_stripped_multiqc.collect()
    file "bbduk_trimmed/*" from bbduk_trimmed_multiqc.collect()
    file "bbduk_trimmed_fastqc/*" from fastqc_bbduk_trimmed.collect()
    file "prokka/*" from annotation_multiqc.collect()
    file quast_evaluation_all from quast_multiqc

    output:
    file("*")

    """
    multiqc --fullnames --config ${params.multiqc_config} . 
    """

}
