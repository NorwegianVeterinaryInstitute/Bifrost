params.reads = "short/*R{1,2}_001.short.fastq.gz"
println params.reads
params.outputdir = "fastqc_results"

Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }

process run_fastqc {
    input:
    set pair_id, file(reads) from read_pairs
    
    output:
    file "$pair_id" into fastqc_results
    
    """
    mkdir ${pair_id}
    fastqc -q ${reads} -o ${pair_id}
    """
}

process run_fastqc_qc {
    publishDir "outputresults"


    input:
    file "testout/*" from fastqc_results.toSortedList()
    
    output:
    file "woow"
    
    """
    python3 $HOME/PycharmProjects/Bifrost/bifrost_py/fastqc_eval.py -d testout -o woow
    
    """
}



