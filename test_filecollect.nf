// First, define the input data that go into input channels
//params.setsize=4
//params.reads = "../testdata/short/*L001_R{1,2}_001.short.fastq.gz"

println params.setsize
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ in_read_pairs }



process collate {
    input:
    set pair_id, file(reads) from in_read_pairs

    output:
    set pair_id, file([${pair_id}_R1.fastq.gz, ${pair_id}_R2.fastq.gz]) into read_pairs

    """
    echo "file set is", ${params.setsize}
    echo cat R1* ${pair_id}_R1.fastq.gz
    echo cat R2* ${pair_id}_R2.fastq.gz
    """
}