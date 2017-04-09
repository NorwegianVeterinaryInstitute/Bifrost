
println params.setsize
println params.file_ending
Channel
    .fromFilePairs( params.reads, size:params.setsize )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ in_read_pairs }

process collate {
    input:
    set pair_id, file(reads) from in_read_pairs

    output:
    set pair_id, file("${pair_id}_R{1,2}${params.file_ending}") into read_pairs_out

    """
    cat ${pair_id}*R1*${params.file_ending} > ${pair_id}_R1${params.file_ending}
    cat ${pair_id}*R2*${params.file_ending} > ${pair_id}_R2${params.file_ending}
    """
}


process test {
    input:
    set one, two, three from read_pairs_out

    """
    echo ${one}, ${two[0]}, ${two[1]}
    """

}