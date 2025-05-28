process LIB_PROCESS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hisat2=2.2.1 bioconda::samtools=1.17 bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/your-container' :
        'your_container_for_lib_process' }"

    input:
    tuple val(meta), path(reads)
    path genome_index
    path genome_annotation
    path transcriptome_fasta
    path transcriptome_fasta_cleaned
    path split_annotation
    path transcripts_split_fasta
    val ht_type
    val ht_id
    val map_transcriptome
    val assign_primary

    output:
    tuple val(meta), path("*.feature_counts.txt")                    , optional: true, emit: feature_counts
    tuple val(meta), path("*.feature_counts_fwd.txt")               , optional: true, emit: feature_counts_fwd
    tuple val(meta), path("*.feature_counts_summary.txt")           , optional: true, emit: feature_counts_summary
    tuple val(meta), path("*.feature_counts_unstranded.txt")        , optional: true, emit: feature_counts_unstranded
    tuple val(meta), path("*.feature_counts_unstranded_summary.txt"), optional: true, emit: feature_counts_unstranded_summary
    tuple val(meta), path("*.feature_counts_split.txt")             , optional: true, emit: feature_counts_split
    tuple val(meta), path("*.mt_counts.txt")                        , optional: true, emit: mt_counts
    tuple val(meta), path("*.mt_counts_fwd.txt")                    , optional: true, emit: mt_counts_fwd
    tuple val(meta), path("*.mt_counts_split.txt")                  , optional: true, emit: mt_counts_split
    tuple val(meta), path("*.mt_statsfile.txt")                     , optional: true, emit: mt_statsfile
    tuple val(meta), path("*.bam")                                  , optional: true, emit: bam
    tuple val(meta), path("*.bam.bai")                              , optional: true, emit: bai
    tuple val(meta), path("*.flagstats.txt")                        , optional: true, emit: flagstats
    tuple val(meta), path("*.bw")                                   , optional: true, emit: bigwig
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = meta.strandedness ?: 'unstranded'
    """
    # Extract genome index
    tar -xzf $genome_index

    # Align reads with HISAT2
    hisat2 \\
        -x genome_index_files/genome_ref \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --threads $task.cpus \\
        $args \\
        | samtools sort -@ $task.cpus -o ${prefix}.bam -

    # Index BAM file
    samtools index ${prefix}.bam

    # Generate flagstats
    samtools flagstat ${prefix}.bam > ${prefix}.flagstats.txt

    # Count features with featureCounts
    if [ -f "$genome_annotation" ]; then
        featureCounts \\
            -a $genome_annotation \\
            -o ${prefix}.feature_counts.txt \\
            -t $ht_type \\
            -g $ht_id \\
            -T $task.cpus \\
            ${prefix}.bam

        # Stranded counts
        featureCounts \\
            -a $genome_annotation \\
            -o ${prefix}.feature_counts_fwd.txt \\
            -t $ht_type \\
            -g $ht_id \\
            -s 1 \\
            -T $task.cpus \\
            ${prefix}.bam

        # Unstranded counts
        featureCounts \\
            -a $genome_annotation \\
            -o ${prefix}.feature_counts_unstranded.txt \\
            -t $ht_type \\
            -g $ht_id \\
            -s 0 \\
            -T $task.cpus \\
            ${prefix}.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(echo \$(hisat2 --version 2>&1) | sed 's/^.*hisat2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        subread: \$(echo \$(featureCounts -v 2>&1) | sed 's/^.*featureCounts v//; s/ .*\$//')
    END_VERSIONS
    """
}