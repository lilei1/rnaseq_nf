/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'

//
// MODULE: Local modules
//
include { INDEX_GENOME           } from '../modules/local/index_genome'
include { SPLIT_REFERENCE_5P3P   } from '../modules/local/split_reference'
include { LIB_PROCESS            } from '../modules/local/lib_process'
include { GATHER_COUNTS          } from '../modules/local/gather_counts'
include { DGE_ANALYSIS           } from '../modules/local/dge_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Index genome if provided
    //
    ch_genome_index = Channel.empty()
    if (params.genome_fasta) {
        INDEX_GENOME (
            Channel.fromPath(params.genome_fasta)
        )
        ch_genome_index = INDEX_GENOME.out.index
        ch_versions = ch_versions.mix(INDEX_GENOME.out.versions)
    }

    //
    // MODULE: Split reference sequences
    //
    ch_genome_annotation = params.genome_annotation ? Channel.fromPath(params.genome_annotation) : Channel.empty()
    ch_transcriptome_fasta = params.transcriptome_fasta ? Channel.fromPath(params.transcriptome_fasta) : Channel.empty()

    SPLIT_REFERENCE_5P3P (
        ch_genome_annotation,
        ch_transcriptome_fasta,
        params.ht_type,
        params.ht_id,
        params.map_transcriptome
    )
    ch_versions = ch_versions.mix(SPLIT_REFERENCE_5P3P.out.versions)

    //
    // MODULE: Process each library
    //
    LIB_PROCESS (
        ch_samplesheet,
        ch_genome_index,
        ch_genome_annotation,
        ch_transcriptome_fasta,
        SPLIT_REFERENCE_5P3P.out.cleaned_transcriptome_fasta,
        SPLIT_REFERENCE_5P3P.out.split_annotation,
        SPLIT_REFERENCE_5P3P.out.transcripts_split_fasta,
        params.ht_type,
        params.ht_id,
        params.map_transcriptome,
        params.assign_primary
    )
    ch_versions = ch_versions.mix(LIB_PROCESS.out.versions.first())

    //
    // MODULE: Gather counts from all libraries
    //
    GATHER_COUNTS (
        LIB_PROCESS.out.feature_counts.collect(),
        LIB_PROCESS.out.feature_counts_fwd.collect(),
        LIB_PROCESS.out.feature_counts_summary.collect(),
        LIB_PROCESS.out.feature_counts_unstranded.collect(),
        LIB_PROCESS.out.feature_counts_unstranded_summary.collect(),
        LIB_PROCESS.out.feature_counts_split.collect(),
        LIB_PROCESS.out.mt_counts.collect(),
        LIB_PROCESS.out.mt_counts_fwd.collect(),
        LIB_PROCESS.out.mt_counts_split.collect(),
        LIB_PROCESS.out.mt_statsfile.collect(),
        params.input,
        params.partial_counts_results ? Channel.fromPath(params.partial_counts_results) : Channel.empty(),
        params.map_transcriptome,
        params.stranded_reads,
        params.geneid_string
    )
    ch_versions = ch_versions.mix(GATHER_COUNTS.out.versions)

    //
    // MODULE: Differential gene expression analysis
    //
    if (!params.skip_dge) {
        DGE_ANALYSIS (
            GATHER_COUNTS.out.lib_counts_tarball,
            Channel.fromPath(params.input),
            ch_genome_annotation,
            params.ht_type,
            params.ht_id,
            params.kegg_ids ? Channel.fromPath(params.kegg_ids) : Channel.empty(),
            params.kegg_pathway_ko_list ? Channel.fromPath(params.kegg_pathway_ko_list) : Channel.empty(),
            params.kegg_pathway_metadata ? Channel.fromPath(params.kegg_pathway_metadata) : Channel.empty(),
            params.ko_terms ? Channel.fromPath(params.ko_terms) : Channel.empty(),
            params.map_transcriptome,
            params.group_counts_by,
            params.stranded_reads,
            params.skip_dge,
            params.skip_rep_analysis,
            params.geneid_string
        )
        ch_versions = ch_versions.mix(DGE_ANALYSIS.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', value: it))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/