include { check_samplesheet } from '../modules/check_samplesheet'

include { FILTER } from "../modules/filter.nf"
include { MERGE_DATASETS } from "../modules/merge_datasets.nf"
include { IDENTIFY_HVGS } from "../modules/identify_hvgs.nf"

workflow PREPROCESSING {
    take:
        ch_samplesheet

    main:
        ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))
        FILTER(ch_samples)
        MERGE_DATASETS(FILTER.out.flatMap{ meta, adata -> adata }.collect())

        ch_preprocessed = MERGE_DATASETS.out
            .map{ adata -> [[id: "preprocessed"], adata] }

        IDENTIFY_HVGS(
            ch_preprocessed,
            params.integration_hvgs
        )

    emit:
        simple = ch_preprocessed
        hvgs = IDENTIFY_HVGS.out
}