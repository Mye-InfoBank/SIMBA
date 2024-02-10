include { check_samplesheet } from '../modules/check_samplesheet'

include { RDS_TO_H5AD } from "../modules/rds_to_h5ad.nf"
include { FILTER } from "../modules/filter.nf"
include { GENES_UPSET } from "../modules/genes_upset.nf"
include { MERGE_DATASETS } from "../modules/merge_datasets.nf"
include { COMPOSITION } from "../modules/composition.nf"
include { DISTRIBUTION } from "../modules/distribution.nf"
include { IDENTIFY_HVGS } from "../modules/identify_hvgs.nf"

workflow PREPROCESSING {
    take:
        ch_samplesheet

    main:
        ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))
            .branch { 
                h5ad: it[2] == "h5ad"
                rds: it[2] == "rds"
            }
        
        RDS_TO_H5AD(ch_samples.rds.map{ meta, rds, format -> [meta, rds]})

        FILTER(ch_samples.h5ad.map{ meta, adata, format -> [meta, adata]}.mix(RDS_TO_H5AD.out))
        GENES_UPSET(FILTER.out.map{ meta, adata -> adata }.collect())
        MERGE_DATASETS(FILTER.out.flatMap{ meta, adata -> adata }.collect())

        ch_adata_integration = MERGE_DATASETS.out.integration
            .map{ adata -> [[id: "integration"], adata] }
        
        ch_adata_intersection = MERGE_DATASETS.out.intersection
            .map{ adata -> [[id: "intersection"], adata] }

        ch_adata_counts = MERGE_DATASETS.out.counts
            .map{ adata -> [[id: "counts"], adata] }

        ch_transfer = MERGE_DATASETS.out.transfer.flatten()
            .map{ adata -> [[id: adata.simpleName], adata]}

        COMPOSITION(ch_adata_intersection)
        DISTRIBUTION(ch_adata_intersection)

        IDENTIFY_HVGS(
            ch_adata_integration,
            params.integration_hvgs
        )

    emit:
        integration = ch_adata_integration
        intersection = ch_adata_intersection
        counts = ch_adata_counts
        transfer = ch_transfer
        hvgs = IDENTIFY_HVGS.out
}