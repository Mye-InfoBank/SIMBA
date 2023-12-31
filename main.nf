#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

// Modules
include { CELLTYPIST } from "./modules/celltypist.nf"
include { SOLO } from "./modules/solo.nf"
include { CELL_CYCLE } from "./modules/cell_cycle.nf"
include { MERGE } from "./modules/merge.nf"

// Workflows
include { PREPROCESSING } from "./workflows/preprocessing.nf"
include { COUNTS } from "./workflows/counts.nf"
include { INTEGRATION } from "./workflows/integration.nf"
include { CLUSTERING } from "./workflows/clustering.nf"
include { BENCHMARKING } from "./workflows/benchmarking.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }
if (!params.celltypist_model) { exit 1, 'CellTypist model not specified!' }

workflow {
    PREPROCESSING(ch_samplesheet)

    ch_preprocessed = PREPROCESSING.out

    COUNTS(ch_preprocessed, params.normalization_method)

    CELLTYPIST(
        ch_preprocessed,
        params.celltypist_model
    )

    INTEGRATION(
        ch_preprocessed,
        Channel.from(params.integration_methods)
    )

    if (params.benchmark) {
        BENCHMARKING(
            ch_preprocessed,
            INTEGRATION.out.integrated_types,
            params.benchmark_hvgs
        )
    }

    SOLO(
        ch_preprocessed,
        INTEGRATION.out.scanvi_model
    )

    CELL_CYCLE(
        ch_preprocessed,
        "human"
    )

    ch_resolutions = Channel.from(params.clustering_resolutions)

    CLUSTERING(
        INTEGRATION.out.integrated,
        ch_resolutions,
        CELLTYPIST.out
    )

    MERGE(
        ch_preprocessed,
        CLUSTERING.out.map { meta, adata -> meta.integration }.collect(),
        CLUSTERING.out.map { meta, adata -> adata }.collect(),
        SOLO.out,
        COUNTS.out,
        CELLTYPIST.out,
        CELL_CYCLE.out,
        ch_resolutions.collect()
    )

}
