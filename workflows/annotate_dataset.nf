
include { CELLTYPIST }  from "../modules/local/celltypist.nf"
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_FINE }  from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_EPI }  from "../modules/local/jupyternotebook/main.nf"
include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_CELL_TYPES } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as EXPORT_ATLAS }  from "../modules/local/jupyternotebook/main.nf"
include { SCVI } from "../modules/local/scvi/main.nf"
include { SCANVI } from "../modules/local/scvi/main.nf"

/**
 * Annotate cell-types of the lung cancer atlas.
 *   - perform a coarse-grained annotation
 *   - perform DE analysis on some clusters
 *   - manually annotate sub-clusters to obtain a fine-grained cell-type annotation.
 */
workflow annotate_dataset {
    take:
        adata_integrated // id, rep, adata

    main:
    CELLTYPIST(
        adata_integrated,
        true,
        "Cells_Intestinal_Tract.pkl"
    )
    ch_adata_annotated = CELLTYPIST.out
    
    /*
    TODO: Find out if this analysis makes sens
    SPLIT_ANNDATA(
        ch_adata_annotated.map{ it -> [it.baseName, it]},
        "cell_type"
    )
    NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
        SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it] },
        "X_scANVI",
        Channel.from(0.5, 0.75, 1.0, 1.5)
    )
    */
    
    /*
    EXPORT_ATLAS(
        Channel.value([
            [id: "export_atlas"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/35_export_atlas.py")
        ]),
        [
            "adata_annotated_fine": "adata_annotated_fine.h5ad",
            "adata_epi": "adata_epithelial.h5ad",
            "platform_metadata": "sequencing_platforms.csv",
            "patient_metadata": "patient_metadata_corrected.xlsx"
        ],
        ANNOTATE_CELL_TYPES_FINE.out.artifacts.mix(
            ANNOTATE_CELL_TYPES_EPI.out.artifacts
        ).mix(
            Channel.fromPath("$baseDir/tables/additional_patient_metadata/patient_metadata_corrected.xlsx")
        ).mix(
            Channel.fromPath("$baseDir/tables/additional_patient_metadata/sequencing_platforms.csv")
        ).collect()
    )
    ch_atlas = EXPORT_ATLAS.out.artifacts.flatten().filter{ it -> it.baseName.equals("full_atlas_annotated") }

    // re-compute model on final atlas, to be used as reference for scArches
    SCVI(
        ch_atlas.map{ it -> ['full_atlas', it]},
        1,
        ["sample", "dataset", null]
    )
    SCANVI(
        SCVI.out.adata.join(SCVI.out.scvi_model),
        "sample",
        "cell_type"
    )

    emit:
        final_atlas = ch_atlas
        scanvi_model = SCANVI.out.scvi_model.map{ id, model -> model }
        scanvi_h5ad = SCANVI.out.adata.map{ id, adata -> adata }
    */
}
