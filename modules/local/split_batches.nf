process SPLIT_BATCHES {
    container "bigdatainbiomedicine/sc-rpy"
    input:
        tuple val(meta), path(input_adata)

    output:
        tuple val(meta), path("*.h5ad"), emit: adata

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad

    adata = ad.read_h5ad("${input_adata}")

    for batch in adata.obs.batch.unique():
        adata[adata.obs.batch == batch].write_h5ad(f"{batch}.h5ad")
    """
}
