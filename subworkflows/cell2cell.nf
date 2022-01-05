include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"


process SQUIDPY {
    publishDir = [ path: { "${params.squidpy_outDir}"}, mode: 'link' ]

    input:
    tuple val(id), val(in_file)

    output:
    path "*.pkl", emit: out_file, optional: true
    path "*.log", emit: log

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    squidpy_cpdb.py -i ${in_file} -o ./ > ${id}.log 2>&1 $ignore_exit_code
    """
}


workflow cell2cell {
    take: adata_annotated

    main:
    ch_adata_annotated = Channel.value([adata_annotated.baseName, adata_annotated])
    SPLIT_ANNDATA(ch_adata_annotated, "sample")
    SQUIDPY(SPLIT_ANNDATA.out.adata.flatten().map{it -> [it.baseName, it]})
}
