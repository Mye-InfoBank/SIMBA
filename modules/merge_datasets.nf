process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  path(adatas)
  
  output:
  path("datasets.integration.h5ad"), emit: integration
  path("datasets.counts.h5ad"), emit: counts
  path("datasets.intersection.h5ad"), emit: intersection
  path("datasets.transfer.h5ad"), emit: transfer, optional: true

  when:
  task.ext.when == null || task.ext.when
  
  script:
  """
  merge_datasets.py --input ${adatas} --output_transfer datasets.transfer.h5ad --output_intersection datasets.intersection.h5ad --output_integration datasets.integration.h5ad --output_counts datasets.counts.h5ad
  """
}