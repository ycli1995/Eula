
library(Seurat)

mat <- Read10X("/develop/10x_data/rna/pbmc10k/filtered_feature_bc_matrix")

RhpcBLASctl::blas_get_num_procs()
RhpcBLASctl::blas_set_num_threads(1)

bench::mark(
  mat2 <- Eula.matrix:::dgcmatrix_rc_norm_arma(mat, rep.int(10000, ncol(mat)))
)

bench::mark(
  mat3 <- Eula.matrix:::dgcmatrix_rc_norm(mat, rep.int(10000, ncol(mat)))
)

bench::mark(
  mat4 <- Seurat:::RelativeCounts(mat, scale.factor = 10000)
)