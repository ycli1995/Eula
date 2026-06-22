#ifndef FAST_STATS
#define FAST_STATS

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

std::vector<double> dgcmatrix_rowsums(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
);

std::vector<double> dgcmatrix_rowmeans(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
);

std::vector<double> dgcmatrix_rowvars_know_means(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    const std::vector<double>& mu
);

std::vector<double> dgcmatrix_rowvars(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
);

std::vector<double> dgcmatrix_rowmeans_expm1(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
);

std::vector<double> dgcmatrix_rowvars_expm1_know_means(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    const std::vector<double>& mu
);

std::vector<double> dgcmatrix_rowvars_expm1(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
);

std::vector<double> matrix_rowmeans_expm1(Eigen::Map<Eigen::MatrixXd> mat);

#endif  // FAST_STATS
