#include "fast_stats.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowsums(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
) {
    std::vector<double> out(mat.rows(), 0.0);
    for (int k = 0; k < mat.nonZeros(); ++k) {
        out[mat.innerIndexPtr()[k]] += mat.valuePtr()[k];
    }
    return out;
}

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowmeans(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
) {
    std::vector<double> out = dgcmatrix_rowsums(mat);
    for (std::vector<double>::iterator it = out.begin(); it != out.end();
         ++it) {
        *it /= mat.cols();
    }
    return out;
}

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowvars_know_means(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    const std::vector<double>& mu
) {
    std::vector<double> out(mat.rows(), 0.0);
    std::vector<int> nzero(mat.rows(), mat.cols());

    for (int k = 0; k < mat.nonZeros(); ++k) {
        out[mat.innerIndexPtr()[k]] +=
            std::pow(mat.valuePtr()[k] - mu[mat.innerIndexPtr()[k]], 2);
        nzero[mat.innerIndexPtr()[k]] -= 1;
    }
    for (int k = 0; k < mat.rows(); ++k) {
        out[k] = (out[k] + (std::pow(mu[k], 2) * nzero[k])) / (mat.cols() - 1);
    }
    return out;
}

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowvars(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
) {
    std::vector<double> mu = dgcmatrix_rowmeans(mat);
    std::vector<double> out = dgcmatrix_rowvars_know_means(mat, mu);
    return out;
}

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowmeans_expm1(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
) {
    std::vector<double> out(mat.rows(), 0.0);
    for (int k = 0; k < mat.nonZeros(); ++k) {
        out[mat.innerIndexPtr()[k]] += std::expm1(mat.valuePtr()[k]);
    }
    for (std::vector<double>::iterator it = out.begin(); it != out.end();
         ++it) {
        *it /= mat.cols();
    }
    return out;
}

// [[Rcpp::export]]
std::vector<double> matrix_rowmeans_expm1(Eigen::Map<Eigen::MatrixXd> mat) {
    std::vector<double> out(mat.rows(), 0.0);
    for (int i = 0; i < mat.cols(); ++i) {
        for (int j = 0; j < mat.rows(); ++j) {
            out[j] += std::expm1(mat(j, i));
        }
    }
    for (std::vector<double>::iterator it = out.begin(); it != out.end();
         ++it) {
        *it /= mat.cols();
    }
    return out;
}

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowvars_expm1_know_means(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    const std::vector<double>& mu
) {
    std::vector<double> out(mat.rows(), 0.0);
    std::vector<int> nzero(mat.rows(), mat.cols());

    for (int k = 0; k < mat.nonZeros(); ++k) {
        out[mat.innerIndexPtr()[k]] += std::pow(
            std::expm1(mat.valuePtr()[k]) - mu[mat.innerIndexPtr()[k]],
            2
        );
        nzero[mat.innerIndexPtr()[k]] -= 1;
    }
    for (int k = 0; k < mat.rows(); ++k) {
        out[k] = (out[k] + (std::pow(mu[k], 2) * nzero[k])) / (mat.cols() - 1);
    }
    return out;
}

// [[Rcpp::export]]
std::vector<double> dgcmatrix_rowvars_expm1(
    Eigen::Map<Eigen::SparseMatrix<double>> mat
) {
    std::vector<double> mu = dgcmatrix_rowmeans_expm1(mat);
    std::vector<double> out = dgcmatrix_rowvars_expm1_know_means(mat, mu);
    return out;
}
