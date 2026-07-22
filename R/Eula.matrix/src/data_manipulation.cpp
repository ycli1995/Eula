#include "data_manipulation.h"

#include "fast_stats.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd dgcmatrix_row_scale(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    bool scale = true,
    bool center = true,
    double scale_max = 10
) {
    std::vector<double> mu, sigma;
    if (center) {
        mu = dgcmatrix_rowmeans(mat);
        if (scale) {
            sigma = dgcmatrix_rowvars_know_means(mat, mu);
        }
    } else {
        if (scale) {
            sigma = dgcmatrix_rowvars(mat);
        }
    }
    Eigen::MatrixXd scaled_mat =
        dgcmatrix_row_scale2(mat, mu, sigma, scale, center, scale_max);
    return scaled_mat;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd dgcmatrix_row_scale2(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    std::vector<double>& mu,
    std::vector<double>& sigma,
    bool scale = true,
    bool center = true,
    double scale_max = 10
) {
    Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
    for (int i = 0; i < mat.cols(); ++i) {
        scaled_mat.col(i) = mat.col(i);
        for (int j = 0; j < mat.rows(); ++j) {
            if (center) {
                scaled_mat(j, i) -= mu[j];
            }
            if (scale) {
                scaled_mat(j, i) /= sigma[j];
            }
            if (scaled_mat(j, i) > scale_max) {
                scaled_mat(j, i) = scale_max;
            }
        }
    }
    return scaled_mat;
}

// [[Rcpp::export(rng = false)]]
std::vector<double> dgcmatrix_std_rowvar(
    Eigen::Map<Eigen::SparseMatrix<double>> mat,
    std::vector<double>& mu,
    std::vector<double>& sd,
    double clip
) {
    std::vector<double> out(mat.rows(), 0.0);
    std::vector<int> nzero(mat.rows(), mat.cols());

    Eigen::SparseMatrix<double>::StorageIndex i;
    Eigen::SparseMatrix<double>::Scalar v;
    for (int k = 0; k < mat.nonZeros(); ++k) {
        i = mat.innerIndexPtr()[k];
        v = mat.valuePtr()[k];
        if (sd[i] == 0) continue;
        out[i] += std::pow(std::min(clip, (v - mu[i]) / sd[i]), 2);
        nzero[i] -= 1;
    }
    for (int i = 0; i < mat.rows(); ++i) {
        if (sd[i] == 0) continue;
        out[i] += pow((0 - mu[i]) / sd[i], 2) * nzero[i];
        out[i] /= (mat.cols() - 1);
    }
    return out;
}

// [[Rcpp::export(rng = false)]]
std::vector<double> matrix_std_rowvar(
    Eigen::Map<Eigen::MatrixXd> mat,
    std::vector<double>& mu,
    std::vector<double>& sd,
    double clip
) {
    std::vector<double> out(mat.rows(), 0.0);
    for (int i = 0; i < mat.rows(); ++i) {
        if (sd[i] == 0) continue;
        for (int j = 0; j < mat.cols(); ++j) {
            out[i] += std::pow(std::min(clip, (mat(i, j) - mu[i]) / sd[i]), 2);
        }
        out[i] /= (mat.cols() - 1);
    }
    return out;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> dgcmatrix_log_norm(
    Eigen::SparseMatrix<double> data,
    const std::vector<double>& scale_factor
) {
    for (int k = 0; k < data.outerSize(); ++k) {
        double colsum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            colsum += it.value();
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            it.valueRef() =
                std::log1p(double(it.value()) / colsum * scale_factor[k]);
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> dgcmatrix_rc_norm(
    Eigen::SparseMatrix<double> data,
    const std::vector<double>& scale_factor
) {
    for (int k = 0; k < data.outerSize(); ++k) {
        double colsum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            colsum += it.value();
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            it.valueRef() = double(it.value()) / colsum * scale_factor[k];
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> dgcmatrix_clr_norm(
    Eigen::SparseMatrix<double> data
) {
    for (int k = 0; k < data.outerSize(); ++k) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            if (it.value() > 0) {
                sum += std::log1p(it.value());
            }
        }
        sum = std::exp(sum / data.rows());
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            it.valueRef() = std::log1p(double(it.value()) / sum);
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd
matrix_log_norm(Eigen::MatrixXd data, const std::vector<double>& scale_factor) {
    for (int i = 0; i < data.cols(); ++i) {
        double sum = data.col(i).sum();
        for (int j = 0; j < data.rows(); ++j) {
            data(j, i) = std::log1p(data(j, i) / sum * scale_factor[i]);
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd
matrix_rc_norm(Eigen::MatrixXd data, const std::vector<double>& scale_factor) {
    for (int i = 0; i < data.cols(); ++i) {
        double sum = data.col(i).sum();
        for (int j = 0; j < data.rows(); ++j) {
            data(j, i) = data(j, i) / sum * scale_factor[i];
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd matrix_clr_norm(Eigen::MatrixXd data) {
    for (int i = 0; i < data.cols(); ++i) {
        double sum = 0.0;
        for (int j = 0; j < data.rows(); ++j) {
            if (data(j, i) > 0) {
                sum += std::log1p(data(j, i));
            }
        }
        sum = std::exp(sum / data.rows());
        for (int j = 0; j < data.rows(); ++j) {
            data(j, i) = std::log1p(data(j, i) / sum);
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd matrix_column_standardize(Eigen::MatrixXd data) {
    for (int i = 0; i < data.cols(); ++i) {
        double mu = data.col(i).mean();
        double sdev = std::sqrt(
            (data.col(i).array() - mu).square().sum() / (data.rows() - 1)
        );
        for (int j = 0; j < data.rows(); j++) {
            data(j, i) = (data(j, i) - mu) / sdev;
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd matrix_l2norm(Eigen::MatrixXd data, uint32_t margin) {
    if (margin == 1) {
        for (int i = 0; i < data.rows(); ++i) {
            double sum = data.row(i).array().pow(2).sum();
            sum = std::sqrt(sum);
            if (sum == 0) {
                data.row(i).fill(0);
                continue;
            }
            for (int j = 0; j < data.cols(); j++) {
                data(i, j) /= sum;
            }
        }
        return data;
    }
    for (int i = 0; i < data.cols(); ++i) {
        double sum = data.col(i).array().pow(2).sum();
        sum = std::sqrt(sum);
        if (sum == 0) {
            data.col(i).fill(0);
            continue;
        }
        for (int j = 0; j < data.rows(); j++) {
            data(j, i) /= sum;
        }
    }
    return data;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double>
dgcmatrix_l2norm(Eigen::SparseMatrix<double> data, uint32_t margin) {
    if (margin == 1) {
        std::vector<double> rsum(data.rows(), 0.0);
        for (int k = 0; k < data.nonZeros(); ++k) {
            rsum[data.innerIndexPtr()[k]] += std::pow(data.valuePtr()[k], 2);
        }
        for (int k = 0; k < data.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it;
                 ++it) {
                if (rsum[it.index()] == 0) {
                    continue;
                }
                it.valueRef() = it.value() / std::sqrt(rsum[it.index()]);
            }
        }
        return data;
    }
    for (int k = 0; k < data.outerSize(); ++k) {
        double csum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            csum += std::sqrt(it.value());
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it) {
            it.valueRef() = it.value() / std::sqrt(csum);
        }
    }
    return data;
}
