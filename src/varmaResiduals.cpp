#include <Rcpp.h>
#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]
//' Calculate the residuals for a VARMA model
//'
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd varmaResiduals(const Eigen::Map<Eigen::MatrixXd> &zt,
                               const Eigen::Map<Eigen::VectorXd> &Ph0,
                               const Eigen::Map<Eigen::MatrixXd> &PH,
                               const Eigen::Map<Eigen::MatrixXd> &TH,
                               int p,
                               int q,
                               int include_mean,
                               const Eigen::Map<Eigen::MatrixXd> &beta) {
    int nT = zt.rows();
    int k = zt.cols();
    int pqmax = std::max(p,q);
    Eigen::MatrixXd at = zt;
    at.row(0) -= Ph0.transpose();
    if (pqmax > 1) {
        for (int t = 1; t < pqmax; t++) {
            Eigen::VectorXd tmp = zt.row(t) - Ph0.transpose();
            if (p > 0) {
                for (int j = 0; j < p; j++) {
                    if ((t - j) > 0) {
                        int jdx = j * k;
                        Eigen::VectorXd tmp1 = zt.row(t-j-1) * PH.middleRows(jdx, k);
                        tmp = tmp - tmp1;
                    }
                }
            }
            if (q > 0) {
                for (int j = 0; j < q; j++) {
                    int jdx = j * k;
                    if ((t - j) > 0) {
                        Eigen::VectorXd tmp2 = at.row(t-j-1) * TH.middleRows(jdx, k);
                        tmp = tmp - tmp2;
                    }
                }
            }
            at.row(t) = tmp.transpose();
        }
    }

    for (int t = pqmax; t < nT; t++) {
        Eigen::VectorXd Past((p+q)*zt.cols() + include_mean);
        Past(0) = 1;
        if (p > 0) {
            for (int j = 0; j < p; j++) {
                Past.segment(include_mean+j*zt.cols(), zt.cols()) = zt.row(t - j - 1);
            }
        }
        if (q > 0) {
            for (int j = 0; j < q; j++) {
                Past.segment(include_mean+j+p*zt.cols(), at.cols()) = at.row(t - j - 1);
            }
        }
        Eigen::VectorXd tmp = Past.transpose() * beta;
        Eigen::VectorXd tmp3 = zt.row(t) - tmp.transpose();
        at.row(t) = tmp3.transpose();
    }
    return at.bottomRows(nT-pqmax);
}
