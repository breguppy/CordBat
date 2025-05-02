// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//------------------------------------------------------------------------------
// Helper: sum but ignore NaNs (like na.rm = TRUE)
//------------------------------------------------------------------------------
static double safe_sum(const arma::vec& v) {
  double s = 0.0;
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    double x = v[i];
    if (!std::isnan(x)) s += x;
  }
  return s;
}

//------------------------------------------------------------------------------
// Helper: dot but ignore NaNs
//------------------------------------------------------------------------------
static double safe_dot(const arma::vec& x, const arma::vec& y) {
  double s = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double xi = x[i], yi = y[i];
    if (!std::isnan(xi) && !std::isnan(yi)) s += xi * yi;
  }
  return s;
}


//' @name updateCorrectCoefCpp
//' @rdname updateCorrectCoefCpp
//' @title   High‐performance C++ core for update.CorrectCoef
//' @description
//'   RcppArmadillo implementation of the “a”/“b” coefficient updater used by CordBat.
//' @param X0_glist   List of G reference‐batch matrices (each n0×p).
//' @param X1_glist   List of G non‐reference batch matrices (each n1×p).
//' @param ThetaList  List of G precision (p×p) matrices.
//' @param a_i        Numeric vector length p of current “a” coefficients.
//' @param b_i        Numeric vector length p of current “b” coefficients.
//' @param penal_ksi    Numeric penalty on “a” (soft‐threshold λ for a‐update).
//' @param penal_gamma  Numeric penalty on “b” (soft‐threshold λ for b‐update).
//' @return A list with components \code{coef.a} and \code{coef.b} (each numeric length p).
//' @export
// [[Rcpp::export]]
List updateCorrectCoefCpp(
    List X0_glist,
    List X1_glist,
    List ThetaList,
    NumericVector a_i,
    NumericVector b_i,
    double penal_ksi,
    double penal_gamma
) {
  int G = X0_glist.size();
  int p = a_i.size();
  
  // Copy inputs into Armadillo objects
  arma::vec a     = as<arma::vec>(a_i);
  arma::vec b     = as<arma::vec>(b_i);
  arma::vec new_a = a;
  arma::vec new_b = b;
  
  // Pre‐load all matrices
  std::vector<arma::mat> X0(G), X1(G), Theta(G);
  for (int g = 0; g < G; ++g) {
    X0[g]    = as<arma::mat>(X0_glist[g]);
    X1[g]    = as<arma::mat>(X1_glist[g]);
    Theta[g] = as<arma::mat>(ThetaList[g]);
  }
  
  // Loop over each feature j
  for (int j = 0; j < p; ++j) {
    // === 1) a‐update ===
    arma::vec tmp1(G), tmp2(G);
    for (int g = 0; g < G; ++g) {
      int n0 = X0[g].n_rows;
      int n1 = X1[g].n_rows;
      int Ng = n0 + n1;
      
      // Corrected X1 with current a, b
      arma::mat A   = arma::diagmat(a);
      arma::mat B   = arma::repmat(b.t(), n1, 1);
      arma::mat X1c = X1[g] * A + B;
      
      // Stack and compute mean/sd
      arma::mat Xall = arma::join_cols(X0[g], X1c);
      arma::rowvec Mu  = arma::mean(Xall, 0);
      arma::rowvec Sig = arma::stddev(Xall, 0, 0);
      
      // Warn & fix zero/NA
      for (int k = 0; k < p; ++k) {
        if (std::isnan(Sig[k]) || std::abs(Sig[k]) < 1e-6) {
          Sig[k] = 1e-6;
          Rcpp::warning("Sigma_g[%d] was zero or NA for group %d; replaced with 1e-6",
                        k+1, g+1);
        }
      }
      
      // Scale
      arma::mat Xsca = Xall;
      Xsca.each_row() -= Mu;
      Xsca.each_row() /= Sig;
      
      // Extract non‐ref rows and override column j
      arma::mat Yg = Xsca.rows(n0, Ng-1);
      Yg.col(j).fill((b[j] - Mu[j]) / Sig[j]);
      
      // Compute gradients with NaN‐safe dot
      tmp1[g] = 2.0/(Ng * Sig[j]) * safe_dot(X1[g].col(j), Yg * Theta[g].col(j));
      tmp2[g] = 2.0 * Theta[g](j,j)/(Ng * Sig[j]*Sig[j])
        * safe_dot(X1[g].col(j), X1[g].col(j));
    }
    
    // Soft‐threshold update for a[j]
    double denomA = safe_sum(tmp2);
    if (std::abs(denomA) > 1e-6) {
      double z = - safe_sum(tmp1) - safe_sum(tmp2);
      double s = std::copysign(std::max(std::abs(z) - penal_ksi, 0.0), z);
      new_a[j]  = 1.0 + s / denomA;
    }
    a[j] = new_a[j];  // update for next step
    
    // === 2) b‐update ===
    arma::vec tmp3(G), tmp4(G);
    for (int g = 0; g < G; ++g) {
      int n0 = X0[g].n_rows;
      int n1 = X1[g].n_rows;
      int Ng = n0 + n1;
      
      // Re‐correct X1 with updated a
      arma::mat A2   = arma::diagmat(a);
      arma::mat B2   = arma::repmat(b.t(), n1, 1);
      arma::mat X1c2 = X1[g] * A2 + B2;
      
      arma::mat Xall2 = arma::join_cols(X0[g], X1c2);
      arma::rowvec Mu2  = arma::mean(Xall2, 0);
      arma::rowvec Sig2 = arma::stddev(Xall2, 0, 0);
      
      for (int k = 0; k < p; ++k) {
        if (std::isnan(Sig2[k]) || std::abs(Sig2[k]) < 1e-6) {
          Sig2[k] = 1e-6;
          Rcpp::warning("Sigma_g[%d] was zero or NA for group %d; replaced with 1e-6",
                        k+1, g+1);
        }
      }
      
      arma::mat Xsca2 = Xall2;
      Xsca2.each_row() -= Mu2;
      Xsca2.each_row() /= Sig2;
      
      arma::mat Zg    = Xsca2.rows(n0, Ng-1);
      Zg.col(j)       = (X1[g].col(j) * a[j] - Mu2[j]) / Sig2[j];
      
      tmp3[g] = 2.0/(Ng * Sig2[j]) * safe_sum(Zg * Theta[g].col(j));
      tmp4[g] = 2.0 * n1 * Theta[g](j,j)/(Ng * Sig2[j]*Sig2[j]);
    }
    
    // Soft‐threshold update for b[j]
    double denomB = safe_sum(tmp4);
    if (std::abs(denomB) > 1e-6) {
      double z2 = - safe_sum(tmp3);
      double s2 = std::copysign(std::max(std::abs(z2) - penal_gamma, 0.0), z2);
      new_b[j]   = s2 / denomB;
    }
  }
  
  // Wrap and drop dim to get plain vectors
  NumericVector coefA = wrap(new_a);
  NumericVector coefB = wrap(new_b);
  coefA.attr("dim") = R_NilValue;
  coefB.attr("dim") = R_NilValue;
  
  return List::create(
    _["coef.a"] = coefA,
    _["coef.b"] = coefB
  );
}