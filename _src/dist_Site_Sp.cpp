#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat dist_Site_Sp(const arma::mat SiteCoord, const arma::mat SpCoord) { 

  // Number of points
  int nsite = SiteCoord.n_rows;
  int nsp = SpCoord.n_rows;
  
  // Initialize with zeros the matrix to store results
  arma::mat distSiteSp; distSiteSp.zeros(nsite, nsp);
  
  // Loop on all points 
  for (int i = 0; i < nsite; i++) {
    arma::vec s = SiteCoord.row(i).t(); // fix a site
    // Loop to calculate the distances between this point and the next ones  
    for (int j = 0; j < nsp; j++) {
      arma::vec sp = SpCoord.row(j).t(); // fix a species
      arma::vec diff = s - sp; // (x0-x1, y0-y1, etc.)
      double squared_diff = as_scalar(diff.t() * diff); // (x0-x1)² + (y0-y1)² + ...
      // Fill the distance matrix with the square root of precedent value
      distSiteSp(i, j) = sqrt(squared_diff);
    }
  }
  return distSiteSp;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// SiteCoord <- matrix(runif(9), ncol=3)
// SpCoord <- matrix(runif(9), ncol=3)
// dist_Site_Sp(SiteCoord, SpCoord)
// */
