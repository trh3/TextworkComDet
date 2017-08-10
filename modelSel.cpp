//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List modelSel(NumericVector y_s, NumericMatrix x_s, LogicalMatrix design_s){
  arma::mat x(x_s.begin(), x_s.nrow(), x_s.ncol(), false);
  arma::vec y(y_s.begin(), y_s.size(), false);
  arma::imat design_i(design_s.begin(), design_s.nrow(), design_s.ncol(), false);
  arma::umat design = arma::conv_to<arma::umat>::from(design_i);
  arma::vec r2(design_s.ncol());
  arma::vec fs(design_s.ncol());
  arma::vec bf(design_s.ncol());
  arma::vec ypred(y_s.size(), arma::fill::zeros);
  int n = y_s.size();
  arma::vec ycen = y - arma::sum(y)/n;
  for(int i = 0; i < design.n_rows; i++){
    int p = arma::sum(design.row(i));
    arma::mat Xsub = x.cols(design.row(i));
    arma::mat XsubT = Xsub.t();
    arma::mat XX = XsubT * Xsub;
    arma::mat XXi = XX.i();
    arma::mat hat = Xsub * XXi *XsubT;
    
    arma::vec resid = (arma::eye<arma::mat>(n,n) - hat)*y;
    r2(i) = arma::as_scalar(resid.t()*resid)/arma::as_scalar(ycen.t()*ycen);
    
    double f = (r2(i)/p)/((1-r2(i))/(n-1-p)) - 1;

    if(f < 0)
      f = 0;
    
    bf(i) = pow((1+f),(n-p-1)/2) * pow((1+f*(1-r2(i))), -(n-1)/2); 
    
  }
  
  bf = bf/arma::sum(bf);
    for(int i = 0; i < design.n_rows; i++){
      arma::mat Xsub = x.cols(design.row(i));
      arma::mat XsubT = Xsub.t();
      arma::mat XX = XsubT * Xsub;
      arma::mat XXi = XX.i();
      arma::mat hat = Xsub * XXi *XsubT;
    
      ypred = ypred + bf(i)*hat*y;
    }
  
  
 return List::create(Named("y") = ypred);
}

double clamp(double d, double min, double max) {
  const double t = d < min ? min : d;
  return t > max ? max : t;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

