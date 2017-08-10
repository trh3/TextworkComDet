//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List DomainInterestProcess(NumericMatrix domainInterest_s, NumericVector control){
	
  
   arma::mat domainInterest(domainInterest_s.begin(), domainInterest_s.nrow(), domainInterest_s.ncol(), false);
  arma::mat data = arma::zeros<arma::mat>(control[1],control[1]);
  

  
  
  
  
  

    for(int domain = 0; domain < control[1]; domain++){
     for(int domainTwo = 0; domainTwo < domain; domainTwo++){

         data(domain, domainTwo) = dot(domainInterest.col(domain), domainInterest.col(domainTwo));
         data(domainTwo, domain) = dot(domainInterest.col(domain), domainInterest.col(domainTwo));
     } 
    }
  
  return List::create(Named("Links") = data);
  
  
	
}