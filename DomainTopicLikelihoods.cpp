//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List domainTopicLikelihoods(NumericMatrix domainInterest_s, NumericMatrix activation_s,  NumericVector exciteParams_s, 
NumericVector domainPostRates_s, NumericVector control){
	

  arma::cube likelihoods = arma::zeros<arma::cube>(control[2], control[1], control[0]); 
  arma::mat domainInterest(domainInterest_s.begin(), control[2], control[1], false);
  arma::mat activation(activation_s.begin(), control[2], control[0], false);
  arma::vec exciteParams(exciteParams_s.begin(), control[2], false);
  arma::vec domainPostRates(domainPostRates_s.begin(), control[1], false);
  for(int day = 1; day <= control[0]; day++){
    for(int domain = 1; domain <= control[1]; domain++){
      for(int topic = 1; topic <= control[2]; topic++){ 
        likelihoods(topic-1, domain-1, day-1) = (domainInterest(topic-1, domain-1)*domainPostRates(domain-1) +
                                          activation(topic-1, day-1)*domainInterest(topic-1, domain-1)*exciteParams(topic-1));
      }
      
    }
    

    arma::mat subview = likelihoods.slice(day-1);
    arma::mat nsubview = arma::normalise(subview, 1,0);

    likelihoods.slice(day-1) = nsubview;  
  }
  likelihoods = log(likelihoods);
  return List::create(Named("Mat") = likelihoods);
  
  
	
}