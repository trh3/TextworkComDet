//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
arma::mat dailyActivationSampler(NumericMatrix fullData_s, NumericVector exciteParams_s, NumericVector domainRates, NumericVector control, double priors){
  arma::mat fullData(fullData_s.begin(), fullData_s.nrow(), fullData_s.ncol(),false); 
  arma::mat dailyActivation = arma::zeros<arma::mat>(control[2], control[0]);
  
  arma::vec days = fullData.col(4);

  for(int day = 1; day <= control[0]; day++){

  arma::mat dailySubset = fullData.rows(arma::find(days == day));

	for(int topic = 1; topic <= control[2]; topic++){
		arma::mat subset = dailySubset.rows(arma::find(dailySubset.col(5) == topic));
		double activeLikelihood = 0;
		double inactiveLikelihood = 0;
	
		for(int row = 0; row < subset.n_rows; row++){
  		activeLikelihood += log(R::dpois(subset(row, 0),  (subset(row,1)*domainRates(subset(row, 3)-1) + subset(row,1)*exciteParams_s(topic-1)) ,0));
			inactiveLikelihood += log(R::dpois(subset(row, 0),  (subset(row,1))*domainRates(subset(row, 3)-1) ,0));
		}
		
		activeLikelihood += log(priors);
		inactiveLikelihood += log(1-priors);
		Rcout << "Test1" << std::endl;
		if(subset(1,6) == 1){
			double a = inactiveLikelihood - activeLikelihood;
			
			if(log(Rcpp::runif(1,0,1)(0)) < a){
				dailyActivation(topic-1,day-1) = 0;
			}else{
				dailyActivation(topic-1,day-1) = 1;
			}
		
		
		}else{
			double a = activeLikelihood - inactiveLikelihood;
			
			if(log(Rcpp::runif(1,0,1)(0)) < a){
				dailyActivation(topic-1,day-1) = 1;
			}else{
				dailyActivation(topic-1,day-1) = 0;
			}
		}
	}
	}
		return dailyActivation;
	}
