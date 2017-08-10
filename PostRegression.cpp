//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List postRegressionData(NumericVector assign_s,NumericVector domains_s , NumericVector day_s, 
NumericMatrix domainInterest_s, NumericMatrix activation_s, NumericVector control){
	
  int count = 0;
  arma::vec assign(assign_s.begin(), assign_s.size(), false);
  arma::vec domains(domains_s.begin(), domains_s.size(), false);
  arma::vec days(day_s.begin(), day_s.size(), false);
  arma::mat data = arma::zeros<arma::mat>(control[0]*control[1]*control[2],7); 
  arma::mat domainInterest(domainInterest_s.begin(), control[2], control[1], false);
  arma::mat activation(activation_s.begin(), control[2], control[0], false);
  for(int day = 1; day <= control[0]; day++){
   // Rcout << day << std::endl;
    arma::colvec assignsubset = assign(arma::find(days == day));
   
   
    arma::colvec domainsubset = domains(arma::find(days == day));
    for(int domain = 1; domain <= control[1]; domain++){
              arma::colvec domaintopicsubset = assignsubset(arma::find(domainsubset == domain));
              for(int topic = 1; topic <= control[2]; topic++){
                arma::uvec domaintopictopicsubset = arma::find(domaintopicsubset == topic);
                int numPosts = domaintopictopicsubset.n_elem;
                
//                if(numPosts > 0){
//                  Rcout << domain<< " " << topic <<" "<<numPosts<<std::endl;
//                }
                data(count, 0) = domaintopictopicsubset.n_elem;

                data(count, 1) = domainInterest(topic-1, domain-1);

                data(count, 2) = domainInterest(topic-1, domain-1) * activation_s(topic-1, day-1);
                data(count, 3) = domain;
                data(count, 4) = day;
                data(count, 5) = topic;
                data(count, 6) = activation(topic-1, day-1);
                count++;
                
                
              }
    }

    
    
  }
  
  return List::create(Named("Mat") = data);
  
  
	
}