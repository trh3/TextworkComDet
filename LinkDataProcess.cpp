//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List LinkDataProcess(NumericMatrix linkData_s, NumericMatrix domainInterest_s, NumericVector blockMembership_s, NumericVector lagLink, NumericVector type,  NumericVector control){
	

  arma::cube linkCube = arma::zeros<arma::cube>(control[1], control[1], control[0]);
  NumericMatrix linkData = Rcpp::clone(linkData_s);
  
  arma::mat LinkData(linkData.begin(), linkData.nrow(), linkData.ncol(), false);
  arma::mat domainInterest(domainInterest_s.begin(), domainInterest_s.nrow(), domainInterest_s.ncol(), false);
  arma::mat data = arma::zeros<arma::mat>(control[1]*control[1]*control[0], 4);
  arma::mat auxdata = arma::zeros<arma::mat>(control[1]*control[1]*control[0], 3);
  
  int nDocs = linkData_s.nrow();
  arma::cube lagCube = arma::zeros<arma::cube>(control[1], control[1], control[0]);
  arma::mat lagMat = arma::zeros<arma::mat>( control[1], control[0]);
 
  for(int doc = 0; doc < nDocs; doc++){
    linkCube.slice(LinkData(doc,0)-1).row(LinkData(doc,1)-1) += LinkData(arma::span(doc), arma::span(2, linkData.ncol()-1));
  }
  for(int day = 1; day < control[0]; day++){
    lagCube.slice(day) = lagCube.slice(day-1) + linkCube.slice(day-1);
    if(day-lagLink(0) >= 0){
      lagCube.slice(day) -= linkCube.slice(day-lagLink(0));
      
    }
  lagMat.col(day) = arma::sum(lagCube.slice(day), 1);
  }
  
  Rcout << "Done Cubing" << std::endl;
//  int count = 0;
//  
//  arma::vec lagged = arma::zeros<arma::vec>(control[1]);
//  
//  for(int day = 0; day < control[0]; day++){
//    Rcout << day+1 << std::endl;
//    lagged -= 1;
//    
//    lagged = arma::clamp(lagged, 0, lagLink[0]);
//    
//    for(int domain = 0; domain < control[0]; domain++){
//     for(int domainTwo = 0; domainTwo < control[0]; domainTwo++){
//       if(lagged(domainTwo) > 0){
//         data(count, 3) = 1; 
//       }
//       arma::mat daySlice = linkCube.slice(day);
//       if(daySlice(domain, domainTwo) > 1){
//          data(count,0) = 1;
//       }
//       if(data(count,0) > 0){
//         lagged(domainTwo) = lagLink[0];
//       }
//       data(count,1) = domain +1;
//       data(count,2) = domainTwo+1;
//      
//       auxdata(count,0) = domain +1;
//       auxdata(count,1) = domainTwo +1;
//       if(blockMembership_s(domain) == blockMembership_s(domainTwo)){
//         auxdata(count, 2) = 1;
//         
//       }else{
//         auxdata(count, 2) = dot(domainInterest.col(domain), domainInterest.col(domainTwo));
//         
//       }
//       count++;
//     } 
//    }
//  }
  
  return List::create(Named("Links") = linkCube, Named("Lags") = lagMat);
}