//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List BlockAssignment(NumericMatrix fullLinkData_s, NumericVector currentParameters, NumericMatrix domainBlockLikelihoods, NumericMatrix TopicInterestDyads, 
NumericVector blockMemb_s, NumericVector control, int nBlocks){
	
  NumericMatrix fullLinkData_c = Rcpp::clone(fullLinkData_s);
  NumericVector blockMemb_c = Rcpp::clone(blockMemb_s);
  
  arma::mat fullLinkData(fullLinkData_c.begin(), fullLinkData_c.nrow(), fullLinkData_c.ncol());
  arma::vec domainIndex = arma::linspace(1, control[1], control[1]);
  NumericVector domainIndex = NumericVector(topicIndex.begin(), topicIndex.end());
  NumericVector domainOrder = RcppArmadillo::sample(TopicIndex , 1, false);
  arma::mat currParams(currentParameters, currentParameters.begin(),false);
  for(int i = 0; i < control[1]; i++){
    Rcout << i+1 << std::endl;
    NumericVector currBlock = blockMemb_c(domainOrder(i));
    NumericVector propBlock = RcppArmadillo::sample(TopicIndex, 1, false);
    while(propBlock(0) == currBlock(0)){
    NumericVector propBlock = RcppArmadillo::sample(TopicIndex, 1, false);
    }
    
    arma::uvec indices = arma::find((fullLinkData.col(4) == domainOrder(i)) | (fullLinkData.col(5) == domainOrder(i)));
    
    arma::mat dataCurr = fullLinkData.row(indices);
    
    arma::mat dataProp = dataCurr;
    
    for(int j = 0; j < dataProp.nrows(); j ++){
      
      
      
      if(((blockMemb(dataProp(j,4)) == propBlock | (blockMemb(dataProp(j,5)) == propBlock)))){
        dataProp(j,3) = 1;
      }else{
        if(dataProp(j, 4) == domainOrder(i)){
          dataProp(j,3) = TopicInterestDyads(domainOrder(i), dataProp(j, 5)); 
        }else{
          dataProp(j,3) = TopicInterestDyads(domainOrder(i), dataProp(j, 4));  
        } 
      }  
    }
    
    double currLik = 0;
    double propLik = 0;
    
    for(int j = 0; j < dataProp.nrows(); j ++){
      double currLikNum = arma::exp(dataCurr(j, 1)*currParams(0) +dataCurr(j, 2)*currParams(1) + dataCurr(j, 3)*currParams(2));
      double currLikDem = 1 + currLikNum;
      
      currLikNum = arma::log(currLikNum)*dataCurr(j,0)*dataCurr(j, 6);
      currLikDem = arma::log(currLikdem)*(dataCurr(j, 6) + dataCurr(j,7));
    
      currLik += currLikNum - currLikDem;

      double propLikNum = arma::exp(dataProp(j, 1)*currParams(0) +dataProp(j, 2)*currParams(1) + dataProp(j, 3)*currParams(2));
      double propLikDem = 1 + propLikNum;
      
      propLikNum = arma::log(propLikNum)*dataProp(j,0)*dataProp(j, 6);
      propLikDem = arma::log(propLikdem)*(dataProp(j, 6) + dataProp(j,7));
    
      propLik += propLikNum - propLikDem;
    }
    
    double a = (propLik + domainBlockLikelihoods(domainOrder(i), propBlock)) - (currLik + domainBlockLikelihoods(domainOrder(i), currBlock));
    
    if(arma::log(Rcpp::runif(1,0,1)) < a){
      blockMemb(domainOrder(i)) = propBlock;
      fullLinkData.row(indices) = dataProp;
      Rcout << "Change!" << std::endl;
    }
    
    
  }


  return List::create(Named("updatedData") = fullLinkData, Named("updatedBlockMemb") = blockMemb);
}