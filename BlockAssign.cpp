//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// List BlockAssignment(NumericMatrix fullLinkData_s, NumericVector currentParameters, NumericMatrix domainBlockLikelihoods, NumericMatrix TopicInterestDyads, 
// NumericVector blockMemb_s, NumericVector control, int nBlocks){
// 	
//   NumericMatrix fullLinkData_c = Rcpp::clone(fullLinkData_s);
//   NumericVector blockMemb_c = Rcpp::clone(blockMemb_s);
//   
//   arma::vec blockMemb(blockMemb_c.begin(), control[1], false);
//   arma::mat fullLinkData(fullLinkData_c.begin(), fullLinkData_c.nrow(), fullLinkData_c.ncol());
//   arma::vec topicIndex = arma::linspace(1, control[1], control[1]);
//   NumericVector domainIndex = NumericVector(topicIndex.begin(), topicIndex.end());
//   NumericVector domainOrder = RcppArmadillo::sample(domainIndex , control[1], true);
//   arma::vec blockIndex = arma::linspace(1,nBlocks, nBlocks);
//   NumericVector blockINDEX = NumericVector(blockIndex.begin(), blockIndex.end());
//   arma::vec currParams(currentParameters.begin(), 3,false);
//   for(int i = 0; i < control[1]; i++){
//     NumericVector currBlock = blockMemb_c(domainOrder(i)-1);
// 
//     arma::vec CurrBlock(currBlock.begin(),1,false);
//     NumericVector propBlock = RcppArmadillo::sample(blockINDEX, 1, false);
//     arma::vec PropBlock(propBlock.begin(),1,false);
//     
//     while(propBlock(0) == currBlock(0)){
//     NumericVector propBlock = RcppArmadillo::sample(blockINDEX, 1, false);
//     arma::vec PropBlock(propBlock.begin(),1,false);
//     }
// 
//     arma::uvec indexProto = (fullLinkData.col(4) == domainOrder(i)) + (fullLinkData.col(5) == domainOrder(i));
//     
//     arma::uvec indices = arma::find(indexProto);
//     
//     arma::mat dataCurr = fullLinkData.rows(indices);
//     
//     arma::mat dataProp = dataCurr;
//     
//     for(int j = 0; j < dataProp.n_rows; j ++){
//         if(((blockMemb(dataProp(j,4)-1) == PropBlock(0)) + (blockMemb(dataProp(j,5)-1) == PropBlock(0))) > 0){
//         dataProp(j,3) = 1;
//       }else{
//         if(dataProp(j, 4) == domainOrder(i)){
//           dataProp(j,3) = TopicInterestDyads(domainOrder(i)-1, dataProp(j, 5)-1); 
//         }else{
//           dataProp(j,3) = TopicInterestDyads(domainOrder(i)-1, dataProp(j, 4)-1);  
//         } 
//       }  
//     }
//     
//     double currLik = 0;
//     double propLik = 0;
//     
//     for(int j = 0; j < dataProp.n_rows; j ++){
//       double currLikNum = exp(dataCurr(j, 1)*currParams(0) +dataCurr(j, 2)*currParams(1) + dataCurr(j, 3)*currParams(2));
//       double currLikDem = 1 + currLikNum;
//       
//       currLikNum = log(currLikNum)*dataCurr(j,0)*dataCurr(j, 6);
//       currLikDem = log(currLikDem)*(dataCurr(j, 6));
//     
//       currLik += currLikNum - currLikDem;
// 
//       double propLikNum = exp(dataProp(j, 1)*currParams(0) +dataProp(j, 2)*currParams(1) + dataProp(j, 3)*currParams(2));
//       double propLikDem = 1 + propLikNum;
//       
//       propLikNum = log(propLikNum)*dataProp(j,0)*dataProp(j, 6);
//       propLikDem = log(propLikDem)*(dataProp(j, 6));
//     
//       propLik += propLikNum - propLikDem;
//     }
//     
//     double a = (propLik + domainBlockLikelihoods(domainOrder(i)-1, propBlock(0)-1)) - (currLik + domainBlockLikelihoods(domainOrder(i)-1, currBlock(0)-1));
//     
//     if(log(Rcpp::runif(1,0,1)(0)) < a){
//       blockMemb(domainOrder(i)-1) = PropBlock(0);
//       fullLinkData.rows(indices) = dataProp;
//     }
//     
//     
//   }
// 
// 
//   return List::create(Named("updatedData") = fullLinkData, Named("updatedBlockMemb") = blockMemb);
// }

//[[Rcpp::export]]
List BlockAssignment(NumericMatrix fullLinkData_s, NumericVector currentParameters, NumericVector domainBlockLikelihood, NumericVector consideredBlocks , NumericMatrix TopicInterestDyads, 
                     NumericVector blockMemb_s, int domain){
   
  NumericMatrix fullLinkData_c = Rcpp::clone(fullLinkData_s);
  NumericVector blockMemb_c = Rcpp::clone(blockMemb_s);
  NumericVector domainBlockLikelihood_c = Rcpp::clone(domainBlockLikelihood);
  
  arma::vec blockMemb(blockMemb_c.begin(), blockMemb_c.size(), false);
  arma::mat fullLinkData(fullLinkData_c.begin(), fullLinkData_c.nrow(), fullLinkData_c.ncol());
  arma::vec topicIndex = arma::linspace(1, consideredBlocks.size() ,consideredBlocks.size());
  NumericVector domainIndex = NumericVector(topicIndex.begin(), topicIndex.end());
  arma::vec blockLik(domainBlockLikelihood_c.begin(),domainBlockLikelihood_c.size() ,false);
  arma::vec currParams(currentParameters.begin(), 5,false);
  
  for(int i = 0; i < consideredBlocks.size(); i++){
    //Rcout << "Test 1" << std::endl;
    int PropBlock = consideredBlocks(i);
    
    arma::uvec indexProto = (fullLinkData.col(4) == domain) + (fullLinkData.col(5) == domain);
    
    arma::uvec indices = arma::find(indexProto);
    
    arma::mat dataCurr = fullLinkData.rows(indices);
    
    arma::mat dataProp = dataCurr;
   // Rcout << "Test 2" << std::endl;
    for(int j = 0; j < dataProp.n_rows; j ++){
      if(((blockMemb(dataProp(j,4)-1) == PropBlock) + (blockMemb(dataProp(j,5)-1) == PropBlock) > 0)){
        dataProp(j,3) = 1;
      }else{
        if(dataProp(j, 4) == domain){
          dataProp(j,3) = TopicInterestDyads(domain-1, dataProp(j, 5)-1); 
        }else{
          dataProp(j,3) = TopicInterestDyads(domain-1, dataProp(j, 4)-1);  
        } 
      }  
    }
   // Rcout << "Test 3" << std::endl;
    double propLik = 0;
    //Rcout << currParams << std::endl;
    for(int j = 0; j < dataProp.n_rows; j ++){

      double propLikNum = exp(dataProp(j, 1)*currParams(0) +dataProp(j, 2)*currParams(1) + dataProp(j, 3)*currParams(2) +dataProp(j,7)*currParams(3)+dataProp(j, 8)*currParams(4));
      double propLikDem = 1 + propLikNum;
      //Rcout << propLikNum << std::endl;
      propLikNum = log(propLikNum)*dataProp(j,0)*dataProp(j, 6);
      propLikDem = log(propLikDem)*dataProp(j, 6);
      //Rcout << propLikDem << std::endl;
      propLik += propLikNum - propLikDem;
     // Rcout << propLik << std::endl;
    }
   
    //Rcout << "Test 4" << std::endl;
   blockLik(i) += propLik;
  }
  
  blockLik -= blockLik.max();
  blockLik = exp(blockLik);
  blockLik = blockLik/sum(blockLik);
  

  NumericVector blockIndex =  RcppArmadillo::sample(NumericVector(topicIndex.begin(), topicIndex.end()), 1, false, NumericVector(blockLik.begin(), blockLik.end()));
  blockIndex(0) = blockIndex(0) - 1;
  
  
  
  int PropBlock = consideredBlocks(blockIndex(0));
  
  arma::uvec indexProto = (fullLinkData.col(4) == domain) + (fullLinkData.col(5) == domain);
  
  arma::uvec indices = arma::find(indexProto);
  
  arma::mat dataCurr = fullLinkData.rows(indices);
  
  arma::mat dataProp = dataCurr;
  
  for(int j = 0; j < dataProp.n_rows; j ++){
    if(((blockMemb(dataProp(j,4)-1) == PropBlock) + (blockMemb(dataProp(j,5)-1) == PropBlock) > 0)){
      dataProp(j,3) = 1;
    }else{
      if(dataProp(j, 4) == domain){
        dataProp(j,3) = TopicInterestDyads(domain-1, dataProp(j, 5)-1); 
      }else{
        dataProp(j,3) = TopicInterestDyads(domain-1, dataProp(j, 4)-1);  
      } 
    }  
  }
  
  fullLinkData.rows(indices) = dataProp;
  
  return List::create(Named("updatedData") = fullLinkData, Named("updatedBlockMemb") = consideredBlocks(blockIndex(0)));
  
}