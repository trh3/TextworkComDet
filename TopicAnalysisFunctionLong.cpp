//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List likeliLong(NumericMatrix docs_s, NumericMatrix wordDist_s, NumericMatrix origWordDist_s, NumericVector membNum_s, NumericVector origMemNum_s,
NumericVector assign_s, NumericVector domains, NumericMatrix topicLikeli_s, NumericVector priors){
	NumericMatrix docs = Rcpp::clone(docs_s);
	NumericMatrix wordDist = Rcpp::clone(wordDist_s);
	NumericMatrix WordDist2 = Rcpp::clone(wordDist_s);
	NumericMatrix origWordDist2 = Rcpp::clone(origWordDist_s);
	NumericVector assign = Rcpp::clone(assign_s);
	NumericVector membNum = Rcpp::clone(membNum_s);
	NumericVector membNum2 = Rcpp::clone(membNum_s);
	NumericVector membNumLag = Rcpp::clone(origMemNum_s);
	int nDocs = docs.nrow(), nWords = docs.ncol(), nTopics = wordDist.nrow(), nDomains = topicLikeli_s.ncol();
	
	arma::mat Docs(docs.begin(), nDocs, nWords, false);
	arma::mat Topics(wordDist.begin(), nTopics, nWords, false);
	arma::mat TopicsLag(origWordDist2.begin(), nTopics, nWords, false);
	arma::mat TopicsOrig(WordDist2.begin(), nTopics, nWords, false);
	arma::vec MembNum(membNum.begin(), nTopics, false);
	arma::vec MembNumOrig(membNum2.begin(), nTopics, false);
	arma::vec MembNumLag(membNumLag.begin(), nTopics, false);
	arma::vec Priors(priors.begin(), 5, false);
	arma::vec topicIndex = arma::linspace(1, nTopics, nTopics);
	arma::mat topicLikeli(topicLikeli_s.begin(), nTopics,  nDomains, false);
  
	NumericVector TopicIndex = NumericVector(topicIndex.begin(), topicIndex.end());
	
	arma::vec Assign(assign.begin(), nDocs, false);
	for(int i = 0; i < nDocs; i++){
		Topics.row(Assign(i)-1) += Docs.row(i);
		MembNum(Assign(i)-1) += 1;
	}
	
	for(int i = 0; i < nDocs; i++){
		int currAssign = Assign(i);
		Topics.row(currAssign-1) -= Docs.row(i);
		MembNum(currAssign-1) -= 1;
		arma::vec doc = arma::vectorise(Docs.row(i));
		double nDocWords = arma::sum(doc);
		arma::vec prob = arma::zeros<arma::vec>(nTopics);
		arma::vec front = arma::log(MembNum + Priors(0)) - log(nDocs - 1 + Priors(0)*nTopics);
		for(int j = 0; j < nTopics; j++){
			 double denom = 0.0;
			 double numTopicWords = arma::sum(Topics.row(j));
			for(int words = 1; words <= nDocWords; words++){
				denom += log(numTopicWords + nWords*Priors(1) + words -1) ;
			}
			double num = 0.0;
			for(int w = 0; w < nWords; w++){

				if(doc(w) > 0){
					for(int wordCount = 1; wordCount <= doc(w); wordCount++){
					num += log(Topics(j,w) + Priors(1) + wordCount - 1);
					}
				}
			}
		prob(j) = num - denom;

			
			
		}
    arma::vec uprob = prob;
    prob += topicLikeli.col(domains[i]-1);
	prob += front;
		double probMax = prob.max();
		double uprobMax = uprob.max();
		for(int p = 0; p < nTopics; p++){
			prob(p) = arma::trunc_exp(prob(p) - probMax);
				uprob(p) = arma::trunc_exp(uprob(p) - uprobMax);
		} 
		
	NumericVector ass; 

	arma::vec nprob = prob/arma::sum(prob);
	arma::vec unprob = uprob/arma::sum(uprob);
 	
	NumericVector propa = NumericVector(nprob.begin(), nprob.end());
    if(arma::is_finite(nprob)){
		ass = RcppArmadillo::sample(TopicIndex , 1, false, NumericVector(nprob.begin(), nprob.end()));
      
    }else{
		NumericVector uproba = NumericVector(unprob.begin(), unprob.end());
		ass =  RcppArmadillo::sample(TopicIndex , 1, false, uproba);
	}

   
    
		int assi = ass[0];
		Assign(i) = assi;
    
//    if(assi != currAssign){
//      
//      Rcout << "Switch!" << std::endl;
//    }
    Topics.row(assi-1) += Docs.row(i);
		MembNum(assi-1) += 1;
	
	}
	
	arma::vec MembNumReturn = MembNum - MembNumLag;
    arma::vec MembNumDayReturn = MembNum - MembNumOrig;

	arma::mat TopicsReturn = Topics - TopicsLag;
    arma::mat TopicsDayReturn = Topics - TopicsOrig;
    
		for(int i = 0; i < nTopics; i++){
			if(MembNumReturn(i) == 0){
				TopicsReturn.row(i) == TopicsOrig.row(i);
			}
		}
		
    List forReturn;
    forReturn["assign"] = Assign;
    forReturn["Topic.Matrix"] = TopicsReturn;
    forReturn["MemberNum"] = MembNumReturn;
    forReturn["DayTopicMatrix"] = TopicsDayReturn;
    forReturn["DayMembNum"] = MembNumDayReturn;
        
		return forReturn;
	
}