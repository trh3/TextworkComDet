//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List likeliLongInit(NumericMatrix docs_s, NumericMatrix wordDist_s, NumericVector assign_s, NumericVector membNum_s, NumericVector priors){
	NumericMatrix docs = Rcpp::clone(docs_s);
	NumericMatrix wordDist = Rcpp::clone(wordDist_s);
	NumericMatrix wordDist2 = Rcpp::clone(wordDist_s);
	NumericVector assign = Rcpp::clone(assign_s);
	NumericVector membNum = Rcpp::clone(membNum_s);
	NumericVector membNum2 = Rcpp::clone(membNum_s);
	int nDocs = docs.nrow(), nWords = docs.ncol(), nTopics = wordDist.nrow();
	

	

	arma::mat Docs(docs.begin(), nDocs, nWords, false);
	
		arma::mat wordSums = arma::sum(Docs, 0);
	arma::umat wordSumsBool = (wordSums > 0);
	
	int realNWords = arma::as_scalar(arma::sum(wordSumsBool,1));
	arma::mat Topics(wordDist.begin(), nTopics, nWords, false);
	arma::mat TopicsOrig(wordDist2.begin(), nTopics, nWords, false);
	arma::vec MembNum(membNum.begin(), nTopics, false);
	arma::vec MembNumOrig(membNum2.begin(), nTopics, false);
	arma::vec Priors(priors.begin(), 5, false);
	arma::vec topicIndex = arma::linspace(1, nTopics, nTopics);
	NumericVector TopicIndex = NumericVector(topicIndex.begin(), topicIndex.end());
	
	NumericVector assignRand = RcppArmadillo::sample(TopicIndex , nDocs, true);
	arma::vec Assign(assignRand.begin(), nDocs, false);
	for(int i = 0; i < nDocs; i++){
		Topics.row(Assign(i)-1) += Docs.row(i);
		MembNum(Assign(i)-1) += 1;
	}
	
	for(int iter = 1; iter <= priors(3); iter++){
	for(int i = 0; i < nDocs; i++){
		int currAssign = Assign(i);
		Topics.row(currAssign-1) -= Docs.row(i);
		MembNum(currAssign-1) -= 1;
		arma::vec doc = arma::vectorise(Docs.row(i));
		int nDocWords = arma::sum(doc);
		arma::vec prob = arma::zeros<arma::vec>(nTopics);
		arma::vec front = (MembNum + Priors(0)) / (nDocs - 1 + Priors(0)*Priors(2));
		for(int j = 0; j < nTopics; j++){
			double denom = 0.0;
			int numTopicWords = arma::sum(Topics.row(j));
			for(int words = 1; words <= nDocWords; words++){
				denom += log(numTopicWords + realNWords*Priors(1) + words -1) ;
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
		for(int p = 0; p < nTopics; p++){
			prob(p) += log(front(p));
		}
		double probMax = prob.max();
		
		for(int p = 0; p < nTopics; p++){
			prob(p) = exp(prob(p) - probMax);
		} 
		 
		arma::vec nprob = prob/arma::sum(prob);
    NumericVector ass = RcppArmadillo::sample(TopicIndex , 1, false, NumericVector(nprob.begin(), nprob.end()));
		int assi = ass[0];
		Assign(i) = assi;
    Topics.row(assi-1) += Docs.row(i);
		MembNum(assi-1) += 1;
	
	}
	}
	
		arma::vec MembNumReturn = MembNum - MembNumOrig;
		arma::mat TopicsReturn = Topics;
		
		for(int i = 0; i < nTopics; i++){
			if(MembNumReturn(i) == 0){
				TopicsReturn.row(i) == TopicsOrig.row(i);
			}
		}
		
		return List::create(Named("assign") = Assign,
Named("Topic.Matrix") = TopicsReturn, Named("MemberNum") = MembNumReturn);
	
}