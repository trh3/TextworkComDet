//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//[[Rcpp::export]]
List likeliLongInit(NumericMatrix docs, NumericMatrix wordDist, NumericVector assign, NumericVector membNum, NumericVector priors){
	
	int nDocs = docs.nrow(), nWords = docs.ncol(), nTopics = wordDist.nrow();
	
	arma::mat Docs(docs.begin(), nDocs, nWords, false);
	arma::mat Topics(wordDist.begin(), nTopics, nWords, false);
	arma::mat TopicsOrig(wordDist.begin(), nTopics, nWords, false);
	arma::vec Assign(assign.begin(), nDocs, false);
	arma::vec MembNum(membNum.begin(), nTopics, false);
	arma::vec MembNumOrig(membNum.begin(), nTopics, false);
	arma::vec Priors(priors.begin(), 3, false);
	arma::vec topicIndex = arma::linspace(1, nTopics, nTopics);
	NumericVector TopicIndex = NumericVector(topicIndex.begin(), topicIndex.end());
	for(int i = 0; i < nDocs; i++){
		arma::vec doc = arma::vectorise(Docs.row(i));
		int nDocWords = arma::sum(doc);
		arma::vec prob = arma::zeros<arma::vec>(nTopics);
		arma::vec front = (MembNum + Priors(0)) / (nDocs - 1 + Priors(0)*Priors(2));
		for(int j = 0; j < nTopics; j++){
			double denom = 0.0;
			int numTopicWords = arma::sum(Topics.row(j));
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
	
		return List::create(Named("assign") = Assign,
Named("Topic.Matrix") = Topics - TopicsOrig, Named("MemberNum") = MembNum - MembNumOrig);
	
}