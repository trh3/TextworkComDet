######Daily Activation Scripts
sourceCpp("DailyActivationSampler.cpp")

dailyActivationFullSampler <- function(processedData, exciteParams, domainRates, priors, n.days, n.domains, n.topics){

	dailyActivation <- dailyActivationSampler(processedData, exciteParams, domainRates, c(n.days, n.domains,n.topics),  priors )

	return(dailyActivation)
}



