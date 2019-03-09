function prob = logNormal_prior1(x,mu,sigma)
prob = log( normpdf(x(3),mu,sigma) );