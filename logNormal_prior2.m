function prob = logNormal_prior2(x,mu1,sigma1,mu2,sigma2)
prob = log( normpdf(x(3),mu1,sigma1)) + log( normpdf(x(1),mu2,sigma2) );