function prob = logNormal_prior3(x,mu1,sigma1,mu2,sigma2,mu3,sigma3)
prob = log( normpdf(x(3),mu1,sigma1)) + log( normpdf(x(1),mu2,sigma2) ) + log( normpdf(x(4),mu3,sigma3) );