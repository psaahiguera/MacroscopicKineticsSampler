clc,clearvars,close all

% Parameter sampling
% RT  = 8.314*293.15/1e3;
% DGr = -10;            % kJ/mol
Keq = 3e8/6.5e5;       % From Brown & Cooper 1993 (~ -15 kJ/mol)

% Based on Bar-Even paper (kcat 1/s, Km M)
ub = [1e+4,1e+4,1e-3,1e-3,Keq];
lb = [1e+0,1e+0,1e-6,1e-6,Keq];

% Define convex space (constraints)
Aeq = [1,-1,-1,1,-1];
beq = 0;
% Aineq = [];
% bineq = [];
Aineq = [1,0,-1,0,0;...
         0,1,0,-1,0];
bineq = log([1e9;1e9]);

% Determine bounding box
[box,xinit,A,b] = calculateBoundingBox(Aeq,beq,Aineq,bineq,log(ub),log(lb));

% %% Prior selection
mu1    = log(16e-6);
sigma1 = 0.2;
mu2    = log(4827);
sigma2 = 0.2;
mu3    = log(560e-6);
sigma3 = 0.2;

% Choices
% priorFxn = @(point) uniform_prior(point);               
% priorFxn = @(point) logNormal_prior1(point,mu1,sigma1);
% priorFxn = @(point) logNormal_prior2(point,mu1,sigma1,mu2,sigma2);
priorFxn = @(point) logNormal_prior3(point,mu1,sigma1,mu2,sigma2,mu3,sigma3);

%% Sampling
% General Hit And Run Sampler
nSamples = 2e4;
nSteps   = 2e2;
nDiscard = 1e5;
points   = generalHR(A,b,box,xinit,priorFxn,nSamples,nSteps,nDiscard);

%% Plot
corrplot(points(1:4,:)','varNames',{'Lnkc+','Lnkc-','LnKs','LnKp'})