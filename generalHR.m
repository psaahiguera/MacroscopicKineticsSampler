function points = generalHR(Aeq,beq,box,xinit,priorFxn,nSamples,nSteps,nDiscard)

% Determine the dimension of the feasible space
A  = Aeq;
UB = box(:,2);
LB = box(:,1);
ix_fixed = (abs(UB-LB)==0);
if any(ix_fixed)
    Aeq(:,ix_fixed) = 0;
end
[~,colsInd] = rref(Aeq);
colsDep = find(~ismember(1:size(Aeq,2),colsInd)&~all(Aeq==0,1));
fVars   = numel(colsDep);
invAind = pinv(A(:,colsInd));
linearMap = @(x) invAind*(beq - A(:,colsDep)*x);

% General Hit And Run Sampler
colA   = size(A,2);
uTol   = 1e-9;
bTol   = 1e-9;
currPoint = xinit;                 % Initial point
Ln_currProb = priorFxn(currPoint);   % Log prior
points = zeros(colA,nSamples);
count = 0;
iter  = 0;
while (count<=nSamples)
    iter = iter+1;
    
    while true
        
        % Draw random direction
        udir = zeros(colA,1);
        udir(colsDep) = randn(fVars,1);
        udir(colsInd) = linearMap(udir(colsDep));
        udir = udir/norm(udir);                         % unitary direction
        
        % Figure out max distance
        posDir      = (udir>uTol);
        negDir      = (udir<-uTol);
        posStepTemp = (UB-currPoint)./udir;
        negStepTemp = (LB-currPoint)./udir;
        posStep     = min([posStepTemp(posDir);negStepTemp(negDir)]);
        negStep     = max([negStepTemp(posDir);posStepTemp(negDir)]);
        
        % Draw new direction if too close to the boundaries
        Lcord = posStep-negStep;
        flag  = (Lcord < bTol) | (posStep < 0) | (negStep > 0);
        if ~flag; break; end
    end
    
    % Sample step size and perform random jump
    lambda    = negStep + rand(1)*Lcord;
    nextPoint = currPoint + lambda*udir;
    
    % Accept/Reject proposal using MCMC kernel
    Ln_nextProb = priorFxn(nextPoint);
    Ln_acceptProb  = Ln_nextProb - Ln_currProb;
    if (Ln_acceptProb < log(rand(1)))
        nextPoint = currPoint;
    else
        Ln_currProb = Ln_nextProb;
    end
    
    % Update current point
    currPoint = nextPoint;
    
    % Save next point
    if (nDiscard <= iter) && ~mod(iter-nDiscard,nSteps)
        count = count + 1;
        points(:,count) = nextPoint;
    end
end