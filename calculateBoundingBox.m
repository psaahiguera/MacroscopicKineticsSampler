function [box,xinit,Aeq,beq] = calculateBoundingBox(Aeq,beq,Aineq,bineq,ub,lb)

% Determine if there are inequalities. If so, add slack variables
if ~isempty(Aineq) && (size(Aineq,1)==numel(bineq))
    numEqConstraints   = size(Aeq,1);
    numIneqConstraints = size(Aineq,1);
    Aeq = [Aeq,zeros(numEqConstraints,numIneqConstraints);...
            Aineq,eye(numIneqConstraints)];
    beq = [beq;bineq];
    ub  = [ub(:);1e2*ones(numIneqConstraints,1)];   % add bounds on the slacks (if appropriate)
    lb  = [lb(:);zeros(numIneqConstraints,1)];
end

% Determine bounding box and initial point for the Hit-And-Run
options.Display = 'off';
box  = zeros(size(Aeq,2),2);
fobj = zeros(size(Aeq,2),1);
x0   = [];
for ix = 1:size(Aeq,2)
    f = fobj;
    f(ix) = 1;
    xmin  = linprog(f,[],[],Aeq,beq,lb,ub,[],options);
    f(ix) = -1;
    xmax  = linprog(f,[],[],Aeq,beq,lb,ub,[],options);
    box(ix,:) = [xmin(ix),xmax(ix)];
    x0 = [x0,xmin,xmax];
end
xinit = mean(x0,2);
