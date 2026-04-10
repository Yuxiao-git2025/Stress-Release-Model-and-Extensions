function [pNew,AICNew]=LSRM_MLE(time,mag,reg,param0,Mmin)
Opt = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'iter', ...
    'MaxIterations', 150, ...
    'StepTolerance', 1e-5, ...
    'ConstraintTolerance', 1e-5);
time=time(:);
mag=mag(:);
reg=reg(:);
n=numel(param0);
J=round((-2+sqrt(4+4*n))/2);   % solves J^2 + 2J - n = 0
FuncHandel=@(param) LSRM_AIC(time,mag,reg,param,Mmin);
% Set bounds
lb=[-10 * ones(J,1);   1e-10 * ones(J,1);  -1 * ones(J*J,1)];
ub=[ 10 * ones(J,1);    1    * ones(J,1);   1 * ones(J*J,1)];
% It is convenient, if ignoring aftershocks, to set θii=1
[pNew,AICNew]=fmincon(FuncHandel,param0,[],[],[],[],lb,ub,[],Opt);
pNew=pNew(:)';
disp('# New parameters');
disp(pNew);

end
