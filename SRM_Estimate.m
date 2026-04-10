function out=SRM_Estimate(time,mag,Mmin,param,IntTime)
% Stress-release model conditional intensity or its integral
% Modified at: 2026/03/23
global evalNums
time=time(:); 
mag=mag(:);
A=param(1);
B=param(2);
C=param(3);
if nargin<5
    IntTime=false;
end
% evalpts: evaluating vector of times
evalpts=linspace(min(time),max(time),evalNums)';
evalTimes=evalpts(:);
eta=0.75;
Si=10.^(eta.*(mag(:)-Mmin));   % S(i)=10^(eta*(M(i)-Mmin))

if ~IntTime
    % >> Conditional intensity at eval times 
    % This is aconvenient compromise between time-predictable and purely 
    % random (Poisson) processes
    Ind=(time(:)<evalTimes(:).');
    St=Si.'*Ind;
    out=exp(A+B.*(evalTimes.'-C .* St));
    out=out(:);
else
    % >> Integral the intensity
    T1=time(1);
    T2=time(end);
    ti=time;
    Ind=time < ti.';         % N x nt
    St=Si.'*Ind;         
    if abs(B)<1e-12
        % limitation: lambda(t)=exp(A), integral=exp(A)*(T2-T1)
        out2 = exp(A) * (T2 - T1);
    else % integral for exponential-type function exp(B*t)
        prev = [T1; ti(1:end-1)];
        out2 = sum( exp(A - B*C .* St) ./ B .* (exp(B .* ti.') - exp(B .* prev.')) );
    end
    out1=sum(A + B .* (ti.' - C .* St));
    % Maximizing the log-likelihood (MLE)
    loglik=out1-out2;
    % AIC value
    out=-2*loglik+2*length(param);
    return;
end
end
