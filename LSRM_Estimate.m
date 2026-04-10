% =========================================================================
% LSRM model and its conditional intensity
% Refer to: (Bebbington & Harte 2003)
% params:  [a1..aJ, b1..bJ, vec(C)] where C is JxJ
% eta:     default 0.75
% =========================================================================
function [outi]=LSRM_Estimate(time, mag, reg, params, Mmin)
global evalNums
if isempty(evalNums)
    evalNums=200;
end
time=time(:);
mag=mag(:);
reg=reg(:);
[J,A,B,C]=LSRM_UnpackingParam(params);

eta=0.75;
Si=10.^(eta.*(mag-Mmin));
cmp=@lt; % i.e. < rather than <=

evalTime=linspace(min(time),max(time),evalNums).';
Ne = numel(evalTime);
if isempty(time) || isempty(mag)
    outi = exp(A + B .* evalTime.');
    return;
else
    % St(j,k) = S^(j)(evalTime_k)
    St=zeros(J, Ne);
    for j=1:J
        idx=(reg==j);
        if ~any(idx), continue; end
        tj=time(idx);
        Sj=Si(idx);
        I=cmp(tj,evalTime.');   % Num x Ne
        St(j,:)=(Sj.'*I);
    end
    % Calculate the intensity of different regional conditions using matrix
    % operations:
    % coupling(i,k)=sum_j[C(i,j) * St(j,k)]
    coupling=C*St;  % Transfered stress
    % lambda(i,k)=exp[ a(i)+b(i)*(ti(k)-coupling(i,k)) ]
    outi=exp(A+B.*(evalTime.'-coupling));  % conditional intensity
    return;
end
end
