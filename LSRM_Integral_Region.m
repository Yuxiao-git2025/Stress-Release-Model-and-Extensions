function [CumInt]=LSRM_Integral_Region(time, mag, reg, params, Mmin)
% Return integrated intensity evaluated at each region's own event times
T1 = min(time);
time = time(:); mag = mag(:); reg = reg(:);
% sort globally by time (important for piecewise integration)
[time, ord] = sort(time, 'ascend');
mag = mag(ord);
reg = reg(ord);
[J, A, B, C] = LSRM_UnpackingParam(params);
A = A(:); 
B = B(:);
eta = 0.75;
Si  = 10.^(eta .* (mag - Mmin));
cmp = @lt;

N = numel(time);
prev = [T1; time(1:end-1)];
dt   = time - prev;
% ---- 1) Stress histories at each global event time: St(j,k)=S^j(time_k)
St = zeros(J, N);
for j = 1:J
    idx = (reg == j);
    if ~any(idx), continue; end
    tj = time(idx);
    Sj = Si(idx);
    I  = cmp(tj, time.');
    St(j,:) = (Sj.' * I);
end

% coupling(i,k)=sum_j C(i,j)*St(j,k)
coupling = C * St; 

% ---- 2) Compute per-interval integrated increments for each region
dInt_all = zeros(J, N);

for i = 1:J
    ai = A(i);
    bi = B(i);

    % K(i,k) = exp(ai - bi*coupling(i,k))
    K = exp(ai - bi .* coupling(i,:));

    if abs(bi) < 1e-10
        % If bi==0, lambda_i(t)=exp(ai) (coupling drops out because multiplied by bi)
        dInt_all(i,:) = exp(ai) .* dt.';
    else
        % ∫_{prev}^{t} K*exp(bi*u) du = K*exp(bi*prev) * expm1(bi*dt)/bi
        dInt_all(i,:) = (K .* exp(bi .* prev.')) .* (expm1(bi .* dt.')./bi);
    end
end

% ---- 3) Convert to cumulative integral on the global grid
CumInt_all = cumsum(dInt_all, 2);          % J x N

% ---- 4) Pick only each region's own event times, return variable-length cells
CumInt = cell(J,1);
dTau   = cell(J,1);
tR      = cell(J,1);

for i = 1:J
    idxi = (reg == i);
    t_i  = time(idxi);
    tR{i} = t_i;

    if isempty(t_i)
        CumInt{i} = [];
        dTau{i}   = [];
        continue;
    end
    CI_i = CumInt_all(i, idxi).';
    CumInt{i} = CI_i;
    % increments between successive events in region i (for residual tests)
    dTau{i} = [CI_i(1); diff(CI_i)];
end
end
