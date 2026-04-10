function AIC=LSRM_AIC(time,mag,reg,params,Mmin)
% out: AIC value = -2*loglik + 2*p, with p = 2J + J^2
time=time(:);
mag=mag(:);
reg=reg(:);
[J,A,B,C]=LSRM_UnpackingParam(params);
eta=0.75;
Si = 10.^(eta .* (mag - Mmin));
cmp=@lt;
% ================================
% Part 1: event log-intensity sum
% ================================
N = numel(time);
St = zeros(J, N);
for j = 1:J
    idx = (reg == j);
    if ~any(idx), continue; end
    tj = time(idx);
    Sj = Si(idx);
    I = cmp(tj, time.');          % (#events) x N
    St(j,:) = (Sj.' * I);   % 1 x N
end
% <Method 1> 
coupling=zeros(N,1);
for i=1:J
    idxi=(reg==i);
    if ~any(idxi), continue; end
    coupling(idxi)=(C(i,:)*St(:,idxi)).';
end
% lambda at each observed event
lambdai=exp( A(reg)+B(reg).*(time-coupling) );
if any(~isfinite(lambdai)) || any(lambdai <= 0)
    AIC=Inf;
    return;
end

% <Method 2> 
% sumlogi=zeros(J,1);
% for i=1:J
%     idxi=(reg==i);
%     if ~any(idxi), continue; end
%     coupling_i=(C(i,:)*St(:, idxi)).';
%     lambda_i=exp( A(i)+B(i).*(time(idxi)-coupling_i) );
%     if any(~isfinite(lambda_i)) || any(lambda_i <= 0)
%         out=Inf;
%         return;
%     end
%     sumlogi(i)=sum(log(lambda_i));
% end




% ================================
% Part 2: integrated intensity sum over regions
% ================================
T1 = min(time);
T2 = max(time);
within = (time >= T1) & (time < T2);
ti = [time(within); T2];
ti = ti(:);
nt = numel(ti);
prev = [T1; ti(1:end-1)];

% St2(j,m) = S^(j)(ti_m)
St2 = zeros(J, nt);
for j = 1:J
    idx = (reg == j);
    if ~any(idx), continue; end
    tj = time(idx);
    Sj = Si(idx);
    I = cmp(tj, ti.');          % (#events_in_j) x nt
    St2(j,:) = (Sj.' * I);      % 1 x nt
end

% integrate each region i
Int = zeros(J,1);
for i = 1:J
    ai = A(i);
    bi = B(i);
    if abs(bi) < 1e-10
        Int(i) = exp(ai) * (T2 - T1);
    else
        Kterm = exp( ai - bi * (C(i,:) * St2) );       % 1 x nt
        seg   = (exp(bi*ti.') - exp(bi*prev.'));       % 1 x nt
        Int(i) = sum( (Kterm./bi) .* seg );
    end
end
loglik=sum(log(lambdai))-sum(Int);
% number of parameters
lenparam=2*J + J^2;
AIC=-2*loglik + 2*lenparam;
end
