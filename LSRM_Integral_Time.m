function [CumInt]=LSRM_Integral_Time(time,mag,reg,params,Mmin)
% Compute lambda_i(t_k) and integrated intensity ∫ lambda_i over time
    T1 = min(time);
    time = time(:);
    mag  = mag(:);
    reg  = reg(:);
    [time, ord]=sort(time, 'ascend');
    mag = mag(ord);
    reg = reg(ord);

    [J, A, B, C]=LSRM_UnpackingParam(params);
    A = A(:); 
    B = B(:);
    eta = 0.75;
    Si = 10.^(eta .* (mag - Mmin));
    cmp = @lt;
    N = numel(time);
    prev = [T1; time(1:end-1)];

    % -----------------------
    % 1) Stress histories at each event time: St(j,k)=S^j(time_k)
    % -----------------------
    St = zeros(J, N);
    for j = 1:J
        idx = (reg == j);
        if ~any(idx), continue; end
        tj = time(idx);
        Sj = Si(idx);
        I  = cmp(tj, time.');
        St(j,:) = (Sj.' * I); 
    end
    % coupling(i,k)=sum_j [C(i,j)*St(j,k)]
    coupling = C * St;

    % -----------------------
    % 2) Lambda at each event time (per region): J x N
    % -----------------------
    Lambda = exp( A + B .* (time.' - coupling) );
    % -----------------------
    % 3) Piecewise integral over [prev(k), time(k)] for each region i
    %    On interval k
    % -----------------------
    dt=(time-prev);
    dInt=zeros(J,N);
    for i=1:J
        ai=A(i);
        bi=B(i);
        % K(i,k) = exp(ai - bi*coupling(i,k))  (1 x N)
        K=exp(ai-bi .* coupling(i,:));
        if abs(bi)<1e-10
            % lambda_i(t) = exp(ai - bi*coupling) * exp(bi*t) ~ exp(ai) if bi==0
            dInt(i,:)=exp(ai).*dt.';    % 1 x N
        else
            % Stable form using expm1 to reduce cancellation for small bi*dt
            % ∫_{prev}^{t} K * exp(bi*u) du = K * exp(bi*prev) * (expm1(bi*dt))/bi
            dInt(i,:)=(K .* exp(bi .* prev.')) .* (expm1(bi .* dt.')./bi);
        end
    end
    % cumulative integral
    CumInt=cumsum(dInt, 2);
end
