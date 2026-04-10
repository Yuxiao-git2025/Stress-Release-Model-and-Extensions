function [EventTime,EventMag,E,OutEnergy,OutSumEnergy]=...
    ASR_SimCatalog_Barani(seed)
if nargin<1 || isempty(seed)
    seed = randi([1 1e4]);
end
rng(seed);
global Mmin Mmax Bvalue Tlim A B C eta
maxEvents = 1e6;
EventTime = nan(1,maxEvents);
EventMag  = nan(1,maxEvents);
OutEnergy = nan(1,maxEvents);   
E         = nan(1,maxEvents); 

% ---- initial ----
t  = 0.0;
S  = 0.0;
n  = 0;
OutSumEnergy = 0.0;
beta = Bvalue*log(10);
K = 1/(1 - exp(-beta*(Mmax - Mmin)));

% ---- simulate event by event ----
while t <= Tlim && n < maxEvents
    % state-dependent rate parameter
    A0 = A + B*(t - C*S);
    R0 = exp(A0)/B;

    % sample waiting time tau using inverse transform
    U1 = rand();
    tau = (1/B) * log( 1 + (-log(U1))/R0 );

    % advance time
    t = t + tau;
    if t > Tlim
        break
    end

    % sample magnitude from truncated GR
    M  = -1/beta * log(1 - rand()/K) + Mmin;

    % relative stress release
    Si = 10^(eta*(M - Mmin));

    % update cumulative stress
    S = S + Si;

    % store
    n = n + 1;
    EventTime(n) = t;
    EventMag(n)  = M;
    OutEnergy(n) = Si;
    E(n)         = S;
    OutSumEnergy = OutSumEnergy + Si;
end

% trim outputs
EventTime = EventTime(1:n);
EventMag  = EventMag(1:n);
OutEnergy = OutEnergy(1:n);
E         = E(1:n);


fprintf('# events=%d, OutSumEnergy=%1.2e\n', n, OutSumEnergy);
end
