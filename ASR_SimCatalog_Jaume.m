% Have not finished yet
function Results=ASR_SimCatalog_Jaume(seed,p)
% Jaume & Bebbington (2004) stress-release model simulation with ASR design
% Output cat fields:
%   .t    event times
%   .M    magnitudes
%   .X    level just AFTER each event
%   .Si   stress release (level drop) used in X update
%   .eps  Benioff strain release (if generated)

if nargin<1 || isempty(seed)
    seed = randi([1 1e6]); 
end
rng(seed);

% ---- defaults ----
if nargin<2, p = struct(); end
p = setDefault(p,'Tlim',1e10);
p = setDefault(p,'X0',0);
p = setDefault(p,'r',0.8);

p = setDefault(p,'rateType','selfcorrecting'); % or 'poisson'
p = setDefault(p,'m',-0.47);
p = setDefault(p,'n',1.0);

p = setDefault(p,'sizeModel','TP'); % 'TP' tapered Pareto or 'TGR'
p = setDefault(p,'alpha',1.0);
p = setDefault(p,'y0',10^(2.4+0.75*4.0)); % consistent with Mmin=4 idea
p = setDefault(p,'m0',4.0);

p = setDefault(p,'gType','const'); % 'const'|'exp'|'pow'
p = setDefault(p,'gConst',7.0);
p = setDefault(p,'u0',5.0);
p = setDefault(p,'u1',0.1);
p = setDefault(p,'u2',2.0);

p = setDefault(p,'Mmin',4.0);
p = setDefault(p,'Mmax',8.5);
p = setDefault(p,'bValue',1.0);

p = setDefault(p,'levelDropType','benioff'); % 'benioff' or 'moment'

% ---- prealloc ----
maxEvents = 2e6;
t  = nan(1,maxEvents);
M  = nan(1,maxEvents);
X  = nan(1,maxEvents);
Si = nan(1,maxEvents);
epsi = nan(1,maxEvents);

% ---- init ----
ti = 0.0;
Xi = p.X0;
k  = 0;

while ti <= p.Tlim && k < maxEvents
    % 1) sample waiting time to next event
    switch lower(p.rateType)
        case 'poisson'
            lambda0 = exp(p.m);
            tau = exprnd(1/lambda0);
        case 'selfcorrecting'
            % Between events, X(t) = Xi + r*s, s in [0,tau]
            % lambda(s)=exp(m+n*(Xi+r*s)) = lambda_i * exp(n*r*s)
            lambda_i = exp(p.m + p.n*Xi);
            nr = p.n * p.r;
            U1 = rand();
            if abs(nr) < 1e-12
                % approx constant intensity
                tau = -log(U1)/lambda_i;
            else
                tau = (1/nr) * log( 1 + (nr/lambda_i) * (-log(U1)) );
            end
        otherwise
            error('Unknown par.rateType');
    end

    % advance to event time, compute level just before event
    ti = ti + tau;
    if ti > p.Tlim, break; end
    X_before = Xi + p.r * tau;

    % 2) decide g(X) (upper turning magnitude / max magnitude control)
    g = compute_g(p, X_before);

    % 3) sample event size
    switch upper(p.sizeModel)
        case 'TP'
            % Generate epsilon ~ tapered Pareto by: epsilon = min(Pareto, shifted-Exp)
            % Pareto: P(Yp>y)= (y/y0)^(-alpha), y>=y0
            U2 = rand();
            Yp = p.y0 * (1 - U2)^(-1/p.alpha);

            % Exponential taper: exp(-(y-y0)/U) => y = y0 - U*log(U3)
            % Need U (in epsilon units). Convert g (magnitude) -> U via (8):
            % M = m0 + (4/3) log10(eps) => eps = 10^((3/4)(M-m0))
            Ueps = 10.^((3/4)*(g - p.m0));  % "equivalent magnitude" mapping
            U3 = rand();
            Ye = p.y0 - Ueps * log(U3);

            eps_event = min(Yp, Ye);
            Mi = p.m0 + (4/3)*log10(eps_event);

        case 'TGR'
            % Truncated GR on [Mmin, mmax], with mmax = g
            mmax = min(g, p.Mmax);
            Mi = sampleTruncGR(p.Mmin, mmax, p.bValue);
            eps_event = 10.^((3/4)*(Mi - p.m0)); % if you want an epsilon proxy
        otherwise
            error('Unknown par.sizeModel');
    end

    % 4) compute stress release Si (level drop)
    switch lower(p.levelDropType)
        case 'benioff'
            % Si = 10^(2.4 + 0.75*M)
            Si_i = 10^(2.4 + 0.75*Mi);
        case 'moment'
            % Si = 10^(9.0 + 1.5*M)
            Si_i = 10^(9.0 + 1.5*Mi);
        otherwise
            error('Unknown par.levelDropType');
    end

    % 5) update level after event (jump down)
    Xi = X_before - Si_i;

    % store
    k = k + 1;
    t(k) = ti;
    M(k) = Mi;
    X(k) = Xi;
    Si(k) = Si_i;
    epsi(k) = eps_event;
end

% trim
Results.t   = t(1:k);
Results.M   = M(1:k);
Results.X   = X(1:k);
Results.Si  = Si(1:k);
Results.eps = epsi(1:k);
Results.seed = seed;
Results.par  = p;
end

% ---------- SubFunctions ----------
function par = setDefault(par, name, val)
if ~isfield(par,name) || isempty(par.(name))
    par.(name) = val;
end
end

function g = compute_g(par, X)
switch lower(par.gType)
    case 'const'
        g = par.gConst;
    case 'exp'
        g = par.u0 + par.u1 * exp(X);                  % (9)
    case 'pow'
        Xp = max(X,0);
        g = par.u0 + par.u1 * (1 + Xp)^(par.u2);       % (10)
    otherwise
        error('Unknown par.gType');
end
end

function M = sampleTruncGR(Mmin, Mmax, b)
beta = b*log(10);
K = 1 / (1 - exp(-beta*(Mmax - Mmin)));
M = Mmin - (1/beta) * log(1 - rand() / K);
end
