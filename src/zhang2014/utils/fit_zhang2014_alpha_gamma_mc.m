function fit = fit_zhang2014_alpha_gamma_mc(obs, vdMethod, wsMethod, varargin)
%FIT_ZHANG2014_ALPHA_GAMMA_MC Fit alpha and gamma with measurement-error MC.
%
% The reported uncertainty is interpreted as a symmetric 95% confidence
% half-width on wd. Each MC realization perturbs all observations, then fits
% alpha and gamma simultaneously using an unweighted log-space objective.
%
% Required table variables:
%   case_key, ustar_mps, z0_m, zr_m, Dp_m, wd_obs_mps,
%   sigma_mps, T_K, p_Pa, L_m, rhop
%
% Optional Name/Value:
%   'nMC'          default 500
%   'nAcceptedTarget' default [] (when set, continue until this many fits
%                     converge successfully)
%   'maxTrials'    default [] (10*nAcceptedTarget when a target is set)
%   'alphaBounds'  default [1e-3 1e2]
%   'gammaBounds'  default [1e-3 1e1]
%   'rngSeed'      default 1
%   'verbose'      default true

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'nMC',500,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'nAcceptedTarget',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0));
addParameter(p,'maxTrials',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0&&mod(x,1)==0));
addParameter(p,'alphaBounds',[1e-3 1e2],@(x)isnumeric(x)&&numel(x)==2&&all(x>0)&&x(2)>x(1));
addParameter(p,'gammaBounds',[1e-3 1e1],@(x)isnumeric(x)&&numel(x)==2&&all(x>0)&&x(2)>x(1));
addParameter(p,'rngSeed',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'verbose',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;

required = ["case_key","ustar_mps","z0_m","zr_m","Dp_m","wd_obs_mps", ...
            "sigma_mps","T_K","p_Pa","L_m","rhop"];
missing = setdiff(required,string(obs.Properties.VariableNames));
if ~isempty(missing)
    error('Missing required observation variables: %s',strjoin(missing,', '));
end

good = isfinite(obs.wd_obs_mps) & obs.wd_obs_mps>0 & ...
       isfinite(obs.Dp_m) & obs.Dp_m>0;
obs = obs(good,:);

if isempty(obs)
    error('No valid observations supplied.');
end

rng(opt.rngSeed,'twister');

ab = opt.alphaBounds(:)';
gb = opt.gammaBounds(:)';

pack = @(a,g) [bounded_inverse(a,ab); bounded_inverse(g,gb)];
unpack = @(x) deal(bounded_forward(x(1),ab), bounded_forward(x(2),gb));

alpha0 = min(max(3,ab(1)*(1+1e-6)),ab(2)*(1-1e-6));
gamma0 = min(max(0.5,gb(1)*(1+1e-6)),gb(2)*(1-1e-6));
x0 = pack(alpha0,gamma0);

optimOptions = optimset( ...
    'Display','off', ...
    'TolX',1e-2, ...
    'TolFun',1e-4, ...
    'MaxIter',500, ...
    'MaxFunEvals',2000);

objNom = @(x)objective_log_median(x,obs,vdMethod,wsMethod,unpack);
[xNom,sseNom,exitNom] = fminsearch(objNom,x0,optimOptions);
[aNom,gNom] = unpack(xNom);

fit = struct();
fit.nominal = struct('alpha',aNom,'gamma',gNom,'SSE',sseNom, ...
                     'exitflag',exitNom,'nObs',height(obs), ...
                     'nDiam',numel(unique(round(obs.Dp_m*1e6,6))));
fit.options = opt;

if isempty(opt.nAcceptedTarget)
    nAcceptedTarget = [];
    maxTrials = opt.nMC;
else
    nAcceptedTarget = opt.nAcceptedTarget;
    if isempty(opt.maxTrials)
        maxTrials = max(10*nAcceptedTarget,nAcceptedTarget+100);
    else
        maxTrials = opt.maxTrials;
    end
    if maxTrials<nAcceptedTarget
        error('maxTrials must be at least nAcceptedTarget.');
    end
end

aMC = nan(maxTrials,1);
gMC = nan(maxTrials,1);
sseMC = nan(maxTrials,1);
exitMC = nan(maxTrials,1);
trialIndex = (1:maxTrials)';

% Convert reported 95% half-width to 1-sigma uncertainty.
sigma1 = obs.sigma_mps ./ 1.96;
validSigma = isfinite(sigma1) & sigma1>=0;
sigma1(~validSigma) = 0;

% Lognormal perturbation with median equal to the reported wd.
sigmaLog = sqrt(log1p((sigma1 ./ obs.wd_obs_mps).^2));
sigmaLog(~isfinite(sigmaLog)) = 0;

xStart = xNom;

imc = 0;
nAccepted = 0;
while imc<maxTrials
    imc = imc+1;
    obsMC = obs;
    obsMC.wd_obs_mps = exp(log(obs.wd_obs_mps) + sigmaLog .* randn(height(obs),1));

    objMC = @(x)objective_log_median(x,obsMC,vdMethod,wsMethod,unpack);
    [xHat,fval,eflag] = fminsearch(objMC,xStart,optimOptions);
    [aMC(imc),gMC(imc)] = unpack(xHat);
    sseMC(imc) = fval;
    exitMC(imc) = eflag;

    acceptedThisFit = isfinite(aMC(imc)) && isfinite(gMC(imc)) && ...
        isfinite(sseMC(imc)) && exitMC(imc)>0;
    nAccepted = nAccepted+acceptedThisFit;

    % Warm start from previous fit.
    xStart = xHat;

    if opt.verbose && (imc==1 || mod(imc,100)==0 || ...
            (~isempty(nAcceptedTarget) && nAccepted>=nAcceptedTarget) || ...
            (isempty(nAcceptedTarget) && imc==maxTrials))
        if isempty(nAcceptedTarget)
            fprintf('MC fit %d of %d\n',imc,maxTrials);
        else
            fprintf('MC attempts %d; accepted %d of %d\n', ...
                imc,nAccepted,nAcceptedTarget);
        end
    end

    if ~isempty(nAcceptedTarget) && nAccepted>=nAcceptedTarget
        break
    end
end

aMC = aMC(1:imc);
gMC = gMC(1:imc);
sseMC = sseMC(1:imc);
exitMC = exitMC(1:imc);
trialIndex = trialIndex(1:imc);

fit.mc = table(trialIndex,aMC,gMC,sseMC,exitMC, ...
    'VariableNames',{'trialIndex','alpha','gamma','SSE','exitflag'});

accepted = isfinite(aMC) & isfinite(gMC) & isfinite(sseMC) & exitMC>0;
fit.accepted = accepted;
fit.nAccepted = nnz(accepted);
fit.nAttempted = imc;

if ~isempty(nAcceptedTarget) && fit.nAccepted<nAcceptedTarget
    error(['Only %d of the requested %d fits were accepted after %d attempts. ' ...
        'Increase maxTrials or inspect the failed optimizations.'], ...
        fit.nAccepted,nAcceptedTarget,imc);
end

if any(accepted)
    fit.summary = struct();
    fit.summary.alphaMedian = median(aMC(accepted),'omitnan');
    fit.summary.alpha95 = prctile(aMC(accepted),[2.5 97.5]);
    fit.summary.gammaMedian = median(gMC(accepted),'omitnan');
    fit.summary.gamma95 = prctile(gMC(accepted),[2.5 97.5]);
    fit.summary.SSEMedian = median(sseMC(accepted),'omitnan');
else
    fit.summary = struct();
end

end

function SSE = objective_log_median(x,obs,vdMethod,wsMethod,unpack)
[alpha,gamma] = unpack(x);

wdModel = nan(height(obs),1);
caseKeys = unique(obs.case_key,'stable');

for ic = 1:numel(caseKeys)
    idx = obs.case_key == caseKeys(ic);
    ii = find(idx,1,'first');

    out = deposition_velocity( ...
        obs.ustar_mps(ii),obs.T_K(ii),obs.p_Pa(ii),obs.L_m(ii), ...
        obs.zr_m(ii),obs.Dp_m(idx), ...
        'vdMethod',vdMethod,'wsMethod',wsMethod, ...
        'z0',obs.z0_m(ii),'rhop',obs.rhop(ii), ...
        'alpha',alpha,'gamma',gamma);

    wdModel(idx) = out.wd(:);
end

good = isfinite(wdModel) & wdModel>0 & ...
       isfinite(obs.wd_obs_mps) & obs.wd_obs_mps>0;

if ~any(good)
    SSE = realmax;
    return
end

r = log(wdModel(good)) - log(obs.wd_obs_mps(good));
Dp_um = obs.Dp_m(good)*1e6;
[~,~,iu] = unique(round(Dp_um,6));

rMed = accumarray(iu,r,[],@median,NaN);
rMed = rMed(isfinite(rMed));

if isempty(rMed)
    SSE = realmax;
else
    SSE = sum(rMed.^2);
end
end

function y = bounded_forward(x,b)
s = 1./(1+exp(-x));
y = b(1) + (b(2)-b(1))*s;
end

function x = bounded_inverse(y,b)
f = (y-b(1))/(b(2)-b(1));
f = min(max(f,1e-10),1-1e-10);
x = log(f/(1-f));
end
