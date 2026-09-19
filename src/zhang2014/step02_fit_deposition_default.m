function results = step02_fit_deposition_default(vdMethod,wsMethod,varargin)
%STEP02_FIT_DEPOSITION Fit Zhang14 alpha, gamma, and Cw using measurement-error MC.
%
% USAGE
%   results = step02_fit_deposition();
%   results = step02_fit_deposition('MCzhang01','bb');
%
% Default methods:
%   vdMethod = 'MCzhang01'
%   wsMethod = 'bb'
%
% Optional Name/Value:
%   'nMC'          default 500
%   'nAcceptedTarget' default []
%   'maxTrials'    default []
%   'alphaBounds'  default [1e-3 1e2]
%   'gammaBounds'  default [1e-3 1e1]
%   'rngSeed'      default 1
%   'verbose'      default true
%   'outputFile'   default '' (uses the historical default filename)
%
% This code applies the one-vote-per-diameter fitting routine
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projectRoot = fileparts(fileparts(fileparts(thisFile)));

% Local Zhang fitter and shared deposition/settling functions in ../.
addpath(fullfile(thisDir,'utils'));
addpath(fileparts(thisDir));

if nargin<1 || isempty(vdMethod)
    vdMethod = 'zhang2001'; % 'MCzhang01';
end
if nargin<2 || isempty(wsMethod)
    wsMethod = 'cliftgauvin';
end

vdMethod = char(string(vdMethod));
wsMethod = char(string(wsMethod));

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'nMC',500,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'nAcceptedTarget',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0));
addParameter(p,'maxTrials',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0&&mod(x,1)==0));
addParameter(p,'alphaBounds',[1e-3 1e2],@(x)isnumeric(x)&&numel(x)==2&&all(x>0));
addParameter(p,'gammaBounds',[1e-3 1e1],@(x)isnumeric(x)&&numel(x)==2&&all(x>0));
addParameter(p,'rngSeed',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'verbose',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'outputFile','',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opt = p.Results;

dataFile = fullfile(projectRoot,'data','zhang2014','derived','Zhang14.mat');
if ~isfile(dataFile)
    error('Expected dataset not found: %s',dataFile);
end

S = load(dataFile,'Zhang14');
T = S.Zhang14.T;

% Fit non-sticky mineral surfaces only.
isWood = contains(lower(string(T.surface)),'wood');
isWater = contains(lower(string(T.surface)),'water');
T = T(~(isWood|isWater),:);

obs = table();
obs.source = repmat("Zhang14",height(T),1);
obs.surface = string(T.surface);
obs.case_key = string(T.case_key);
obs.ustar_mps = T.ustar_mps;
obs.z0_m = T.z0_m;
obs.zr_m = T.zm_m; % displacement height is zero
obs.Dp_m = T.Dp_m;
obs.wd_obs_mps = T.wd_mps;
obs.sigma_mps = T.sig_mps;
obs.T_K = 298*ones(height(T),1);
obs.p_Pa = 101300*ones(height(T),1);
obs.L_m = zeros(height(T),1);
obs.rhop = 2200*ones(height(T),1);

good = isfinite(obs.wd_obs_mps) & obs.wd_obs_mps>0 & ...
       isfinite(obs.Dp_m) & obs.Dp_m>0;
obs = obs(good,:);

if opt.verbose
    fprintf('Zhang14 non-sticky observations: %d\n',height(obs));
    fprintf('Unique experimental cases: %d\n',numel(unique(obs.case_key)));
    fprintf('Fitting with vdMethod=%s | wsMethod=%s\n',vdMethod,wsMethod);
end

fit = fit_zhang2014_alpha_gamma_mc(obs,vdMethod,wsMethod, ...
    'nMC',opt.nMC, ...
    'nAcceptedTarget',opt.nAcceptedTarget, ...
    'maxTrials',opt.maxTrials, ...
    'alphaBounds',opt.alphaBounds, ...
    'gammaBounds',opt.gammaBounds, ...
    'rngSeed',opt.rngSeed, ...
    'verbose',opt.verbose);

results = struct();
results.obs = obs;
results.options = opt;
results.vdMethod = vdMethod;
results.wsMethod = wsMethod;
results.nonsticky = fit;

resultsDir = fullfile(projectRoot,'data','zhang2014','results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

if strlength(string(opt.outputFile))==0
    tag = 'default_cliftgauvin_zhang01';
    outFile = fullfile(resultsDir,['zhang14_fit_' tag '.mat']);
else
    outFile = char(string(opt.outputFile));
end
save(outFile,'results','-v7.3');

fprintf('\nNominal fit:\n');
fprintf('  alpha = %.6g\n',fit.nominal.alpha);
fprintf('  gamma = %.6g\n',fit.nominal.gamma);
fprintf('  SSE   = %.6g\n',fit.nominal.SSE);

if isfield(fit,'summary') && ~isempty(fieldnames(fit.summary))
    fprintf('\nAccepted MC fits: %d of %d attempts\n', ...
        fit.nAccepted,fit.nAttempted);
    fprintf('  alpha median = %.6g; 95%% range = %.6g to %.6g\n', ...
        fit.summary.alphaMedian,fit.summary.alpha95(1),fit.summary.alpha95(2));
    fprintf('  gamma median = %.6g; 95%% range = %.6g to %.6g\n', ...
        fit.summary.gammaMedian,fit.summary.gamma95(1),fit.summary.gamma95(2));
    fprintf('  median SSE   = %.6g\n',fit.summary.SSEMedian);
end

fprintf('Saved fit results to: %s\n',outFile);

end
