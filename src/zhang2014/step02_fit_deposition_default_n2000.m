function results = step02_fit_deposition_default_n2000(varargin)
%STEP02_FIT_DEPOSITION_DEFAULT_N2000 Run both default fits to 2000 accepts.
%
% This is the one-vote-per-diameter fit used by
% step02_fit_deposition_default. It fits the mass-conserving and electric-
% analogue deposition formulations separately with Clift & Gauvin settling.
% Both runs use the same seed so trialIndex identifies the same perturbed
% observation realization in the two output files.
%
% Optional Name/Value:
%   'nAccepted' default 2000
%   'maxTrials' default 20000
%   'rngSeed'   default 1
%   'verbose'   default true

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'nAccepted',2000,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
addParameter(p,'maxTrials',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
addParameter(p,'rngSeed',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'verbose',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;

if opt.maxTrials<opt.nAccepted
    error('maxTrials must be at least nAccepted.');
end

thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(thisDir));
resultsDir = fullfile(projectRoot,'data','zhang2014','results');
mcFile = fullfile(resultsDir, ...
    sprintf('zhang14_fit_default_cliftgauvin_n%d.mat',opt.nAccepted));
eaFile = fullfile(resultsDir, ...
    sprintf('zhang14_fit_default_cliftgauvin_zhang01_n%d.mat',opt.nAccepted));

common = {'nMC',opt.nAccepted, ...
    'nAcceptedTarget',opt.nAccepted, ...
    'maxTrials',opt.maxTrials, ...
    'rngSeed',opt.rngSeed, ...
    'verbose',opt.verbose};

fprintf('\nMass-conserving fit (%d accepted realizations)\n',opt.nAccepted)
results.massConserving = step02_fit_deposition_default( ...
    'MCzhang01','cliftgauvin',common{:},'outputFile',mcFile);

fprintf('\nElectric-analogue fit (%d accepted realizations)\n',opt.nAccepted)
results.electricAnalogue = step02_fit_deposition_default( ...
    'zhang2001','cliftgauvin',common{:},'outputFile',eaFile);

results.files = struct('massConserving',mcFile,'electricAnalogue',eaFile);
fprintf('\nCompleted both fits.\n')
fprintf('Mass-conserving: %s\n',mcFile)
fprintf('Electric analogue: %s\n',eaFile)
end
