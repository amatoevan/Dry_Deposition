function ws = settling_velocity(T, p, Dp, method, varargin)
%SETTLING_VELOCITY  Terminal settling velocity for aerosol/dust particles.
%
%   ws = settling_velocity(T, p, Dp, method, 'Name', value, ...)
%
% Inputs
%   T       Temperature [K] (scalar or array)
%   p       Pressure    [Pa] (same size as T, or scalar)
%   Dp      Particle diameter(s) [m] (scalar or vector)
%   method  One of:
%           'stokes'        : Stokes settling with Cunningham slip correction
%           'cliftgauvin'   : Clift & Gauvin drag correlation (solve for Re)
%           'stokes_asph'   : Stokes * asphericity factor gamma_shape
%           'prolate_ro'    : randomly oriented prolate spheroid (solve for Re)
%
% Name-value options (defaults match your snippet where applicable)
%   'rhop'          Particle density [kg/m^3] (default 2500)
%   'g'             Gravity [m/s^2]           (default 9.81)
%   'Rd'            Gas constant dry air      (default 287)
%   'kb'            Boltzmann constant        (default 1.38065e-23)
%   'dm'            Air molecule diameter [m] (default 3.2e-10)
%   'mu0'           Sutherland mu0 [Pa s]     (default 1.812e-5)
%   'Tp'            Sutherland Tref [K]       (default 288.15)
%   'S'             Sutherland S [K]          (default 120)
%
%   Slip correction (same as your Cc_fun):
%   'CcA'           (default 1.257)
%   'CcB'           (default 0.40)
%   'CcC'           (default 1.10)
%
%   Asphericity factor options (only used for 'stokes_asph')
%   'gamma_shape'   Scalar multiplicative factor gamma (default 1)
%
%   Prolate options (only used for 'prolate_ro')
%   'ARp'           Prolate aspect ratio L/W (default 1.5)
%
%   Root finding options (used for Re-solves)
%   'ReBracket'     [Re_min Re_max] bracket (default [1e-12 1e3])
%   'ReTolX'        fzero TolX (default 1e-10)
%   'ReMaxIter'     fzero MaxIter (default 200)
%
% Output
%   ws  Settling velocity [m/s], size [numel(Dp) x numel(T)] after reshaping T(:).
%       If T and p are scalars, ws is [numel(Dp) x 1].
%
% Notes
%   - Includes buoyancy via (rhop - rhoa), consistent with your Stokes block.
%   - Air viscosity and mean free path follow your snippet exactly.
%   - For 'cliftgauvin' and 'prolate_ro', Re is solved per (Dp,T,p) point.
%
% Example (your ISA profile arrays T(z), p(z) and Dp vector):
%   wsS  = settling_velocity(T, p, Dp, 'stokes',        'rhop', 2500);
%   wsCG = settling_velocity(T, p, Dp, 'cliftgauvin',   'rhop', 2500);
%   wsA  = settling_velocity(T, p, Dp, 'stokes_asph',   'gamma_shape', 0.8);
%   wsP  = settling_velocity(T, p, Dp, 'prolate_ro',    'ARp', 1.5);
%
% Evan / UC Dust workflow: outputs are arranged ws(iDp, jState) like your loops.

% -----------------------------
% Parse inputs
% -----------------------------
method = lower(string(method));
validMethods = ["stokes","cliftgauvin","bb","huang","stokes_asph"];
if ~any(method == validMethods)
    error('Unknown method "%s". Choose: %s', method, strjoin(validMethods, ', '));
end

pIn = inputParser;
pIn.FunctionName = 'settling_velocity';

% Physical constants / air props (match your snippet defaults)
addParameter(pIn,'rhop',2500,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'g',9.81,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'Rd',287,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'kb',1.38065e-23,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'dm',3.2e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'mu0',1.812e-5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'Tp',288.15,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'S',120,@(x)isnumeric(x)&&isscalar(x)&&x>0);

% Cunningham slip correction params
addParameter(pIn,'CcA',1.257,@(x)isnumeric(x)&&isscalar(x));
addParameter(pIn,'CcB',0.40,@(x)isnumeric(x)&&isscalar(x));
addParameter(pIn,'CcC',1.10,@(x)isnumeric(x)&&isscalar(x));

% Prolate
addParameter(pIn,'ARp',1.7,@(x)isnumeric(x)&&isscalar(x)&&x>0);

% Root solve options
addParameter(pIn,'ReBracket',[1e-12 1e3],@(x)isnumeric(x)&&numel(x)==2&&all(x>0)&&x(1)<x(2));
addParameter(pIn,'ReTolX',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'ReMaxIter',200,@(x)isnumeric(x)&&isscalar(x)&&x>0);

parse(pIn,varargin{:});
opt = pIn.Results;

% -----------------------------
% Normalize shapes and validate
% -----------------------------
if ~isnumeric(T) || ~isnumeric(p) || ~isnumeric(Dp)
    error('T, p, and Dp must be numeric.');
end
if any(Dp(:) <= 0)
    error('All Dp must be > 0.');
end
Tvec = T(:);
if isscalar(p)
    pvec = repmat(p, size(Tvec));
else
    pvec = p(:);
    if numel(pvec) ~= numel(Tvec)
        error('If p is not scalar, p must have the same number of elements as T.');
    end
end
if any(Tvec <= 0) || any(pvec <= 0)
    error('All T and p must be > 0.');
end

Dpvec = Dp(:);
nDp = numel(Dpvec);
nS  = numel(Tvec);


% -----------------------------
% Air properties
% -----------------------------
rhoa = pvec ./ (opt.Rd .* Tvec); % kg/m^3

% - confirmed agreement with Fig. 3 in Petroff & Zhang 2010
mu = opt.mu0*(opt.Tp+opt.S)./(Tvec+opt.S).*(Tvec/opt.Tp).^(3/2); % Pa s
nu = mu./rhoa;
lambda = opt.kb*Tvec ./ (pi*sqrt(2.0)*opt.dm^2 .* pvec); % m

% From Carlos - similar
% nu = 1.45e-5;
% mu = nu*rhoa;
% lambda = (1.1e-3*sqrt(T)) ./ (p*1e-2);  % p is in pascal
Cc_fun = @(lam,D) 1 + (2*lam./D).*(opt.CcA + opt.CcB*exp(-(opt.CcC*D)./(2*lam)));


% Precompute repeated state arrays for speed/readability
% We want outputs ws(iDp, jState)
ws = NaN(nDp, nS);
Rey = NaN(nDp, nS);


% -----------------------------
% Particle geometry from Huang
%
% To relate this to the 3-axle ellipsoid:
% AR = S L I^(-2)
% S: shortest length of the particle
% L: longest length of the particle
% I: Intermediate length of the particle, [m]
%
% To relate to HWR from Huang 2020
% H = S
% W = I
% --> HWR = S/I
% --> f = HWR
%
% To relate to AR
% e = I/L 
% --> e = S / (I AR) = HWR/AR
% -----------------------------
HWR_ref = 0.4;
AR_ref = opt.ARp;
f = HWR_ref;
e = 1 / AR_ref;


% -----------------------------
% Compute by method
% -----------------------------
switch method
    case "stokes"
        for i = 1:nDp
            D = Dpvec(i);
            Cc = Cc_fun(lambda, D); % [nS x 1]
            ws(i,:) = (opt.g*(opt.rhop - rhoa)./(18*mu) .* Cc .* D.^2).';
        end

    case "stokes_asph"
        Fs = f*e.^1.3; 
        kS = 0.5*(Fs^(1/3)+Fs.^(-1/3));
        gamma = 1 / kS;
        for i = 1:nDp
            D = Dpvec(i);
            Cc = Cc_fun(lambda, D);
            ws_st = opt.g*(opt.rhop - rhoa)./(18*mu) .* Cc .* D.^2;
            ws(i,:) = (gamma * ws_st).';
        end

    case "cliftgauvin"
        bracket = opt.ReBracket;
        fzopt = optimset('TolX', opt.ReTolX, 'MaxIter', opt.ReMaxIter, 'Display','off');

        for i = 1:nDp
            D = Dpvec(i);
            for j = 1:nS
                Cc = Cc_fun(lambda(j), D);
                Cd_sph = @(Re) (24./Re).*(1+0.15*Re^0.687) + 0.42./(1+42500/(Re^1.16));
                fun = @(Re) (4*Cc*(opt.rhop-rhoa(j))*D^3*opt.g*rhoa(j))./(3*mu(j)^2) - ...
                    Cd_sph(Re) .* Re.^2;
                try
                    Re = fzero(fun, bracket, fzopt);
                    Rey(i,j) = Re;
                    ws(i,j) = mu(j)*Re./(D*rhoa(j));
                catch
                    ws(i,j) = NaN;
                end
            end
        end


    % case "cliftgauvin_asp"
    %     gamma = opt.gamma_shape;      % gamma = ws_asp / ws_sph (e.g., 0.85)
    %     chi   = 1/gamma;              % drag multiplier
    % 
    %     bracket = opt.ReBracket;
    %     fzopt = optimset('TolX', opt.ReTolX, 'MaxIter', opt.ReMaxIter, 'Display','off');
    % 
    %     for i = 1:nDp
    %         D = Dpvec(i);
    %         for j = 1:nS
    %             Cc = Cc_fun(lambda(j), D);
    % 
    %             Cd_sph = @(Re) (24./Re).*(1+0.15*Re^0.687) + 0.42./(1+42500/(Re^1.16));
    %             Cd_asp = @(Re) chi .* Cd_sph(Re);
    % 
    %             fun = @(Re) (8*Cc*(opt.rhop-rhoa(j))*D^3*opt.g*rhoa(j))./(24*mu(j)^2) - ...
    %                         Cd_asp(Re) .* Re.^2;
    % 
    %             try
    %                 Re = fzero(fun, bracket, fzopt);
    %                 ws(i,j) = 2*mu(j)*Re./(D*rhoa(j));   % NO gamma here
    %             catch
    %                 ws(i,j) = NaN;
    %             end
    %         end
    %     end

    % case "prolate_ro"
    % 
    %     % To relate this to the 3-axle ellipsoid:
    %     % AR = S L I^(-2)
    %     % S: shortest length of the particle
    %     % L: longest length of the particle
    %     % I: Intermediate length of the particle, [m]
    %     %
    %     % And the relation to the flatness (f) & elongation (e)
    %     % f = S/I
    %     % e = I/L
    %     %
    %     % For prolate spheroid:
    %     % I = S -> AR = L/I
    %     % --> e = 1/AR
    %     % --> f = 1;
    %     %
    %     % To relate to HWR from Huang 2020
    %     % H = S
    %     % W = I
    %     % --> HWR = S/I
    %     % --> f = HWR
    %     %
    %     % To relate to AR
    %     % e = I/L 
    %     % --> e = S / I AR = HWR/AR
    % 
    %     ARp = opt.ARp;
    %     bracket = opt.ReBracket;
    %     fzopt = optimset('TolX', opt.ReTolX, 'MaxIter', opt.ReMaxIter, 'Display','off');
    % 
    %     for i = 1:nDp
    %         D = Dpvec(i);
    %         for j = 1:nS
    % 
    %             % Mallios - random prolate spheroid (Eqs. 31--35)
    %             rr = opt.rhop./rhoa(j);
    %             a2 = 0.45 + 10*(30 + exp(2.5*log10(rr))).^(-1);
    %             b2 = 1 - 37*(exp(3*log10(rr))+100).^(-1);
    %             kN = 10^(a2.*(-log10(1./ARp)).^b2);
    %             kS = 0.5*ARp^(-0.433) + 0.5*ARp^(0.433);
    % 
    %             Cc = Cc_fun(lambda(j), D);
    % 
    %             fun = @(Re) (8*Cc*(opt.rhop-rhoa(j))*D^3*opt.g*rhoa(j))/(24*mu(j)^2) - ...
    %                 (12./Re.*ARp^(2/3).*kS.*(1+0.1984.*(Re*kN/kS).^(2/3)) + ...
    %                 0.46*kN ./(1+(2665*ARp^(2/3))./(Re * kN/kS))) .* Re.^2;
    % 
    %             try
    %                 Re = fzero(fun, bracket, fzopt);
    %                 Rey(i,j) = Re;
    %                 ws(i,j) = 2*mu(j)*Re./(D*rhoa(j));
    %             catch
    %                 ws(i,j) = NaN;
    %             end
    %         end
    %     end

    case "bb"
        bracket = opt.ReBracket;
        fzopt = optimset('TolX', opt.ReTolX, 'MaxIter', opt.ReMaxIter, 'Display','off');

        for i = 1:nDp
            D = Dpvec(i);
            for j = 1:nS
                % shape factors from Table 8 in B&B 2016
                rr = opt.rhop./rhoa(j);
                b2 = 1 - 37*(exp(3*log10(rr))+100).^(-1);
                a2 = 0.45 + 10*(exp(2.5*log10(rr)+30)).^(-1);
                Fs = f*e.^1.3; 
                Fn = f.^2*e;
                kN = 10.^(a2.*(-log10(Fn)).^b2);
                kS = 0.5*(Fs^(1/3)+Fs.^(-1/3));

                Cc = Cc_fun(lambda(j), D);

                Cd = @(Re) 24.*kS./Re .* (1+0.125*(Re.*kN./kS).^(2/3)) + ...
                    0.46*kN ./ (1+5330./(Re.*kN./kS));

                fun = @(Re) (4*opt.g*Cc*D^3*(opt.rhop-rhoa(j))*rhoa(j))/(3*mu(j)^2) - ...
                     Cd(Re) .* Re.^2;

                try
                    Re = fzero(fun, bracket, fzopt);
                    Rey(i,j) = Re;
                    ws(i,j) = mu(j)*Re./(D*rhoa(j));
                catch
                    ws(i,j) = NaN;
                end
            end
        end

        case "huang"
        bracket = opt.ReBracket;
        fzopt = optimset('TolX', opt.ReTolX, 'MaxIter', opt.ReMaxIter, 'Display','off');

        for i = 1:nDp
            D = Dpvec(i);
            for j = 1:nS
                % B&B Eqs 13, 24, 25

                Fs = f*e.^1.3; 
                kS = 0.5*(Fs^(1/3)+Fs.^(-1/3));

                Cc = Cc_fun(lambda(j), D);

                Cd = @(Re) 24.*kS./Re;

                fun = @(Re) (4*opt.g*Cc*D^3*(opt.rhop-rhoa(j))*rhoa(j))/(3*mu(j)^2) - ...
                     Cd(Re) .* Re.^2;

                try
                    Re = fzero(fun, bracket, fzopt);
                    Rey(i,j) = Re;
                    ws(i,j) = mu(j)*Re./(D*rhoa(j));
                catch
                    ws(i,j) = NaN;
                end
            end
        end

end
end

