function out = deposition_velocity(ustar, T, p, L, z_r, Dp, varargin)
%DEPOSITION_VELOCITY Compute aerosol dry deposition velocity at z_r.
%
% out = deposition_velocity(ustar,T,p,L,z_r,Dp,Name,Value,...)
%
% Computes aerosol dry deposition velocity w_d [m s^-1] at reference height
% z_r [m] for particle diameters Dp [m] (vector). Settling velocity w_s is
% computed by SETTLING_VELOCITY (must be on the MATLAB path).
%
% INPUTS
%   ustar  friction velocity [m s^-1]
%   T      air temperature [K]
%   p      air pressure [Pa]
%   L      Monin-Obukhov length [m] (L<0 unstable, L=0 neutral, L>0 stable)
%   z_r    reference height [m]
%   Dp     particle diameter(s) [m] (vector)
%
% KEY OPTIONS (Name-Value)
%   'vdMethod'  deposition scheme (case-insensitive)
%      'MCzhang01'        mass-conserving Zhang01, rebound (bare)
%      'MCzhang01_water'  mass-conserving Zhang01, no rebound (sticky)
%      'zhang2001'        Zhang et al. (2001) electric analogue
%      'zhang2001R1'      Zhang et al. (2001) w/updated rebound
%      'pleim2022'        Pleim-type (mass-conserving form)
%      'pleim2022_water'  Pleim water variant
%      'petroffzhang2010' Petroff & Zhang (2010)
%
%   'alpha','gamma'  tuning parameters used by MCzhang01*:
%      Eb  = Sc.^(-gamma)
%      EIM = (St/(alpha+St)).^2
%     Defaults: gamma=0.54 (bare) or 0.50 (sticky); alpha=1.
%
% OUTPUT (struct)
%   out.wd, out.ws, out.ra, out.rb and diagnostic fields.

% FROM ZHANG FIT
% x = out.nonsticky.bootstrap.gamma;
% [prctile(x,25) median(x) prctile(x,75)]
% ans = 0.2932    0.3296    0.3516
%
% x = out.nonsticky.bootstrap.alpha;
% [prctile(x,25) median(x) prctile(x,75)]
% ans = 3.0204    3.8342    5.5694

% -----------------------------
% Parse inputs
% -----------------------------
pIn = inputParser;
pIn.FunctionName = 'deposition_velocity';

validateattributes(ustar, {'numeric'},{'scalar','positive'});
validateattributes(T,     {'numeric'},{'scalar','positive'});
validateattributes(p,     {'numeric'},{'scalar','positive'});
validateattributes(L,     {'numeric'},{'scalar'});          % can be neg/pos/0
validateattributes(z_r,   {'numeric'},{'scalar','positive'});
validateattributes(Dp,    {'numeric'},{'vector','positive'}); % [m]
Dp = Dp(:); % enforce column early

% Core options
addParameter(pIn,'vdMethod','MCzhang01',@(s)ischar(s)||isstring(s));
addParameter(pIn,'wsMethod','bb',       @(s)ischar(s)||isstring(s));

% Surface / similarity
addParameter(pIn,'kappa',0.4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'z0',1e-3, @(x)isnumeric(x)&&isscalar(x)&&x>0);  % roughness length [m]

% Vertical integration grid for r_a
addParameter(pIn,'zGrid',[],@(x)isnumeric(x)&&isvector(x)&&all(x>0));

% Air props
addParameter(pIn,'kb',1.38065e-23,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'Rd',287,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'S',120,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'mu0',1.812e-5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'Tp',288.15,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'dm',3.2e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'g',9.81,@(x)isnumeric(x)&&isscalar(x)&&x>0);

% Particle props
addParameter(pIn,'rhop',2500,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'ARp',1.46,@(x)isnumeric(x)&&isscalar(x)&&x>0);

% Cunningham slip correction constants
addParameter(pIn,'CcA',1.257,@(x)isnumeric(x)&&isscalar(x));
addParameter(pIn,'CcB',0.40, @(x)isnumeric(x)&&isscalar(x));
addParameter(pIn,'CcC',1.10, @(x)isnumeric(x)&&isscalar(x));

% Pleim 2022 updated R1 (Preb = A*(1-exp(-B*ws)))
addParameter(pIn,'PleimPrebA',0.95,@(x)isnumeric(x)&&isscalar(x));
addParameter(pIn,'PleimPrebB',1.0, @(x)isnumeric(x)&&isscalar(x));

% MCzhang01 tuning parameters (fit knobs) - initial conditions based on
% fitting to Zhang et al 2014 data
addParameter(pIn,'alpha',3.3,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(pIn,'gamma',0.35, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% If using the electric analogue, the parameters are
%if lower(strtrim(string(opt.vdMethod))) == "zhang2001"
%    addParameter(pIn,'alpha',2.2435,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
%    addParameter(pIn,'gamma',0.3454, @(x)isnumeric(x)&&isscalar(x)&&x>0);
%end

% Petroff & Zhang ground parameter
addParameter(pIn,'CITg',0.14,@(x)isnumeric(x)&&isscalar(x)&&x>=0);

% Return profiles?
addParameter(pIn,'returnProfiles',false,@(x)islogical(x)&&isscalar(x));

parse(pIn,varargin{:});
opt = pIn.Results;

% -----------------------------
% Normalize/validate method strings (robust to case + aliases)
% -----------------------------
vdMethodRaw = lower(strtrim(string(opt.vdMethod)));
wsMethodRaw = lower(strtrim(string(opt.wsMethod)));

% Common aliases / synonyms
switch vdMethodRaw
    case {"zhang2001", "zhang01"}
        vdMethodRaw = "zhang2001";
    case {"zhang2001r1", "zhang01r1"}
        vdMethodRaw = "zhang2001r1";
    case {"pleim2022","pleim"}
        vdMethodRaw = "pleim2022";
    case {"pleim2022water","pleim_water"}
        vdMethodRaw = "pleim2022_water";
    case {"petroffzhang2010","petroff&zhang2010","petroff_zhang_2010","pz2010"}
        vdMethodRaw = "petroffzhang2010";
    case {"mczhang01","mczhang","mczhang01_bare"}
        vdMethodRaw = "MCzhang01";
    case {"mczhang01_water","mczhang_water","mczhang01_sticky"}
        vdMethodRaw = "MCzhang01_water";
end

validVd = ["zhang2001_LUC13","zhang2001","zhang2001_LUC1","pleim2022", ...
    "pleim2022_water","petroffzhang2010","MCzhang01","zhang2001r1", ...
    "MCzhang01_water","MCzhang01rs"];
vdMethod = validatestring(vdMethodRaw, validVd, mfilename, 'vdMethod');

validWs = ["stokes","cliftgauvin","bb","huang","stokes_asph"];
% validWs = ["stokes","cliftgauvin","stokes_asph","prolate_ro","cliftgauvin_asp"];
wsMethod = validatestring(wsMethodRaw, validWs, mfilename, 'wsMethod');

% -----------------------------
% Build z-grid (for r_a integration)
% -----------------------------
% if isempty(opt.zGrid)
%     zmin = max(opt.z0*1.05, 1e-4);
%     % keep a bit of near-surface resolution without starting below zmin
%     zA0  = max(zmin, 0.001);
%     zA1  = min(1, z_r);
%     z1 = linspace(zA0, zA1, 200);
%
%     z2 = [];
%     if z_r > 1
%         z2 = linspace(1.0, z_r, 200);
%     end
%     z = unique([zmin, z1, z2]);
% else
%     z = opt.zGrid(:).';
%     if z(1) <= 0
%         error('zGrid must be strictly positive.');
%     end
%     if z(end) < z_r
%         error('zGrid must extend to at least z_r.');
%     end
% end
z = logspace(log10(opt.z0),log10(z_r),200);


% -----------------------------
% Air properties at (T,p)
% -----------------------------
rhoa   = p/(opt.Rd*T);                                      % kg/m^3
mu     = opt.mu0*(opt.Tp+opt.S)/(T+opt.S)*(T/opt.Tp)^(3/2); % Pa s
nu     = mu/rhoa;                                           % m^2/s
lambda = opt.kb*T/(pi*sqrt(2.0)*opt.dm^2*p);                % m
Cc_fun = @(lam,D) 1 + (2*lam./D).*(opt.CcA + opt.CcB*exp(-(opt.CcC*D)./(2*lam)));

% -----------------------------
% Settling velocity (ws) at (T,p) for each Dp
% -----------------------------
ws = settling_velocity(T, p, Dp, wsMethod, 'rhop', opt.rhop, 'ARp', opt.ARp);
ws = ws(:);
nDp = numel(Dp);

% -----------------------------
% Similarity functions phiM, phiW
% -----------------------------
phiM = ones(size(z));
phiW = ones(size(z));

if L < 0
    phiM = (1 - 16*z/L).^(-1/4);
    phiW = 1.25*(1 - 3*z/L).^(1/3);
elseif L > 0
    phiM = 1 + 5*z/L;
    phiW = 1.25*(1 + 0.2*z/L);
end


% -----------------------------
% Aerodynamic resistance r_a(z): integrate 1/Kp from z0 to z
%   Kp = Km / Sct
%   Km = kappa*u* z / phiM
%   Sct = sqrt(1 + ws^2 / (u*^2 * phiW^2))
% -----------------------------
beta = 1;
kappa = opt.kappa;
Km = kappa * ustar * z ./ phiM;     % [1 x nz]
Sct = sqrt(1 + (beta.^2 .* ws.^2) ./ (ustar^2 * phiW.^2)); % [nDp x nz]
Kp  = Km ./ Sct;                                % [nDp x nz]

ra_prof = NaN(nDp, numel(z));
for i = 1:nDp
    ra_prof(i,:) = cumtrapz(z, 1./Kp(i,:));
end
% Robust z_r index (avoid exact floating-point match failures)
[~, izr] = min(abs(z - z_r));
ra = ra_prof(:, izr);

% -----------------------------
% Quasi-laminar layer inputs
% -----------------------------
Cc = Cc_fun(lambda, Dp);
D_brown = opt.kb*T.*Cc ./ (3*pi*mu.*Dp);    % [m^2/s]
Sc = nu ./ D_brown;                         % Schmidt number
St = ws * ustar^2 / (opt.g * nu);           % Stokes number (your definition)


% -----------------------------
% r_b and w_d by method
% -----------------------------
wd = NaN(nDp,1);
rb = NaN(nDp,1);

switch string(vdMethod)

    case "MCzhang01" % mass conserving version of the zhang model
        alpha = opt.alpha;
        gamma = opt.gamma;
        Eb  = Sc.^(-gamma);
        EIM = (St./(alpha + St)).^2;
        Preb = opt.PleimPrebA .* (1 - exp(-opt.PleimPrebB .* ws));
        R1   = max(1 - Preb, 1e-6);
        rb = 1 ./ (3*ustar*(Eb + EIM).*R1);
        R  = ra + rb;
        wd = ws ./ (1 - exp(-ws.*R));

    % Mass conserving, zhang EIM, explicit rs
    case "MCzhang01rs"
        alpha = opt.alpha;
        gamma = opt.gamma;

        % Molecular diffusivity
        EB = Sc.^(-gamma);

        % Impaction Efficiency
        EIM = (St./(alpha + St)).^2;

        % Rebound Probability
        Preb = opt.PleimPrebA .* (1 - exp(-opt.PleimPrebB .* ws));

        rb = 1./(3*ustar*(EB+EIM));
        rs = Preb./(ustar.*(1-Preb)); % Preb = 1 -> rs = infty; Preb = 0 -> rs = 0

        R  = ra + rb;
        wd = ws ./ (1 - exp(-ws.*R) .* (1-ws.*rs));

    case "MCzhang01_water"
        alpha = opt.alpha;
        gamma = opt.gamma;
        Eb  = Sc.^(-gamma);
        EIM = (St./(alpha + St)).^2;
        R1   = ones(size(EIM));
        rb = 1 ./ (3*ustar*(Eb + EIM).*R1);
        R  = ra + rb;
        wd = ws ./ (1 - exp(-ws.*R));

    case "zhang2001" % This matches petroff & zhang 2010
        alpha = opt.alpha;
        gamma = opt.gamma;
        ra_prof = cumtrapz(z, 1./Km);
        % ra = ra_prof(izr);
        Eb  = Sc.^(-gamma); % 0.54 is for LUC 8
        EIM = (St./(alpha + St)).^2; % 50 is for LUC 8, 500 is me
        R1   = exp(-sqrt(St));
        rb  = 1 ./ (3*ustar*(Eb + EIM).*R1);
        wd  = ws + 1 ./ (ra + rb + ra.*rb.*ws);


    case "zhang2001r1" % This matches petroff & zhang 2010
        alpha = opt.alpha;
        gamma = opt.gamma;
        ra_prof = cumtrapz(z, 1./Km);
        % ra = ra_prof(izr);
        Eb  = Sc.^(-gamma); % 0.54 is for LUC 8
        EIM = (St./(alpha + St)).^2; % 50 is for LUC 8, 500 is me
        Preb = opt.PleimPrebA .* (1 - exp(-opt.PleimPrebB .* ws));
        R1   = max(1 - Preb, 1e-6);
        rb  = 1 ./ (3*ustar*(Eb + EIM).*R1);
        wd  = ws + 1 ./ (ra + rb + ra.*rb.*ws);

    case "zhang2001_LUC1"
        A = 2e-3; % m
        St = ws * ustar^2 / (opt.g*A);
        ra_prof = cumtrapz(z, 1./Km);
        % ra = ra_prof(izr);
        Eb  = Sc.^(-0.56); % exponent for LUC 1
        EIM = (St./(1 + St)).^2; % 1 is for LUC 1
        EIN = (1/2) * (Dp/A).^2; % A = 2 mm for SC1
        R1  = exp(-sqrt(St));
        rb  = 1 ./ (3*ustar*(Eb + EIM + EIN).*R1);
        wd  = ws + 1 ./ (ra + rb + ra.*rb.*ws);

    case "zhang2001_LUC13" % water
        ra_prof = cumtrapz(z, 1./Km);
        % ra = ra_prof(izr);
        Eb  = (nu./D_brown).^(-0.5); % exponent for LUC 13
        EIM = (St./(100 + St)).^2; % 1 is for LUC 13
        R1  = exp(-sqrt(St));
        rb  = 1 ./ (3*ustar*(Eb + EIM).*R1);
        wd  = ws + 1 ./ (ra + rb + ra.*rb.*ws);

    case "pleim2022"
        Eb  = 1/3*Sc.^(-2/3);
        EIM = 10.^(-3./St);
        % EIM = (St./(500 + St)).^2; % 50 is for LUC 8, 500 is me
        Preb = opt.PleimPrebA .* (1 - exp(-opt.PleimPrebB .* ws));
        R1   = max(1 - Preb, 1e-6);
        rb = 1 ./ (ustar*(Eb + EIM).*R1);
        R  = ra + rb;
        wd = ws ./ (1 - exp(-ws.*R));

    case "pleim2022_water"
        U10 = ustar/kappa * log(10/opt.z0);
        a1 = 8.46e-5; a2 = 1.63e-6; a3 = -3.35e-8; b1 = 3.354; b2 = -0.062;
        a = a1 + a2*T + a3*T;
        b = b1 + b2*T;
        fwc = a*(b+U10)^2;
        Eb  = (1-fwc) .* 1/3*Sc.^(-2/3) + fwc .* ustar/U10;
        EIM = 10.^(-3./St);
        R1   = ones(size(EIM)); % no rebound off water or sticky wood
        rb = 1 ./ (ustar*(Eb + EIM).*R1);
        R  = ra + rb;
        wd = ws ./ (1 - exp(-ws.*R));

    case "petroffzhang2010"
        tau_p_plus = (ws./opt.g) .* (ustar.^2 ./ nu);
        F = (Sc.^(1/3)) ./ 2.9;

        term1 = (1/6) * log( ((1+F).^2) ./ (1 - F + F.^2) );
        term2 = (1/sqrt(3)) * atan( (2*F - 1) ./ sqrt(3) );
        term3 = pi / (6*sqrt(3));
        Egb   = (Sc.^(-2/3) ./ 14.5) .* (term1 + term2 + term3).^(-1);

        CITg = opt.CITg;
        Egt  = zeros(size(tau_p_plus));
        mask = (tau_p_plus < 20);
        Egt(mask)  = 2.5e-3 * CITg .* (tau_p_plus(mask)).^2;
        Egt(~mask) = CITg;

        Eg = Egb + Egt;
        rb = 1 ./ (Eg .* ustar);

        wd = ws + 1 ./ (ra + rb);

    otherwise
        error('Unhandled vdMethod "%s".', vdMethod);
end


% -----------------------------
% Package output
% -----------------------------
out = struct();
out.wd     = wd;
out.ws     = ws;
out.ra     = ra;
out.rb     = rb;
out.Eb     = Eb;
out.EIM    = EIM;
out.St     = St;
out.inputs = struct('ustar',ustar,'T',T,'p',p,'L',L,'z_r',z_r,'Dp',Dp, ...
                    'vdMethod',char(vdMethod),'wsMethod',char(wsMethod),'z0',opt.z0);
out.air    = struct('rhoa',rhoa,'mu',mu,'nu',nu,'lambda',lambda);

if opt.returnProfiles
    out.z       = z(:);
    out.ra_prof = ra_prof;
    out.phiM    = phiM(:);
    out.phiW    = phiW(:);
end

end
