function Zhang14 = step01_ingest_zhang2014(csvFile)
%STEP01_INGEST_ZHANG2014  Read the Zhang14 CSV and save standardized .mat
%
%   Zhang14 = step01_ingest_zhang2014()
%   Zhang14 = step01_ingest_zhang2014(csvFile)
%
% If csvFile is not provided, the function uses the hard-coded
% path to 'Zhang14.csv' in the project’s data folder:
%   <projectRoot>/data/zhang2014/raw/Zhang14.csv
%
% The output .mat file is written as:
%   <projectRoot>/data/zhang2014/derived/Zhang14.mat

    % Determine project root based on this file’s location
    thisFile = mfilename('fullpath');
    projectRoot = fileparts(fileparts(fileparts(thisFile)));

    % Default CSV location (hard-coded)
    if nargin < 1 || isempty(csvFile)
        rawDir = fullfile(projectRoot,'data','zhang2014','raw');
        csvFile = fullfile(rawDir,'Zhang14.csv');
    end

    % Validate CSV
    if ~isfile(csvFile)
        error('step01_ingest_zhang2014:CSVNotFound', ...
              'Zhang14 CSV not found: %s', csvFile);
    end

    % Construct output MAT path
    derivedDir = fullfile(projectRoot,'data','zhang2014','derived');
    if ~exist(derivedDir,'dir')
        mkdir(derivedDir);
    end
    [~,name,~] = fileparts(csvFile);
    matFile = fullfile(derivedDir,[name '.mat']);

    % Call your CSV-to-MAT converter
    fprintf('Ingesting Zhang14 CSV from:\n  %s\n', csvFile);
    Zhang14 = Zhang_csv_to_mat(csvFile, matFile);

    fprintf('Saved derived .mat to:\n  %s\n', matFile);

end


function Zhang14 = Zhang_csv_to_mat(csvFile, matFile)
%ZHANG_CSV_TO_MAT  Parse Zhang14.csv into a single long table (standard units).
%
% USAGE
%   Zhang14 = Zhang_csv_to_mat();                            % defaults
%   Zhang14 = Zhang_csv_to_mat("Zhang14.csv");               % writes Zhang14.mat next to CSV
%   Zhang14 = Zhang_csv_to_mat("Zhang14.csv","Zhang14.mat"); % custom output
%
% OUTPUT (saved to MAT as variable 'Zhang14'):
%   Zhang14.meta : struct with provenance and units
%   Zhang14.T    : table (long format) with one row per (surface, u*, diameter)
%
% Units saved:
%   lengths: m
%   velocities: m/s
%   Dp: m
%
% Assumptions:
%   - wd and uncertainty in the CSV are in mm/s (per user)
%   - Dp column is in um (as labeled)
%
% Missing entries are preserved as NaN.

    if nargin < 1 || isempty(csvFile)
        csvFile = "Zhang14.csv";
    end
    if nargin < 2 || isempty(matFile)
        [p,n,~] = fileparts(csvFile);
        if isempty(p), p = "."; end
        matFile = fullfile(p, n + ".mat");
    end

    raw = readcell(csvFile, "Delimiter", ",");

    % Strip UTF-8 BOM if present in first cell
    if ischar(raw{1,1})
        raw{1,1} = regexprep(raw{1,1}, '^\xEF\xBB\xBF', '');
    end

    surfaceNames = { ...
        'Wood Plane', ...
        'Water Surface', ...
        'Sand Surface', ...
        'Sandy Loam Surface', ...
        'Gobi Surface' };

    toNum = @(x) local_toDouble(x);

    % ---- Find surface header rows ----
    headerRows = [];
    headerIdx  = [];
    for r = 1:size(raw,1)
        v = raw{r,1};
        if ischar(v) || isstring(v)
            vv = strtrim(string(v));
            hit = find(strcmp(vv, string(surfaceNames)), 1);
            if ~isempty(hit)
                headerRows(end+1,1) = r; %#ok<AGROW>
                headerIdx(end+1,1)  = hit; %#ok<AGROW>
            end
        end
    end
    if isempty(headerRows)
        error("No surface blocks found (expected surface headers in column 1).");
    end

    [headerRows, order] = sort(headerRows);
    headerIdx = headerIdx(order);

    % ---- First pass: count rows to preallocate ----
    % Each surface has nU cases (2–3) and nDp diameters (expected 9).
    % Total rows = sum_s (nU_s * nDp_s).
    totalRows = 0;
    blockInfo = struct("s",{}, "r1",{}, "rEnd",{}, "nU",{}, "nDp",{});

    for b = 1:numel(headerRows)
        r0 = headerRows(b);
        r1 = r0 + 1;
        if b < numel(headerRows), rEnd = headerRows(b+1) - 1; else, rEnd = size(raw,1); end

        % ustar row is r1+2, cols 2:4
        ustar = [toNum(raw{r1+2,2}), toNum(raw{r1+2,3}), toNum(raw{r1+2,4})];
        nU = sum(~isnan(ustar));

        % find Diameter header row
        diamHeader = [];
        for r = r1:rEnd
            if (ischar(raw{r,1}) || isstring(raw{r,1})) && contains(string(raw{r,1}), "Diameter", "IgnoreCase", true)
                diamHeader = r; break;
            end
        end
        if isempty(diamHeader)
            error("Could not find 'Diameter (um)' row in block %d.", b);
        end

        % count numeric diameter rows
        nDp = 0;
        for r = (diamHeader+1):rEnd
            dp = toNum(raw{r,1});
            if ~isnan(dp), nDp = nDp + 1; end
        end

        totalRows = totalRows + nU * nDp;

        blockInfo(end+1,1) = struct("s",headerIdx(b),"r1",r1,"rEnd",rEnd,"nU",nU,"nDp",nDp); %#ok<AGROW>
    end

    % ---- Preallocate columns ----
    surface     = strings(totalRows,1);
    surface_id  = zeros(totalRows,1,'uint8');
    case_id     = zeros(totalRows,1,'uint8');

    ustar_mps = nan(totalRows,1);
    z0_m      = nan(totalRows,1);
    zd_m      = nan(totalRows,1);
    zm_m      = nan(totalRows,1);

    Dp_m    = nan(totalRows,1);
    wd_mps  = nan(totalRows,1);
    sig_mps = nan(totalRows,1);

    % ---- Fill ----
    k = 0;
    for b = 1:numel(blockInfo)
        sIdx = blockInfo(b).s;
        r1   = blockInfo(b).r1;
        rEnd = blockInfo(b).rEnd;

        sName = surfaceNames{sIdx};

        % meta rows: r1+0..3
        zd    = [toNum(raw{r1+0,2}), toNum(raw{r1+0,3}), toNum(raw{r1+0,4})];
        zm    = [toNum(raw{r1+1,2}), toNum(raw{r1+1,3}), toNum(raw{r1+1,4})];
        ustar = [toNum(raw{r1+2,2}), toNum(raw{r1+2,3}), toNum(raw{r1+2,4})];
        z0    = [toNum(raw{r1+3,2}), toNum(raw{r1+3,3}), toNum(raw{r1+3,4})];

        nU = sum(~isnan(ustar));
        if nU == 0
            warning("Surface '%s' has no u* values; skipping.", sName);
            continue;
        end

        zd    = zd(1:nU);
        zm    = zm(1:nU);
        ustar = ustar(1:nU);
        z0    = z0(1:nU);

        % find Diameter header row
        diamHeader = [];
        for r = r1:rEnd
            if (ischar(raw{r,1}) || isstring(raw{r,1})) && contains(string(raw{r,1}), "Diameter", "IgnoreCase", true)
                diamHeader = r; break;
            end
        end
        dataStart = diamHeader + 1;

        % read diameter list and associated wd/sig in mm/s
        dp_um_list = [];
        wd_mmps_mat  = [];  % [nDp x nU]
        sig_mmps_mat = [];  % [nDp x nU]

        for r = dataStart:rEnd
            dp_um = toNum(raw{r,1});
            if isnan(dp_um), continue; end

            dp_um_list(end+1,1) = dp_um; %#ok<AGROW>

            meas = nan(1,nU);
            sig  = nan(1,nU);
            for j = 1:nU
                meas(j) = toNum(raw{r, 1+j});  % cols 2..4
                sig(j)  = toNum(raw{r, 4+j});  % cols 5..7
            end
            wd_mmps_mat(end+1,1:nU)  = meas; %#ok<AGROW>
            sig_mmps_mat(end+1,1:nU) = sig;  %#ok<AGROW>
        end

        % Convert once
        dp_m_list = dp_um_list * 1e-6;     % um -> m
        wd_mps_mat  = wd_mmps_mat  * 1e-3; % mm/s -> m/s
        sig_mps_mat = sig_mmps_mat * 1e-3;

        nDp = numel(dp_m_list);

        % expand into long rows
        for c = 1:nU
            for i = 1:nDp
                k = k + 1;

                surface(k)    = string(sName);
                surface_id(k) = uint8(sIdx);
                case_id(k)    = uint8(c);

                ustar_mps(k) = ustar(c);
                z0_m(k)      = z0(c);
                zd_m(k)      = zd(c);
                zm_m(k)      = zm(c);

                Dp_m(k)    = dp_m_list(i);
                wd_mps(k)  = wd_mps_mat(i,c);
                sig_mps(k) = sig_mps_mat(i,c);
            end
        end
    end

    % trim if any blocks skipped
    if k < totalRows
        surface     = surface(1:k);
        surface_id  = surface_id(1:k);
        case_id     = case_id(1:k);
        ustar_mps   = ustar_mps(1:k);
        z0_m        = z0_m(1:k);
        zd_m        = zd_m(1:k);
        zm_m        = zm_m(1:k);
        Dp_m        = Dp_m(1:k);
        wd_mps      = wd_mps(1:k);
        sig_mps     = sig_mps(1:k);
    end

    % Build table
    T = table(surface, surface_id, case_id, ...
              ustar_mps, z0_m, zd_m, zm_m, ...
              Dp_m, wd_mps, sig_mps);

    % Add a convenience key for grouping (optional but useful)
    % One case = unique (surface_id, case_id)
    T.case_key = strcat("s", string(T.surface_id), "_c", string(T.case_id));

    % Package output
    Zhang14 = struct();
    Zhang14.meta = struct();
    Zhang14.meta.source_csv = char(csvFile);
    Zhang14.meta.assumed_input_units = struct("Dp_um","um","wd_mmps","mm/s","sig_mmps","mm/s");
    Zhang14.meta.saved_units = struct("length","m","velocity","m/s","Dp_m","m");
    Zhang14.meta.columns = T.Properties.VariableNames;

    Zhang14.T = T;

    save(matFile, "Zhang14", "-v7.3");
    fprintf("Saved %d rows to: %s\n", height(T), matFile);
end

function y = local_toDouble(x)
    if isempty(x), y = NaN; return; end
    if isnumeric(x), y = double(x); return; end
    if islogical(x), y = double(x); return; end
    if isstring(x) || ischar(x)
        s = strtrim(string(x));
        if s == "" || lower(s) == "nan"
            y = NaN;
            return;
        end
        y = str2double(s);
        if isnan(y), y = NaN; end
        return;
    end
    y = NaN;
end
