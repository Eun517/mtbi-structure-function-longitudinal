%% FNC analysis for Brain Communications (rs-fMRI) 
% This script reproduces all functional network connectivity (FNC) analyses:
%   1) Compute subject-level FNC (Fisher Z)
%   2) Edge-wise LME: FNC ~ Age + Sex + Group*Time + (1|Subject)
%   3) FDR-BH Correction
%   4) Post-hoc t-tests (Pat_vs_Con, T0_vs_T1)
%   5) Cognition LME (Patients only): Y ~ FNC*Time + (1|Subject)
%   6) Pearson Correlations (FNC & Cognition)
%   7) Joint LME: Y ~ FNC*Group*Time + Age + Sex + (1|Subject)
%
% IMPORTANT
% - This script assumes that time-series matrices and metadata files are
%   prepared separately and stored under ./inputs (see "Project structure").
%% -------------------------------------------------------------------------
% Project structure (recommended)
%   repo_root/
%     matlab/
%       run_fnc_analysis.m           (this script)
%     inputs/                        (NOT committed; ignored by .gitignore)
%       age_patients.txt
%       age_controls.txt
%       sex_patients.txt
%       sex_controls.txt
%       NetInfo.txt
%       data_T0.mat                  (contains data_T0)
%       data_T1.mat                  (contains data_T1)
%       cognitive_patients.mat       (optional; contains cognitive_patients)
%       cognitive_controls.mat       (optional; contains cognitive_controls)
%     outputs/                       (optional; can be ignored)
%
% If you use a different layout, update "cfg" paths below.

%% Housekeeping
clearvars; clc;

%% -------------------------------------------------------------------------
% 0) Configuration
cfg = struct();

% Root directory of the repository (assumes this script is in repo_root/matlab/)
cfg.repo_root   = fileparts(fileparts(mfilename('fullpath')));
cfg.input_dir   = fullfile(cfg.repo_root, 'inputs');
cfg.output_dir  = fullfile(cfg.repo_root, 'outputs');

% Input files
cfg.age_pat_file   = fullfile(cfg.input_dir, 'age_patients.txt');
cfg.age_con_file   = fullfile(cfg.input_dir, 'age_controls.txt');
cfg.sex_pat_file   = fullfile(cfg.input_dir, 'sex_patients.txt');
cfg.sex_con_file   = fullfile(cfg.input_dir, 'sex_controls.txt');
cfg.netinfo_file   = fullfile(cfg.input_dir, 'NetInfo.txt'); % must contain network name
cfg.data_T0_file   = fullfile(cfg.input_dir, 'data_T0.mat');  % must contain "data_T0"
cfg.data_T1_file   = fullfile(cfg.input_dir, 'data_T1.mat');  % must contain "data_T1"

% Optional cognition files (set cfg.run_cognition = false to skip)
cfg.run_cognition      = true;
cfg.cog_pat_file       = fullfile(cfg.input_dir, 'cognitive_patients.mat'); % "cognitive_patients"
cfg.cog_con_file       = fullfile(cfg.input_dir, 'cognitive_controls.mat'); % "cognitive_controls"
cfg.cog_missing_code   = 999;
cfg.cog_score_idx      = [4 6]; % [WCST non-perseverative error, backward digit span]
cfg.cog_score_names    = {'WCST_nonpersev', 'DigitSpan_backward'};

% Statistics
cfg.alpha_uncorrected  = 0.05;
cfg.alpha_fdr          = 0.05;

% Notes:
% - This script requires fdr_bh.m (Benjamini-Hochberg FDR).
%   If you do not have it, add it to your MATLAB path (e.g., ./matlab/utils).
% - Sex coding is assumed to be 0/1 as used in your dataset.
%   Update comments and/or recode if needed.

%% -------------------------------------------------------------------------
% 1) Load demographics
fprintf('>> Loading data...\n');
age.con = load(cfg.age_con_file); age.pat = load(cfg.age_pat_file);
sex.con = load(cfg.sex_con_file); sex.pat = load(cfg.sex_pat_file);
netNames = importdata(cfg.netinfo_file);
nCon = numel(age.con); nPat = numel(age.pat);


%% -------------------------------------------------------------------------
% 2) Load fMRI time series (time x network x subject)
% Expected variables:
%   data_T0: [T x N x (nCon+nPat)]
%   data_T1: [T x N x (nCon+nPat)]

S0 = load(cfg.data_T0_file, 'data_T0');
S1 = load(cfg.data_T1_file, 'data_T1');

%% -------------------------------------------------------------------------
% 3) Compute FNC (network correlation) and Fisher Z
% Output:
%   zval.con: [nEdges x nCon x 2] 
%   zval.pat: [nEdges x nPat x 2]
%   nEdges = N × (N − 1) / 2
N = numel(netNames);      % number of networks
mask = triu(true(N), 1); % 
nEdges = nnz(mask);

times = 1:2;
zval  = struct();

% Preallocate for reproducibility and speed
zval.con = nan(nEdges, nCon, 2);
zval.pat = nan(nEdges, nPat, 2);


for t = times
    if t == 1
        S = load(cfg.data_T0_file, 'data_T0'); cur_data = S.data_T0;
    else
        S = load(cfg.data_T1_file, 'data_T1'); cur_data = S.data_T1;
    end
    
    % Controls
    for s = 1:nCon
        R = corr(cur_data(:,:,s));
        R = min(max(R, -0.999), 0.999); % avoid Inf
        zval.con(:,s,t) = atanh(R(mask));
    end
    % Patients
    for s = 1:nPat
        R = corr(cur_data(:,:,nCon+s));
        R = min(max(R, -0.999), 0.999);
        zval.pat(:,s,t) = atanh(R(mask));
    end
    clear S cur_data;
end


%% -------------------------------------------------------------------------
% 4) Edge-wise LME: FNC ~ Age + Sex + Group*Time + (1|Subject)
% Build long-format vectors per edge:
%   controls: [nCon x 2] -> stack to [2*nCon x 1]
%   patients: [nPat x 2] -> stack to [2*nPat x 1]
G = nan(nEdges,1);  % p(Group)
T = nan(nEdges,1);  % p(Time)
GT = nan(nEdges,1); % p(Group*Time) interaction (optional but useful)

Group = categorical([repmat({'Control'}, nCon*2, 1); repmat({'Patient'}, nPat*2, 1)]);
Time  = categorical([repmat({'T0'}, nCon, 1); repmat({'T1'}, nCon, 1); ...
    repmat({'T0'}, nPat, 1); repmat({'T1'}, nPat, 1)]);

Subject = [repmat((1:nCon)', 2, 1); repmat((nCon+1:nCon+nPat)', 2, 1)];

Age = [age.con(:); age.con(:); age.pat(:); age.pat(:)];
Sex = [sex.con(:); sex.con(:); sex.pat(:); sex.pat(:)];


for iEdge = 1:nEdges
    z_con = squeeze(zval.con(iEdge,:,:)); % nCon x 2
    z_pat = squeeze(zval.pat(iEdge,:,:)); % nPat x 2

    y = [z_con(:); z_pat(:)]; % (2*nCon + 2*nPat) x 1

    
    tbl = table(y, Age, Sex, Group, Time, categorical(Subject), ...
        'VariableNames', {'FNC','Age','Sex','Group','Time','Subject'});

    lme = fitlme(tbl, 'FNC ~ Age + Sex + Group*Time + (1|Subject)');

    coef = lme.Coefficients;

    % Robust extraction by coefficient name (avoids dependence on row order)
    G(iEdge)  = coef.pValue(strcmp(coef.Name, 'Group_Patient'));
    T(iEdge)  = coef.pValue(strcmp(coef.Name, 'Time_T1'));
    if any(strcmp(coef.Name, 'Group_Patient:Time_T1'))
        GT(iEdge) = coef.pValue(strcmp(coef.Name, 'Group_Patient:Time_T1'));
    end
end

% FDR correction (BH)
[G_q, G_h] = fdr_bh(G,  cfg.alpha_fdr, 'pdep', 'yes');
[T_q, T_h] = fdr_bh(T,  cfg.alpha_fdr, 'pdep', 'yes');
[GT_q, GT_h] = fdr_bh(GT, cfg.alpha_fdr, 'pdep', 'yes');

%% -------------------------------------------------------------------------
% 5) Edge annotation (node names / edge list)
% NetInfo.txt should contain one node/network name per line (Dim = N nodes)
% Upper triangle indices (must match FNC edge ordering)
[iNet, jNet] = find(mask);

% Build edge table (same order as zval.*(:, :, t))
EdgeTable = table( (1:nEdges)', iNet, jNet, ...
                   netNames(iNet), netNames(jNet), ...
    'VariableNames', {'EdgeID','Node1_idx','Node2_idx','Node1_name','Node2_name'});


%% -------------------------------------------------------------------------
% 6) Optional post-hoc tests for edges passing an uncorrected screen
% (Use only if you report these post-hoc contrasts in the manuscript.)
sig_edges_unc = find((G < cfg.alpha_uncorrected) | (T < cfg.alpha_uncorrected) | (GT < cfg.alpha_uncorrected));


sig_edges_unc_G = find(G < cfg.alpha_uncorrected); 
% Uncorrected Group-significant edges
% Saved for downstream use in brain_comm_DTI_metrics_clean.m


posthoc = struct();
posthoc.edge_id = sig_edges_unc;
posthoc.tests = {'Pat_vs_Con_T0', 'Pat_vs_Con_T1', 'Pat_T0_vs_T1', 'Con_T0_vs_T1'};
posthoc.p = nan(numel(sig_edges_unc), 4);
posthoc.t = nan(numel(sig_edges_unc), 4);

for k = 1:numel(sig_edges_unc)
    e = sig_edges_unc(k);
    z_con = squeeze(zval.con(e,:,:)); % nCon x 2
    z_pat = squeeze(zval.pat(e,:,:)); % nPat x 2

    pat_T0 = z_pat(:,1); pat_T1 = z_pat(:,2);
    con_T0 = z_con(:,1); con_T1 = z_con(:,2);

    [~, posthoc.p(k,1), ~, st] = ttest2(pat_T0, con_T0); posthoc.t(k,1) = st.tstat;
    [~, posthoc.p(k,2), ~, st] = ttest2(pat_T1, con_T1); posthoc.t(k,2) = st.tstat;
    [~, posthoc.p(k,3), ~, st] = ttest(pat_T0, pat_T1);  posthoc.t(k,3) = st.tstat;
    [~, posthoc.p(k,4), ~, st] = ttest(con_T0, con_T1);  posthoc.t(k,4) = st.tstat;
end

%% -------------------------------------------------------------------------
% 7) Optional cognition models (patients only)
% Cognitive matrices must be:
%   cognitive_patients: [nPat x nScoresTotal x 2]  (3rd dim: 1=T0, 2=T1)
%   cognitive_controls: [nCon x nScoresTotal x 2]  (loaded only for QC / future use)
%
% Patients-only LME (edge-wise):
%   Y ~ FNC * Time + (1|Subject)
%
% Requires in workspace:
%   zval.pat      : [nEdges x nPat x 2]
%   EdgeTable     : optional (for node labels; must match zval edge ordering)
%   sig_edges_unc : candidate edges (e.g., uncorrected significant from main FNC LME)
%   G, T, GT      : main FNC LME p-values (length = nEdges)

if cfg.run_cognition

    % ---------- Load cognitive data ----------
    assert(isfile(cfg.cog_pat_file), 'Missing file: %s', cfg.cog_pat_file);
    Cpat = load(cfg.cog_pat_file);
    assert(isfield(Cpat,'cognitive_patients'), 'cog_pat_file must contain "cognitive_patients".');
    cognitive_patients = Cpat.cognitive_patients;

    % Controls are optional (QC / future use)
    if isfield(cfg,'cog_con_file') && ~isempty(cfg.cog_con_file) && isfile(cfg.cog_con_file)
        Ccon = load(cfg.cog_con_file);
        if isfield(Ccon,'cognitive_controls')
            cognitive_controls = Ccon.cognitive_controls; %#ok<NASGU>
        end
    end

    % ---------- Sanity checks ----------
    nPat = size(cognitive_patients,1);
    assert(nPat == size(zval.pat,2), 'Mismatch: cognitive_patients rows vs zval.pat subjects.');

    % ---------- Select edges for cognition ----------
    edgesAny = sig_edges_unc(:);
    edgesTimeOnly = find( ~(G < cfg.alpha_uncorrected) & (T < cfg.alpha_uncorrected) & ~(GT < cfg.alpha_uncorrected) );
    edgesForCognition = setdiff(edgesAny, edgesTimeOnly);

    nEdgesUse = numel(edgesForCognition);
    if nEdgesUse == 0
        warning('No edges selected for cognition models. Skipping cognition LME.');
        return;
    end

    % ---------- Settings ----------
    scoreIdx     = cfg.cog_score_idx;      % e.g., [4 6]
    scoreNames   = cfg.cog_score_names;    % e.g., {'WCST_nonpersev','DigitSpan_backward'}
    missing_code = cfg.cog_missing_code;   % e.g., 999
    alpha_fdr    = cfg.alpha_fdr;
    nScores      = numel(scoreIdx);

    % ---------- Preallocate ----------
    beta_FNC  = nan(nEdgesUse, nScores);
    beta_Time = nan(nEdgesUse, nScores);
    beta_Int  = nan(nEdgesUse, nScores);

    p_FNC  = nan(nEdgesUse, nScores);
    p_Time = nan(nEdgesUse, nScores);
    p_Int  = nan(nEdgesUse, nScores);

    % ---------- Edge-wise cognition LME (patients only) ----------
    for s = 1:nScores
        cogIdx = scoreIdx(s);

        y_T0 = cognitive_patients(:, cogIdx, 1);
        y_T1 = cognitive_patients(:, cogIdx, 2);

        for i = 1:nEdgesUse
            edgeID = edgesForCognition(i);

            fnc_T0 = squeeze(zval.pat(edgeID, :, 1))';   % [nPat x 1]
            fnc_T1 = squeeze(zval.pat(edgeID, :, 2))';   % [nPat x 1]

            % Long format
            Y    = [y_T0; y_T1];
            FNC  = [fnc_T0; fnc_T1];
            Time = [zeros(nPat,1); ones(nPat,1)];        % 0=T0, 1=T1
            Subj = [(1:nPat)'; (1:nPat)'];               % same subject ID across time

            valid = (Y ~= missing_code) & ~isnan(Y) & ~isnan(FNC);
            if nnz(valid) < 5
                continue;
            end

            tbl = table(Y(valid), FNC(valid), categorical(Time(valid)), categorical(Subj(valid)), ...
                'VariableNames', {'Y','FNC','Time','Subject'});

            lme = fitlme(tbl, 'Y ~ FNC*Time + (1|Subject)');
            coef = lme.Coefficients;

            k = strcmp(coef.Name,'FNC');
            if any(k)
                beta_FNC(i,s) = coef.Estimate(k);
                p_FNC(i,s)    = coef.pValue(k);
            end

            k = strcmp(coef.Name,'Time_1');
            if any(k)
                beta_Time(i,s) = coef.Estimate(k);
                p_Time(i,s)    = coef.pValue(k);
            end

            k = strcmp(coef.Name,'FNC:Time_1');
            if any(k)
                beta_Int(i,s) = coef.Estimate(k);
                p_Int(i,s)    = coef.pValue(k);
            end
        end
    end

    % ---------- FDR correction (BH), per score and per effect ----------
    q_FNC  = nan(size(p_FNC));  sig_FNC  = false(size(p_FNC));
    q_Time = nan(size(p_Time)); sig_Time = false(size(p_Time));
    q_Int  = nan(size(p_Int));  sig_Int  = false(size(p_Int));

    for s = 1:nScores
        % FNC main
        pvec = p_FNC(:,s); v = ~isnan(pvec);
        if any(v)
            [h, ~, adj_p] = fdr_bh(pvec(v), alpha_fdr, 'pdep', 'no');
            q_FNC(v,s) = adj_p; sig_FNC(v,s) = h;
        end

        % Time main
        pvec = p_Time(:,s); v = ~isnan(pvec);
        if any(v)
            [h, ~, adj_p] = fdr_bh(pvec(v), alpha_fdr, 'pdep', 'no');
            q_Time(v,s) = adj_p; sig_Time(v,s) = h;
        end

        % Interaction
        pvec = p_Int(:,s); v = ~isnan(pvec);
        if any(v)
            [h, ~, adj_p] = fdr_bh(pvec(v), alpha_fdr, 'pdep', 'no');
            q_Int(v,s) = adj_p; sig_Int(v,s) = h;
        end
    end

    % ---------- Build summary table (all effects) ----------
    hasEdgeTable = exist('EdgeTable','var') && istable(EdgeTable) && height(EdgeTable) >= max(edgesForCognition) ...
                   && all(ismember({'Node1_name','Node2_name'}, EdgeTable.Properties.VariableNames));

    nRows = nEdgesUse * nScores;
    Score   = strings(nRows,1);
    EdgeID  = nan(nRows,1);
    Node1   = strings(nRows,1);
    Node2   = strings(nRows,1);

    BetaF = nan(nRows,1); PF = nan(nRows,1); QF = nan(nRows,1);
    BetaT = nan(nRows,1); PT = nan(nRows,1); QT = nan(nRows,1);
    BetaI = nan(nRows,1); PI = nan(nRows,1); QI = nan(nRows,1);

    r = 0;
    for s = 1:nScores
        for i = 1:nEdgesUse
            r = r + 1;
            edgeID = edgesForCognition(i);

            Score(r)  = string(scoreNames{s});
            EdgeID(r) = edgeID;

            if hasEdgeTable
                Node1(r) = string(EdgeTable.Node1_name{edgeID});
                Node2(r) = string(EdgeTable.Node2_name{edgeID});
            end

            BetaF(r) = beta_FNC(i,s); PF(r) = p_FNC(i,s); QF(r) = q_FNC(i,s);
            BetaT(r) = beta_Time(i,s); PT(r) = p_Time(i,s); QT(r) = q_Time(i,s);
            BetaI(r) = beta_Int(i,s);  PI(r) = p_Int(i,s);  QI(r) = q_Int(i,s);
        end
    end

    cog_summary = table(Score, EdgeID, Node1, Node2, ...
        BetaF, PF, QF, BetaT, PT, QT, BetaI, PI, QI, ...
        'VariableNames', {'Score','EdgeID','Node1','Node2', ...
                          'Beta_FNC','P_FNC','Q_FNC', ...
                          'Beta_Time','P_Time','Q_Time', ...
                          'Beta_FNCxTime','P_FNCxTime','Q_FNCxTime'});

    % "Joint edges": interaction survives FDR
    joint_edges = cog_summary(cog_summary.Q_FNCxTime < alpha_fdr, :);

end

%% -------------------------------------------------------------------------
% 8) Correlation between FNC and cognitive measures (optional summary)
% This section appends simple Pearson correlations to cog_summary for:
%   - Patients vs controls
%   - T0 vs T1
% using the same EdgeID and Score rows already present in cog_summary.
%
% Requires in workspace:
%   cog_summary (table with variables: Score, EdgeID)
%   zval.pat, zval.con
%   cognitive_patients, cognitive_controls
%   cfg.cog_score_idx, cfg.cog_score_names, cfg.cog_missing_code



missCode = cfg.cog_missing_code;

% Add output columns once (keeps the table tidy and avoids repeated variables)
outVars = {'R_pat_T0','P_pat_T0','R_pat_T1','P_pat_T1', ...
    'R_con_T0','P_con_T0','R_con_T1','P_con_T1'};
for v = 1:numel(outVars)
    if ~ismember(outVars{v}, cog_summary.Properties.VariableNames)
        cog_summary.(outVars{v}) = nan(height(cog_summary), 1);
    end
end

scoreNames = string(cfg.cog_score_names);
scoreIdx   = cfg.cog_score_idx;

for i = 1:height(cog_summary)

    e = cog_summary.EdgeID(i);
    sName = string(cog_summary.Score(i));

    % Map score name -> original cognitive index (e.g., 4 or 6)
    loc = find(scoreNames == sName, 1);
    if isempty(loc)
        continue;
    end
    s = scoreIdx(loc);

    % ---- Patients ----
    x = squeeze(zval.pat(e,:,1))';  y = squeeze(cognitive_patients(:,s,1));
    v = (y ~= missCode) & ~isnan(y) & ~isnan(x);
    if nnz(v) >= 3, [cog_summary.R_pat_T0(i), cog_summary.P_pat_T0(i)] = corr(x(v), y(v)); end

    x = squeeze(zval.pat(e,:,2))';  y = squeeze(cognitive_patients(:,s,2));
    v = (y ~= missCode) & ~isnan(y) & ~isnan(x);
    if nnz(v) >= 3, [cog_summary.R_pat_T1(i), cog_summary.P_pat_T1(i)] = corr(x(v), y(v)); end

    % ---- Controls ----
    x = squeeze(zval.con(e,:,1))';  y = squeeze(cognitive_controls(:,s,1));
    v = (y ~= missCode) & ~isnan(y) & ~isnan(x);
    if nnz(v) >= 3, [cog_summary.R_con_T0(i), cog_summary.P_con_T0(i)] = corr(x(v), y(v)); end

    x = squeeze(zval.con(e,:,2))';  y = squeeze(cognitive_controls(:,s,2));
    v = (y ~= missCode) & ~isnan(y) & ~isnan(x);
    if nnz(v) >= 3, [cog_summary.R_con_T1(i), cog_summary.P_con_T1(i)] = corr(x(v), y(v)); end
end

%% -------------------------------------------------------------------------
% 9) Joint LME for selected edges (controls + patients) to identify FNC by
% Group by Time effect. 
% Model: Y ~ FNC * Group * Time + Age + Sex + (1|ID)
% - Y: cognitive score (choose index, e.g., 4=WCST non-persev)
% - FNC: edge-wise z connectivity (zval.con/zval.pat)
% - Group: Control vs Patient
% - Time: T0 vs T1
%
% Requires:
%   joint_edges (table) with column 'EdgeID'  OR  jointedges.EdgeID
%   zval.con: [nEdges x nCon x 2], zval.pat: [nEdges x nPat x 2]
%   cognitive_controls: [nCon x 6 x 2], cognitive_patients: [nPat x 6 x 2]
%   age.con, age.pat (nCon/nPat), sex.con, sex.pat (same length; categorical or numeric)
%   EdgeTable optional (for printing labels)
%
% NOTE: missing cognitive values coded as cfg.cog_missing_code (e.g., 999)

% pick which cognitive score to use in this joint model
cogIdx_joint = 4;  % e.g., 4 = WCST non-perseverative error
missCode = cfg.cog_missing_code;

% interaction significant (uncorrected)  -> jointedges
jointedges = cog_summary(cog_summary.P_FNCxTime < 0.05, :);
edgeList = jointedges.EdgeID;



fprintf('\n=== Running Joint LME (FNC*Group*Time) on selected edges ===\n');
fprintf('Number of edges to test: %d\n\n', numel(edgeList));

hasEdgeTable = exist('EdgeTable','var') && istable(EdgeTable) && height(EdgeTable) >= max(edgeList) ...
               && all(ismember({'Node1_name','Node2_name'}, EdgeTable.Properties.VariableNames));

for k = 1:numel(edgeList)

    edgeID = edgeList(k);

    % Identify nodes (optional)
    if hasEdgeTable
        node1 = EdgeTable.Node1_name{edgeID};
        node2 = EdgeTable.Node2_name{edgeID};
    else
        node1 = 'NA'; node2 = 'NA';
    end

    fprintf('--- Joint LME Effects for EdgeID %d (%s — %s) ---\n', edgeID, node1, node2);

    % FNC vectors (stack in this order: con T0, con T1, pat T0, pat T1)
    fnc_con_T0 = squeeze(zval.con(edgeID,:,1))';   % [nCon x 1]
    fnc_con_T1 = squeeze(zval.con(edgeID,:,2))';   % [nCon x 1]
    fnc_pat_T0 = squeeze(zval.pat(edgeID,:,1))';   % [nPat x 1]
    fnc_pat_T1 = squeeze(zval.pat(edgeID,:,2))';   % [nPat x 1]
    FNC = [fnc_con_T0; fnc_con_T1; fnc_pat_T0; fnc_pat_T1];

    % Cognitive Y stacked in the SAME order
    Y = [ ...
        cognitive_controls(:, cogIdx_joint, 1); ...
        cognitive_controls(:, cogIdx_joint, 2); ...
        cognitive_patients(:, cogIdx_joint, 1); ...
        cognitive_patients(:, cogIdx_joint, 2)];

    % Group label (same length as FNC/Y)
    Group = categorical([ ...
        repmat({'Control'}, nCon, 1); ...
        repmat({'Control'}, nCon, 1); ...
        repmat({'Patient'}, nPat, 1); ...
        repmat({'Patient'}, nPat, 1)]);

    % Time label
    Time = categorical([ ...
        repmat({'T0'}, nCon, 1); ...
        repmat({'T1'}, nCon, 1); ...
        repmat({'T0'}, nPat, 1); ...
        repmat({'T1'}, nPat, 1)]);

    % Unique subject IDs across groups (controls 1..nCon, patients nCon+1..nCon+nPat)
    ID = [ ...
        (1:nCon)'; (1:nCon)'; ...
        (nCon + (1:nPat))'; (nCon + (1:nPat))' ];

    % Covariates (repeat per timepoint in the SAME stacking order)
    Age = [age.con(:); age.con(:); age.pat(:); age.pat(:)];
    Sex = [sex.con(:); sex.con(:); sex.pat(:); sex.pat(:)];  % ok if already categorical/numeric

    % Remove missing
    valid = (Y ~= missCode) & ~isnan(Y) & ~isnan(FNC) & ~isnan(Age);

    tbl = table(Y(valid), FNC(valid), Age(valid), Sex(valid), Group(valid), Time(valid), categorical(ID(valid)), ...
        'VariableNames', {'Y','FNC','Age','Sex','Group','Time','ID'});

    lme = fitlme(tbl, 'Y ~ FNC * Group * Time + Age + Sex + (1|ID)');
    coef = lme.Coefficients;

    % Safer extraction (in case a term is missing)
    getP = @(name) (coef.pValue(find(strcmp(coef.Name,name), 1, 'first')));
    pGroup = getP('Group_Patient');
    pTime  = getP('Time_T1');
    pFNC   = getP('FNC');
    pFG    = getP('FNC:Group_Patient');
    pFT    = getP('FNC:Time_T1');
    pGT    = getP('Group_Patient:Time_T1');
    pFGT   = getP('FNC:Group_Patient:Time_T1');

    fprintf('Group main effect:          p = %.4f\n', pGroup);
    fprintf('Time main effect:           p = %.4f\n', pTime);
    fprintf('FNC main effect:            p = %.4f\n', pFNC);
    fprintf('FNC × Group:                p = %.4f\n', pFG);
    fprintf('FNC × Time:                 p = %.4f\n', pFT);
    fprintf('Group × Time:               p = %.4f\n', pGT);
    fprintf('FNC × Group × Time:         p = %.4f\n\n', pFGT);

end




%% Done
if exist('EdgeTable','var')
    save(fullfile(cfg.output_dir,'EdgeTable.mat'),'EdgeTable');
end

if exist('sig_edges_unc_G','var')
    save(fullfile(cfg.output_dir,'sig_edges_unc_G.mat'),'sig_edges_unc_G');
end

disp('FNC analysis complete.');

