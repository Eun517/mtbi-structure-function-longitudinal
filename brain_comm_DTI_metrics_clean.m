%% Associations between DTI metrics (FA, AD, RD, and MD) and cognition 
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
cfg.wmint_file     = fullfile(cfg.input_dir, 'WMint.mat'); % must contain DTI metrics "WMint" 
cfg.EdgeTable    = fullfile(cfg.input_dir, 'EdgeTable.mat'); 
cfg.sig_edges    = fullfile(cfg.input_dir, 'sig_edges_unc_G.mat');  

% WMint is ordered by EdgeTable which was previously defined in
% brain_comm_fnc_clean.m
cfg.netinfo_file   = fullfile(cfg.input_dir, 'NetInfo.txt'); % must contain network name
cfg.cog_pat_file       = fullfile(cfg.input_dir, 'cognitive_patients.mat'); % "cognitive_patients"
cfg.cog_con_file       = fullfile(cfg.input_dir, 'cognitive_controls.mat'); % "cognitive_controls"
cfg.cog_missing_code   = 999;
cfg.cog_score_idx      = [4 6]; % [WCST non-perseverative error, backward digit span]
cfg.cog_score_names    = {'WCST_nonpersev', 'DigitSpan_backward'};

% Statistics
cfg.alpha_uncorrected  = 0.05;
cfg.alpha_fdr          = 0.05;

%% ------------------------------------------------------------------------
%% 1) Load Demographics & ID lists
fprintf('>> Loading demographics...\n');
S = load(cfg.wmint_file);
assert(isfield(S,'WMint'), 'WMint.mat must contain variable "WMint".');
WMint = S.WMint;             % WMint: nSub × metric × edge × time
% WMint is ordered by EdgeTable which showed significant main effects of
% group from the FNC analyses

Ccon = load(cfg.cog_con_file);
cognitive_controls = Ccon.cognitive_controls;
Cpat = load(cfg.cog_pat_file);
cognitive_patients = Cpat.cognitive_patients;

% Basic dimensions
nCon = size(cognitive_controls, 1);
nPat = size(cognitive_patients, 1);
nSub = size(WMint, 1);

idxCon = 1:nCon;
idxPat = (nCon+1):nSub;

% Analysis settings
scoreIdx   = cfg.cog_score_idx(:)';
scoreNames = string(cfg.cog_score_names(:)');
missCode   = cfg.cog_missing_code;
alpha_fdr  = cfg.alpha_fdr;
alpha_unc  = cfg.alpha_uncorrected;

nScores  = numel(scoreIdx);
metricNames = string({'FA','AD','RD','MD'});
nMetrics = numel(metricNames);
nEdges   = size(WMint, 3);


%% 2) load EdgeTable (optional, for node labels)
TT = load(cfg.EdgeTable);
EdgeTable = TT.EdgeTable;
sig_edges = load(cfg.sig_edges);
sig_edges = sig_edges.sig_edges_unc_G;

EdgeTable_sig = EdgeTable(ismember(EdgeTable.EdgeID, sig_edges), :);

%% -------------------------------------------------------------------------
% 3) Patients-only LME: WM × Time
% Model: Y ~ WM * Time + (1|Subject)
fprintf('>> Running patients-only LME for DTI metrics...\n');

scoreIdx   = cfg.cog_score_idx;
scoreNames = cfg.cog_score_names;
nScores    = numel(scoreIdx);

beta_WM   = nan(nEdges, nMetrics, nScores);
beta_Time = nan(nEdges, nMetrics, nScores);
beta_Int  = nan(nEdges, nMetrics, nScores);

p_WM   = nan(nEdges, nMetrics, nScores);
p_Time = nan(nEdges, nMetrics, nScores);
p_Int  = nan(nEdges, nMetrics, nScores);

for s = 1:nScores
    cogIdx = scoreIdx(s);

    y_T0 = cognitive_patients(:, cogIdx, 1);
    y_T1 = cognitive_patients(:, cogIdx, 2);

    for m = 1:nMetrics
        for e = 1:nEdges

            wm_T0 = squeeze(WMint(idxPat, m, e, 1));
            wm_T1 = squeeze(WMint(idxPat, m, e, 2));

            Y    = [y_T0; y_T1];
            WMv  = [wm_T0; wm_T1];
            Time = [zeros(nPat,1); ones(nPat,1)];
            ID   = [(1:nPat)'; (1:nPat)'];

            valid = (Y ~= cfg.cog_missing_code) & ~isnan(Y) & ~isnan(WMv);
            if nnz(valid) < 5, continue; end

            tbl = table(Y(valid), WMv(valid), categorical(Time(valid)), categorical(ID(valid)), ...
                'VariableNames', {'Y','WM','Time','Subject'});

            lme = fitlme(tbl, 'Y ~ WM*Time + (1|Subject)');
            coef = lme.Coefficients;

            if any(strcmp(coef.Name,'WM'))
                beta_WM(e,m,s) = coef.Estimate(strcmp(coef.Name,'WM'));
                p_WM(e,m,s)    = coef.pValue(strcmp(coef.Name,'WM'));
            end
            if any(strcmp(coef.Name,'Time_1'))
                beta_Time(e,m,s) = coef.Estimate(strcmp(coef.Name,'Time_1'));
                p_Time(e,m,s)    = coef.pValue(strcmp(coef.Name,'Time_1'));
            end
            if any(strcmp(coef.Name,'WM:Time_1'))
                beta_Int(e,m,s) = coef.Estimate(strcmp(coef.Name,'WM:Time_1'));
                p_Int(e,m,s)    = coef.pValue(strcmp(coef.Name,'WM:Time_1'));
            end
        end
    end
end

%% -------------------------------------------------------------------------
% 4) FDR correction (patients-only)
fprintf('>> Applying FDR correction...\n');

q_WM   = nan(size(p_WM));
q_Time = nan(size(p_Time));
q_Int  = nan(size(p_Int));

for s = 1:nScores
    for m = 1:nMetrics
        pvec = squeeze(p_Int(:,m,s));
        v = ~isnan(pvec);
        if any(v)
            [~,~,adj_p] = fdr_bh(pvec(v), cfg.alpha_fdr, 'pdep', 'no');
            q_Int(v,m,s) = adj_p;
        end
    end
end

%% -------------------------------------------------------------------------
% 5) Build summary table
fprintf('>> Building summary table...\n');

EdgeID = repmat(sig_edges(:), nMetrics*nScores, 1);
Metric = strings(numel(EdgeID),1);
Score  = strings(numel(EdgeID),1);

BetaWM = nan(numel(EdgeID),1);
P_WM   = nan(numel(EdgeID),1);
P_Int  = nan(numel(EdgeID),1);
Q_Int  = nan(numel(EdgeID),1);

r = 0;
for s = 1:nScores
    for m = 1:nMetrics
        for i = 1:numel(sig_edges)
            r = r + 1;
            Metric(r) = metricNames{m};
            Score(r)  = scoreNames{s};

            BetaWM(r) = beta_WM(i,m,s);
            P_WM(r)   = p_WM(i,m,s);
            P_Int(r)  = p_Int(i,m,s);
            Q_Int(r)  = q_Int(i,m,s);
        end
    end
end

DTI_summary = table(Score, Metric, EdgeID, BetaWM, P_WM, P_Int, Q_Int, ...
    'VariableNames', {'Score','Metric','EdgeID','Beta_WM','P_WM','P_WMxTime','Q_WMxTime'});

%% -------------------------------------------------------------------------
% 6) Correlation between WM metrics and cognition (patients & controls)
fprintf('>> Computing WM–cognition correlations...\n');

nGroups = 2; % 1=patients, 2=controls
nTimes  = 2; % T0, T1

r_WM_cog = nan(nScores, nMetrics, nEdges, nGroups, nTimes);
p_WM_cog = nan(size(r_WM_cog));

for s = 1:nScores
    cogIdx = scoreIdx(s);
    for m = 1:nMetrics
        for e = 1:nEdges
            for g = 1:2
                for t = 1:2

                    if g == 1
                        wm = squeeze(WMint(idxPat, m, e, t));
                        cog = cognitive_patients(:, cogIdx, t);
                    else
                        wm = squeeze(WMint(idxCon, m, e, t));
                        cog = cognitive_controls(:, cogIdx, t);
                    end

                    valid = (cog ~= cfg.cog_missing_code) & ~isnan(cog) & ~isnan(wm);
                    if nnz(valid) >= 3
                        [r_WM_cog(s,m,e,g,t), p_WM_cog(s,m,e,g,t)] = corr(wm(valid), cog(valid));
                    end
                end
            end
        end
    end
end



%% ------------------------------------------------------------------------
% Append WM–cognition correlations (r/p) into DTI_summary
% Requires:
%   DTI_summary (Score, Metric, EdgeID)
%   r_WM_cog, p_WM_cog : [nScores x nMetrics x nEdges x 2(group) x 2(time)]
%   scoreNames, metricNames (same naming used when building DTI_summary)
% -------------------------------------------------------------------------

fprintf('>> Appending WM–cognition correlations into DTI_summary...\n');

% 1) Add columns once
outVars = {'R_pat_T0','P_pat_T0','R_pat_T1','P_pat_T1', ...
           'R_con_T0','P_con_T0','R_con_T1','P_con_T1'};
for v = 1:numel(outVars)
    if ~ismember(outVars{v}, DTI_summary.Properties.VariableNames)
        DTI_summary.(outVars{v}) = nan(height(DTI_summary),1);
    end
end

% 3) Fill row-wise using (Score, Metric, EdgeID)
for i = 1:height(DTI_summary)

    e = find(sig_edges == DTI_summary.EdgeID(i), 1);
    % score index (1..nScores)
    s = find(strcmp(scoreNames, string(DTI_summary.Score(i))), 1);

    % metric index (1..nMetrics)
    m = find(strcmp(metricNames, string(DTI_summary.Metric(i))), 1);

    % group: 1=patients, 2=controls / time: 1=T0, 2=T1
    DTI_summary.R_pat_T0(i) = r_WM_cog(s,m,e,1,1);
    DTI_summary.P_pat_T0(i) = p_WM_cog(s,m,e,1,1);
    DTI_summary.R_pat_T1(i) = r_WM_cog(s,m,e,1,2);
    DTI_summary.P_pat_T1(i) = p_WM_cog(s,m,e,1,2);

    DTI_summary.R_con_T0(i) = r_WM_cog(s,m,e,2,1);
    DTI_summary.P_con_T0(i) = p_WM_cog(s,m,e,2,1);
    DTI_summary.R_con_T1(i) = r_WM_cog(s,m,e,2,2);
    DTI_summary.P_con_T1(i) = p_WM_cog(s,m,e,2,2);
end

%% -------------------------------------------------------------------------
% 7) Save outputs
fprintf('>> Saving results...\n');
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

save(fullfile(cfg.output_dir,'DTI_metrics_LME_summary.mat'), 'DTI_summary');

fprintf('>> DTI metrics analysis completed.\n');
