%% DTI-based structural connectivity: LME + cognition + plots
% This script:
%   1) Loads subject-wise DTI network matrices (14x14) for patients & controls at T0/T1
%   2) Extracts upper-triangle edges (91 edges)
%   3) Runs linear mixed-effects models:
%           Var1 ~ Age + Sex + Group*Time + (1|ID)
%   4) Optionally examines correlations with cognitive scores
%   5) Generates boxplots for selected edges
%
% Required inputs (to be prepared by the user):
%   - Base directory with group-wise subfolders and subject folders
%   - Age / sex text files
%   - DTI network matrices named as: <SubID>_<TimePoint>_fdt_network_matrix.txt
%   - Cognitive data (cognitive_patients, cognitive_controls) if correlation is needed
%
% -------------------------------------------------------------------------
%% Housekeeping
clearvars; clc;

%% -------------------------------------------------------------------------
% 0) Configuration
cfg = struct();

% Root directory of the repository (assumes this script is in repo_root/matlab/)
cfg.repo_root   = fileparts(fileparts(mfilename('fullpath')));
cfg.input_dir   = fullfile(cfg.repo_root, 'inputs');
cfg.output_dir  = fullfile(cfg.repo_root, 'outputs');
cfg.dti_dir     = fullfile(cfg.input_dir, 'dti_network'); % DTI txt 파일들이 있는 곳

% Input files
cfg.age_pat_file   = fullfile(cfg.input_dir, 'age_patients.txt');
cfg.age_con_file   = fullfile(cfg.input_dir, 'age_controls.txt');
cfg.sex_pat_file   = fullfile(cfg.input_dir, 'sex_patients.txt');
cfg.sex_con_file   = fullfile(cfg.input_dir, 'sex_controls.txt');
cfg.netinfo_file   = fullfile(cfg.input_dir, 'NetInfo.txt'); % must contain network name

cfg.cog_pat_file       = fullfile(cfg.input_dir, 'cognitive_patients.mat'); % "cognitive_patients"
cfg.cog_con_file       = fullfile(cfg.input_dir, 'cognitive_controls.mat'); % "cognitive_controls"
cfg.cog_missing_code   = 999;
cfg.cog_score_idx      = [4 6]; % [WCST non-perseverative error, backward digit span]
cfg.cog_score_names    = {'WCST_nonpersev', 'DigitSpan_backward'};

% Statistics
cfg.alpha_uncorrected  = 0.05;
cfg.alpha_fdr          = 0.05;

%% 1) Load Demographics & ID lists
fprintf('>> Loading demographics...\n');
age.con = load(cfg.age_con_file); age.pat = load(cfg.age_pat_file);
sex.con = load(cfg.sex_con_file); sex.pat = load(cfg.sex_pat_file);
netNames = importdata(cfg.netinfo_file);
nCon = numel(age.con); nPat = numel(age.pat);
N = numel(netNames);
mask = triu(true(N), 1); 
nEdges = nnz(mask);



% Subject ID lists
% You can either fill these manually or load from a text file.
% Example below uses manual cell arrays (same as original code).

subj_pat = {'Sub01','Sub02','Sub03','Sub04','Sub05','Sub06','Sub07','Sub09','Sub10','Sub11','Sub12','Sub13','Sub14','Sub15','Sub16','Sub18','Sub19','Sub21','Sub22','Sub23', ...
    'Sub24','Sub25','Sub26','Sub27','Sub28','Sub29','Sub30','Sub31','Sub33','Sub34','Sub35','Sub36','Sub37','SubR01','SubR02','SubR03','SubR04','SubR07','SubR10','SubR12','SubR14'};  % nPat

subj_con = {'Con01','Con02','Con03','Con04','Con05','Con06','Con08','Con10','Con11','Con15','Con16','Con17','Con18','Con19','Con20','Con21','Con22','Con23','Con24','Con25', ...
    'Con26','Con27','Con28','Con30','Con31','Con35','Con36','Con37','Con38','Con40', 'Con41','Con42','N001','N002','N005'};  % nCon
groupName = {'Patients', 'Controls'};
timepoint = {'T0', 'T1'};


%% ------------------------------------------------------------------------
% 3. Load DTI network matrices
%    For each subject and timepoint:
%       <dti_dir>/<groupName>/<SubID>/<TimePoint>/<SubID>_<TimePoint>_fdt_network_matrix.txt
%    Each file: [Dim x Dim] adjacency (e.g., streamline count)
% -------------------------------------------------------------------------

fprintf('>> Loading DTI network matrices...\n');
scval.con = nan(nEdges, nCon, 2);
scval.pat = nan(nEdges, nPat, 2);

for t = 1:2
    tp = timepoint{t};
    % Patients
    for s = 1:nPat
        sid = subj_pat{s};
        fPath = fullfile(cfg.dti_dir, groupName{1}, sid, tp, sprintf('%s_%s_fdt_network_matrix.txt', sid, tp));
        mat = load(fPath);
        scval.pat(:,s,t) = mat(mask); % 복잡한 reshape 없이 바로 mask 적용
    end
    % Controls
    for s = 1:nCon
        sid = subj_con{s};
        fPath = fullfile(cfg.dti_dir, groupName{2}, sid, tp, sprintf('%s_%s_fdt_network_matrix.txt', sid, tp));
        mat = load(fPath);
        scval.con(:,s,t) = mat(mask);
    end
end

%% ------------------------------------------------------------------------
[ix_row, ix_col] = find(mask);   % IMPORTANT: same order as mat(mask)
EdgeTable = table( (1:nEdges)', ix_row, ix_col, ...
    netNames(ix_row), netNames(ix_col), ...
    'VariableNames', {'EdgeID','Node1_idx','Node2_idx','Node1_name','Node2_name'} );

%% ------------------------------------------------------------------------
% 5. Linear Mixed-Effects Model analysis + summary
%     Var1 ~ Age + Sex + Group * Time + (1 | Individuals)
%     net_pat, net_con: [nEdges × nSub × 2(timepoints)]
% -------------------------------------------------------------------------

G   = nan(nEdges,1);   % p(Group)
T   = nan(nEdges,1);   % p(Time)
Int = nan(nEdges,1);   % p(Group×Time)

beta_G   = nan(nEdges,1);
beta_T   = nan(nEdges,1);
beta_Int = nan(nEdges,1);

t_G   = nan(nEdges,1);
t_T   = nan(nEdges,1);
t_Int = nan(nEdges,1);

Group = categorical([repmat({'Control'}, nCon*2, 1); repmat({'Patient'}, nPat*2, 1)]);
Time  = categorical([repmat({'T0'}, nCon, 1); repmat({'T1'}, nCon, 1); ...
                     repmat({'T0'}, nPat, 1); repmat({'T1'}, nPat, 1)]);
Subject = [repmat((1:nCon)', 2, 1); repmat((nCon+1:nCon+nPat)', 2, 1)];

Age = [age.con(:); age.con(:); age.pat(:); age.pat(:)];
Sex = [sex.con(:); sex.con(:); sex.pat(:); sex.pat(:)];



output = [];    % uncorrected p<0.05 (Group or Time) edge index

for e = 1:nEdges
    y_con = squeeze(scval.con(e,:,:)); % nCon x 2
    y_pat = squeeze(scval.pat(e,:,:)); % nPat x 2
    y = [y_con(:); y_pat(:)];

    tbl = table(y, Age, Sex, Group, Time, categorical(Subject), ...
        'VariableNames', {'SC','Age','Sex','Group','Time','Subject'});

    lme  = fitlme(tbl, 'SC ~ Age + Sex + Group*Time + (1|Subject)');
    coef = lme.Coefficients;


    % Coefficient order: (Intercept) Age Sex Group_Patient Time_Time2 Group_Patient:Time_Time2
    G(e)   = coef.pValue(4);   % Group main effect
    T(e)   = coef.pValue(5);   % Time main effect
    Int(e) = coef.pValue(6);   % Group×Time interaction

    beta_G(e)   = coef.Estimate(4);
    beta_T(e)   = coef.Estimate(5);
    beta_Int(e) = coef.Estimate(6);

    t_G(e)   = coef.tStat(4);
    t_T(e)   = coef.tStat(5);
    t_Int(e) = coef.tStat(6);

    if any([G(e), T(e)] < cfg.alpha_uncorrected)
        output = [output; e];
    end
end

%% ------------------------------------------------------------------------
% 5-1. FDR correction (Benjamini–Hochberg)
% -------------------------------------------------------------------------
[hG,  G_fdr]   = fdr_bh(G);     % Group main effect
[hT,  T_fdr]   = fdr_bh(T);     % Time main effect
[hInt, Int_fdr] = fdr_bh(Int);  % Group×Time interaction

%% ------------------------------------------------------------------------
% 5-2. Summary table (원하면 output만 사용하거나, 전체 nEdges 사용 가능)
% -------------------------------------------------------------------------
EdgeID = find( (G < cfg.alpha_uncorrected) | (T < cfg.alpha_uncorrected) | (Int < cfg.alpha_uncorrected) );
Node1  = EdgeTable.Node1_name(EdgeID);
Node2  = EdgeTable.Node2_name(EdgeID);

LME_summary = table(EdgeID, Node1, Node2, beta_G(EdgeID), t_G(EdgeID), G(EdgeID),...
    beta_T(EdgeID),   t_T(EdgeID),   T(EdgeID), beta_Int(EdgeID), t_Int(EdgeID), Int(EdgeID), ...
    'VariableNames', {'EdgeID','Node1','Node2', ...
                      'Beta_Group','t_Group','p_Group', ...
                      'Beta_Time','t_Time','p_Time', ...
                      'Beta_Interaction','t_Interaction','p_Interaction'});

disp(LME_summary)

%% ------------------------------------------------------------------------
% 5-bis. Identify time-only edges and define edgesForCognition
% -------------------------------------------------------------------------
isGroupSig = G < cfg.alpha_uncorrected;   % edges with significant Group main effect
isTimeSig  = T < cfg.alpha_uncorrected;   % edges with significant Time main effect

% Edges where ONLY Time is significant (Group not significant)
edgesTimeOnly = find(~isGroupSig & isTimeSig);

% Edges where either Group or Time was significant (same as "output")
edgesAny = output;

% Final set of edges used for cognition analyses:
%   = significant edges minus "time-only" edges
edgesForCognition = setdiff(edgesAny, edgesTimeOnly);



%% ------------------------------------------------------------------------
% 6. Post-hoc analysis for edgesForCognition
%    Tests (per edge):
%      1) Patient vs Control @ T0  (ttest2)
%      2) Patient vs Control @ T1  (ttest2)
%      3) Within-patient: T0 vs T1 (paired ttest)
%      4) Within-control: T0 vs T1 (paired ttest)
% -------------------------------------------------------------------------

nSig = numel(edgesAny);

posth = struct();
posth.edgeID = edgesAny;
posth.p      = nan(nSig, 4);
posth.t      = nan(nSig, 4);

for k = 1:nSig
    e = edgesAny(k);

    % Structural connectivity for this edge (already [nSub x 1])
    conT0 = squeeze(scval.con(e,:,1))';   % [nCon x 1]
    conT1 = squeeze(scval.con(e,:,2))';   % [nCon x 1]
    patT0 = squeeze(scval.pat(e,:,1))';   % [nPat x 1]
    patT1 = squeeze(scval.pat(e,:,2))';   % [nPat x 1]

    % 1) Patient vs Control @ T0
    [~, p, ~, stats] = ttest2(patT0, conT0);
    posth.p(k,1) = p;
    posth.t(k,1) = stats.tstat;

    % 2) Patient vs Control @ T1
    [~, p, ~, stats] = ttest2(patT1, conT1);
    posth.p(k,2) = p;
    posth.t(k,2) = stats.tstat;

    % 3) Within-patient: T0 vs T1 (paired)
    [~, p, ~, stats] = ttest(patT0, patT1);
    posth.p(k,3) = p;
    posth.t(k,3) = stats.tstat;

    % 4) Within-control: T0 vs T1 (paired)
    [~, p, ~, stats] = ttest(conT0, conT1);
    posth.p(k,4) = p;
    posth.t(k,4) = stats.tstat;
end


%% ------------------------------------------------------------------------
% 7) Cognition models (patients only) - DTI connectivity
% Cognitive matrices:
%   cognitive_patients : [nPat x nScoresTotal x 2] (dim3: 1=T0, 2=T1)
%
% Patients-only LME (edge-wise):
%   Y ~ Conn * Time + (1|Subject)
%
% Requires in workspace:
%   scval.pat         : [nEdges x nPat x 2]   (DTI edge values)
%   edgesForCognition : vector of edge IDs to test
%   EdgeTable         : optional (for node labels; must match scval edge ordering)

edgesForCognition = edgesForCognition(:);
nEdgesUse = numel(edgesForCognition);

% ---- Settings (reuse cfg like fMRI) ----
scoreIdx     = cfg.cog_score_idx;      % [4 6]
scoreNames   = cfg.cog_score_names;    % {'WCST_nonpersev','DigitSpan_backward'}
missing_code = cfg.cog_missing_code;   % 999
alpha_fdr    = cfg.alpha_fdr;
Cpat = load(cfg.cog_pat_file);
cognitive_patients = Cpat.cognitive_patients;
Ccon = load(cfg.cog_con_file);
cognitive_controls = Ccon.cognitive_controls;

nScores = numel(scoreIdx);

% ---- Preallocate ----
beta_Conn = nan(nEdgesUse, nScores);
beta_Time = nan(nEdgesUse, nScores);
beta_Int  = nan(nEdgesUse, nScores);

p_Conn = nan(nEdgesUse, nScores);
p_Time = nan(nEdgesUse, nScores);
p_Int  = nan(nEdgesUse, nScores);

% ---- Edge-wise patients-only LME ----
for s = 1:nScores
    cogIdx = scoreIdx(s);

    y_T0 = cognitive_patients(:, cogIdx, 1);
    y_T1 = cognitive_patients(:, cogIdx, 2);

    for i = 1:nEdgesUse
        edgeID = edgesForCognition(i);

        conn_T0 = squeeze(scval.pat(edgeID, :, 1))';  % [nPat x 1]
        conn_T1 = squeeze(scval.pat(edgeID, :, 2))';  % [nPat x 1]

        % Long format (stack T0/T1)
        Y    = [y_T0; y_T1];
        Conn = [conn_T0; conn_T1];
        Time = [zeros(nPat,1); ones(nPat,1)];         % 0=T0, 1=T1
        Subj = [(1:nPat)'; (1:nPat)'];                % subject ID repeated

        valid = (Y ~= missing_code) & ~isnan(Y) & ~isnan(Conn);
        if nnz(valid) < 5
            continue;
        end

        tbl = table(Y(valid), Conn(valid), categorical(Time(valid)), categorical(Subj(valid)), ...
            'VariableNames', {'Y','Conn','Time','Subject'});

        lme  = fitlme(tbl, 'Y ~ Conn*Time + (1|Subject)');
        coef = lme.Coefficients;

        k = strcmp(coef.Name,'Conn');
        if any(k)
            beta_Conn(i,s) = coef.Estimate(k);
            p_Conn(i,s)    = coef.pValue(k);
        end

        k = strcmp(coef.Name,'Time_1');
        if any(k)
            beta_Time(i,s) = coef.Estimate(k);
            p_Time(i,s)    = coef.pValue(k);
        end

        k = strcmp(coef.Name,'Conn:Time_1');
        if any(k)
            beta_Int(i,s) = coef.Estimate(k);
            p_Int(i,s)    = coef.pValue(k);
        end
    end
end

% ---- FDR correction (BH): Conn / Time / Interaction, per score ----
q_Conn = nan(size(p_Conn)); sig_Conn = false(size(p_Conn));
q_Time = nan(size(p_Time)); sig_Time = false(size(p_Time));
q_Int  = nan(size(p_Int));  sig_Int  = false(size(p_Int));

for s = 1:nScores
    v = ~isnan(p_Conn(:,s));
    if any(v)
        [h, ~, adj_p] = fdr_bh(p_Conn(v,s), alpha_fdr, 'pdep', 'no');
        q_Conn(v,s) = adj_p; sig_Conn(v,s) = h;
    end

    v = ~isnan(p_Time(:,s));
    if any(v)
        [h, ~, adj_p] = fdr_bh(p_Time(v,s), alpha_fdr, 'pdep', 'no');
        q_Time(v,s) = adj_p; sig_Time(v,s) = h;
    end

    v = ~isnan(p_Int(:,s));
    if any(v)
        [h, ~, adj_p] = fdr_bh(p_Int(v,s), alpha_fdr, 'pdep', 'no');
        q_Int(v,s) = adj_p; sig_Int(v,s) = h;
    end
end

% ---- Build summary table (same layout as fMRI cog_summary) ----
hasEdgeTable = exist('EdgeTable','var') && istable(EdgeTable) && height(EdgeTable) >= max(edgesForCognition) && ...
               all(ismember({'Node1_name','Node2_name'}, EdgeTable.Properties.VariableNames));

nRows = nEdgesUse * nScores;
Score  = strings(nRows,1);
EdgeID = nan(nRows,1);
Node1  = strings(nRows,1);
Node2  = strings(nRows,1);

BetaC = nan(nRows,1); PC = nan(nRows,1); QC = nan(nRows,1);
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

        BetaC(r) = beta_Conn(i,s); PC(r) = p_Conn(i,s); QC(r) = q_Conn(i,s);
        BetaT(r) = beta_Time(i,s); PT(r) = p_Time(i,s); QT(r) = q_Time(i,s);
        BetaI(r) = beta_Int(i,s);  PI(r) = p_Int(i,s);  QI(r) = q_Int(i,s);
    end
end

cog_summary = table(Score, EdgeID, Node1, Node2, ...
    BetaC, PC, QC, BetaT, PT, QT, BetaI, PI, QI, ...
    'VariableNames', {'Score','EdgeID','Node1','Node2', ...
                      'Beta_Conn','P_Conn','Q_Conn', ...
                      'Beta_Time','P_Time','Q_Time', ...
                      'Beta_ConnxTime','P_ConnxTime','Q_ConnxTime'});

% "Joint edges": interaction survives FDR
joint_edges = cog_summary(cog_summary.Q_ConnxTime < alpha_fdr, :);



%% -------------------------------------------------------------------------
% 6b) Correlation between structural connectivity and cognitive measures
%
% Requires in workspace:
%   scval.pat : [nEdges x nPat x 2]   (DTI connectivity, e.g., streamline count)
%   scval.con : [nEdges x nCon x 2]
%   cognitive_patients : [nPat x nScoresTotal x 2]
%   cognitive_controls : [nCon x nScoresTotal x 2]
%   edgesForCognition  : [nEdgesCog x 1] edge IDs (subset)
%   cfg.cog_score_idx, cfg.cog_score_names, cfg.cog_missing_code
% -------------------------------------------------------------------------

cogIdxCorr   = cfg.cog_score_idx(:)';          % e.g., [4 6]
scoreNamesDTI = string(cfg.cog_score_names(:)); % e.g., ["WCST_nonpersev","DigitSpan_backward"]
missCode     = cfg.cog_missing_code;

nScoresCorr  = numel(cogIdxCorr);
edgesForCognition = edgesForCognition(:);
nEdgesCog    = numel(edgesForCognition);

r_cog_DTI = nan(nScoresCorr, nEdgesCog, 2, 2);   % score x edge x group x time
p_cog_DTI = nan(nScoresCorr, nEdgesCog, 2, 2);

for j = 1:nEdgesCog
    edgeID = edgesForCognition(j);

    for g = 1:2              % group: 1=patients, 2=controls
        for t = 1:2          % time : 1=T0, 2=T1

            if g == 1
                conn_all = squeeze(scval.pat(edgeID,:,t))';   % [nPat x 1]
                cog_mat  = cognitive_patients(:,:,t);         % [nPat x nScoresTotal]
            else
                conn_all = squeeze(scval.con(edgeID,:,t))';   % [nCon x 1]
                cog_mat  = cognitive_controls(:,:,t);         % [nCon x nScoresTotal]
            end

            for k = 1:nScoresCorr
                sOrig = cogIdxCorr(k);                        % original cognitive index (e.g., 4 or 6)
                score = cog_mat(:, sOrig);

                valid = (score ~= missCode) & ~isnan(score) & ~isnan(conn_all);
                if nnz(valid) >= 3
                    [r_cog_DTI(k,j,g,t), p_cog_DTI(k,j,g,t)] = corr(conn_all(valid), score(valid));
                end
            end
        end
    end
end


%% -------------------------------------------------------------------------
% 6b-bis) Append DTI connectivity–cognition correlations to cog_summary
%
% Requires in workspace:
%   cog_summary (table) with variables: Score, EdgeID
%   r_cog_DTI, p_cog_DTI, edgesForCognition
% -------------------------------------------------------------------------

outVarsDTI = {'Rpat_T0_DTI','Ppat_T0_DTI','Rpat_T1_DTI','Ppat_T1_DTI', ...
              'Rcon_T0_DTI','Pcon_T0_DTI','Rcon_T1_DTI','Pcon_T1_DTI'};
for v = 1:numel(outVarsDTI)
    if ~ismember(outVarsDTI{v}, cog_summary.Properties.VariableNames)
        cog_summary.(outVarsDTI{v}) = nan(height(cog_summary), 1);
    end
end

ScoreCol = string(cog_summary.Score);

for i = 1:height(cog_summary)

    edgeID = cog_summary.EdgeID(i);

    % Map Score -> k (must match cfg.cog_score_names ordering)
    k = find(scoreNamesDTI == ScoreCol(i), 1);
    if isempty(k), continue; end

    % Map EdgeID -> j (position in edgesForCognition)
    j = find(edgesForCognition == edgeID, 1);
    if isempty(j), continue; end

    % Patients (group=1)
    cog_summary.Rpat_T0_DTI(i) = r_cog_DTI(k, j, 1, 1);
    cog_summary.Ppat_T0_DTI(i) = p_cog_DTI(k, j, 1, 1);
    cog_summary.Rpat_T1_DTI(i) = r_cog_DTI(k, j, 1, 2);
    cog_summary.Ppat_T1_DTI(i) = p_cog_DTI(k, j, 1, 2);

    % Controls (group=2)
    cog_summary.Rcon_T0_DTI(i) = r_cog_DTI(k, j, 2, 1);
    cog_summary.Pcon_T0_DTI(i) = p_cog_DTI(k, j, 2, 1);
    cog_summary.Rcon_T1_DTI(i) = r_cog_DTI(k, j, 2, 2);
    cog_summary.Pcon_T1_DTI(i) = p_cog_DTI(k, j, 2, 2);
end
