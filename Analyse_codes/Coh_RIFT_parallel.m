%%%% using hilbert to compute coherence
% 1.select tagging response sensors (only from occipital); 2.get 2D coherence; 3.plot it!
%%% get the foveal coherence for each epochs

function Coh_RIFT_parallel(sid)
%%% set paths
server = 1; % run analysis on the server or local computer
[rootdir] = Get_Paths(server); %get the root dir and the colormap of cbrewer(RdBu)
PPath.ResPath = [rootdir 'Results' filesep 'Coh_parallel' filesep]; %result path
if (~exist(PPath.ResPath,'dir'))
    mkdir(PPath.ResPath);
end
%%% setting
freqrange = 40:1:80; % freq rang to run
freq_halfwidth = 5; % for hilbert filter
SenSelectP = 0.01; % sig threshold for selection tagging sensors
EpochType = {'PreTarg';'Targ';'PosTarg'}; % epochs of interest
loc2tar = [-1 0 1]; % loc2tar for {'Targ','PosTarg'}
%%% get the tag freq for pretarget and target
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
if ExpInfo.FreqVersion(sid) == 1
    tagfreq = [60 65];
else
    tagfreq = [65 60];
end
whichfreq = [2 1 2]; % analyze which tag freq
pds = {'MISC004','MISC005'}; % photodiodes for the first and second tagging freq (i.e., flicker target and flicker postartget)
%%% get the exp information
conds = ExpInfo.CondID;
condnm = [{'all'},ExpInfo.CondName];
cond_col = find(strcmp(ExpInfo.EventHdr,'SentenceCondition'));
tag_col = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
fixdur_col = find(strcmp(ExpInfo.EventHdr,'fixation_duration'));
subjects = ExpInfo.subjects;
sub = subjects{sid};
%%% set output
TagCoh = [];
TagCoh.CondName = condnm;
TagCoh.freq_halfwidth = freq_halfwidth;
TagCoh.EventHdr = ExpInfo.EventHdr;
TagCoh.subs = sub;
TagCoh.Coh_hdr = {'freq*time*sub*cond'};

%% %%%% ============= subjects run! ===============%%%%%%%%%
fprintf(['***** analyzing: s' num2str(sid) '**** \n\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 1. loading data %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPath.SaveData = [rootdir 'Analyse_data' filesep sub filesep];
load([PPath.SaveData 'epoch_WrdOn']); % epoch

%% nouth-filter the line noise around 50Hz
% already include this step in the pre-processing will be deleted for new
% datset
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = [49 51];
epoch = ft_preprocessing(cfg, epoch);

%%% normalize pd and remove trials that without photodies signal (just in case)
% get pd label index
pdi = [];
for p = 1:length(pds)
    pdi = [pdi find(strcmp(epoch.label,pds{p}))];
end
rmtrl = [];
for ttt = 1:length(epoch.trial)
    %%% remove trials that have no photodiode signal
    if max(max(epoch.trial{ttt}(pdi,:))) < 0.005
        rmtrl = [rmtrl; ttt];
    end
    epoch.trial{ttt}(pdi,:) = zscore(epoch.trial{ttt}(pdi,:),0,2);
end
cfg = [];
cfg.trials = setdiff(1:length(epoch.trial),rmtrl);
epoch = ft_selectdata(cfg, epoch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 2. getting significant tagging response sensors %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% permutation test between epoch data and baseline data
%%% run it seperately for different epochs->get different tag sensors
%%% selecting comparing epochs
% load in baseline data
load([PPath.SaveData 'epoch_BL_Cross']);
% run coherence for each occipital sensors(mag+gradio) (RIFT won't be in other brain regions)
load([rootdir 'OccipSens_full'])

for ep = 1:length(EpochType) % epoch index for sensor selection;{'PreTarg','Targ','PosTarg'};
    tf = whichfreq(ep); % the index of which freq/pds to use for this epoch type
    SigTagSens = cell(1,2); % [sensor_label sensor_id]
    % select the epoch data
    flickid     = epoch.trialinfo(:,tag_col)==loc2tar(ep);
    cfg         = [];
    cfg.trials  = find(flickid);
    data        = ft_selectdata(cfg, epoch);
    % remove bad sensors in occipital if there's any
    for i = 1:length(ExpInfo.BadSensor{sid})
        tmp = strcmp(OccipSens_full,['MEG' ExpInfo.BadSensor{sid}{i}]);
        if ~isempty(tmp)
            OccipSens_full(tmp) = [];
        end
    end
    
    % get fourier spctrm - no time domain (1D data) - easier to run statistics "ft_statfun_indepsamplesZcoh"
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.channel = [OccipSens_full;pds{tf}];
    cfg.method  = 'mtmfft';
    cfg.taper   = 'hanning';
    cfg.foi     = tagfreq(tf);
    cfg.pad     = 'nextpow2';
    fourier     = ft_freqanalysis(cfg,data);%=% tfr.powspctrm:3D data: trialnum*channelnum*freqpoint
    fourier_bl  = ft_freqanalysis(cfg,epoch_BL_Cross);
    
    %%% statistical test of coherence
    nchancom = length(OccipSens_full);
    stat_mask = zeros(nchancom,1);
    for i = 1:nchancom % run stat over sensors
        cfg            = [];
        cfg.channel    = {OccipSens_full{i},pds{tf}};
        fourier_tmp    = ft_selectdata(cfg, fourier);
        fourier_bl_tmp = ft_selectdata(cfg, fourier_bl);
        % do stat between tagging_period and baseline for each sensor
        cfg                  = [];
        cfg.parameter        = 'fourierspctrm';
        cfg.statistic        = 'ft_statfun_indepsamplesZcoh';  % take fourierspctrm as input, so no time domain information
        cfg.method           = 'montecarlo';
        cfg.tail             = 1; %% right sided, group1 is bigger than group2
        cfg.alpha            = SenSelectP;
        cfg.numrandomization = 2000;
        ntrl_1 = size(fourier.fourierspctrm,1);
        ntrl_2 = size(fourier_bl.fourierspctrm,1);
        design = zeros(1, ntrl_1 + ntrl_2);
        design(1,1:ntrl_1) = 1;
        design(1,(ntrl_1 +1):(ntrl_1 + ntrl_2))= 2;
        cfg.design = design;
        cfg.ivar   = 1;
        stat = ft_freqstatistics(cfg, fourier_tmp, fourier_bl_tmp);
        stat_mask(i) = any(any(stat.mask));
    end
    
    %%% get the sig occipital tagging sensors
    if any(stat_mask)
        sigchanid = find(stat_mask==1);
        SigTagSens{1,1} = OccipSens_full(stat_mask > 0);
        SigTagSens{1,2} = sigchanid;
    else % no tag sig sensors
        SigTagSens(1,1) = {'No'};
        SigTagSens{1,2} = 0;
        sigchanid = 1; %pre-select an occipital sensor for plotting later
    end
    eval(['TagCoh.' EpochType{ep} '_SigTagSens = SigTagSens;']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% 3. runing 2D coherence for different conditions %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TagCoh.freq = freqrange;
    timerange = epoch.time{1};
    TagCoh.time = timerange;
    TrlNum = []; % get trl numbers for each condition
    trlidx = {}; % get trl index for each condition
    %%%=== select the epoch data for each condition
    for ccc = 1:length(condnm)% conditon loop
        if ccc == 1
            m_trlid = epoch.trialinfo(:,tag_col) == loc2tar(ep); %trl idx for this EpochType
            idx = find(m_trlid);
        else
            idx = find(m_trlid & epoch.trialinfo(:,cond_col)==conds(ccc-1));
        end
        TrlNum = [TrlNum length(idx)];
        trlidx{ccc} = idx;
    end
    % make sure equal trial number across all conditions
    n_min = min(TrlNum(2:end));
    for ccc = 1:length(condnm)
        idx = trlidx{ccc};
        if ccc > 1 % data from conditions
            idx = idx(randperm(TrlNum(ccc)));
            idx = idx(1:n_min);
            trlidx{ccc} = idx;
        end
        cfg        = [];
        cfg.trials = idx;
        cfg.channel = [OccipSens_full;pds{tf}];
        eval(['data_' num2str(ccc) '= ft_selectdata(cfg, epoch);']);
    end
    % save out
    for ccc = 1:length(condnm)
        eval(['TagCoh.' EpochType{ep} '_TrlInfo{1,ccc} = data_' num2str(ccc) '.trialinfo;']);
        eval(['TagCoh.' EpochType{ep} '_RT(1,ccc) = nanmean(data_' num2str(ccc) '.trialinfo(:,fixdur_col));']);
        eval(['TagCoh.' EpochType{ep} '_TrlIdx_equ{1,ccc} = trlidx{ccc};']);
    end
    
    %%%=== get the 2D Coherence for all (grad)sensors
    for ccc = 1:length(condnm) % conditon loop
        %%% just to get a dummy coherence struct
        cfg            = [];
        cfg.output     = 'fourier';
        cfg.method     = 'mtmconvol';
        cfg.taper      = 'hanning';
        cfg.foi        = freqrange;
        cfg.toi        = timerange(1):0.5:timerange(end);
        cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
        cfg.keeptrials = 'yes';
        cfg.pad        = 'nextpow2';
        eval(['fourier = ft_freqanalysis(cfg, data_' num2str(ccc) ');']);
        % get coherence spctrm
        cfg            = [];
        cfg.method     = 'coh';
        cfg.channelcmb = {'all',pds{tf}}; %% all channell combinations calculated together
        coh            = ft_connectivityanalysis(cfg, fourier); %% the raw coherence
        coh.label      = fourier.label(1:length(coh.labelcmb));
        coh.time       = timerange;
        %%% do the real coherence using hilbert complex
        coh.cohspctrm = zeros(length(coh.label),length(coh.freq),length(coh.time));
        for fff = 1:length(freqrange) %freq loop
            % get hilbert filert data
            cfg = [];
            cfg.bpfilter   = 'yes';
            cfg.bpfreq     = [freqrange(fff)-freq_halfwidth freqrange(fff)+freq_halfwidth];
            cfg.hilbert    = 'complex';
            cfg.keeptrials = 'yes';
            eval(['fltdata = ft_preprocessing(cfg,data_' num2str(ccc) ');']);
            for chan = 1:length(fltdata.label)-1 %the last channel is photodiode channel
                for ttt = 1:length(fltdata.trial)
                    sig1(:,ttt) = fltdata.trial{ttt}(chan,:); %time*trl
                    sig2(:,ttt) = fltdata.trial{ttt}(end,:);  %the photodiode channel
                end
                spec1 = nanmean(sig1.*conj(sig1),2);%time*1
                spec2 = nanmean(sig2.*conj(sig2),2);
                specX = abs(nanmean(sig1.*conj(sig2),2)).^2;
                coh.cohspctrm(chan,fff,:) = specX./(spec1.*spec2);%time*1
            end
        end
        % get coh for the sig occipital tagging sensors
        cohcoh = nanmean(coh.cohspctrm(sigchanid,:,:),1);
        eval(['TagCoh.' EpochType{ep} '_Coh(:,:,1,ccc) = squeeze(cohcoh);']);
    end
end
%% %%%%%%%======== get the cohernece during the baseline ==========%%%
tf = 2;
cfg        = [];
cfg.channel = [OccipSens_full;pds{tf}];
epoch_BL_Cross_new = ft_selectdata(cfg,epoch_BL_Cross);
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = freqrange;
cfg.toi        = timerange(1):0.5:timerange(end);
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
cfg.keeptrials = 'yes';
cfg.pad        = 'nextpow2';
fourier = ft_freqanalysis(cfg,epoch_BL_Cross_new);
% get coherence spctrm
cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'all',pds{tf}}; %% all channell combinations calculated together
coh            = ft_connectivityanalysis(cfg, fourier); %% the raw coherence
coh.label      = fourier.label(1:length(coh.labelcmb));
coh.time       = timerange;
%%% do the real coherence using hilbert complex
coh.cohspctrm = zeros(length(coh.label),length(coh.freq),length(coh.time));
for fff = 1:length(freqrange) %freq loop
    % get hilbert filert data
    cfg = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [freqrange(fff)-freq_halfwidth freqrange(fff)+freq_halfwidth];
    cfg.hilbert    = 'complex';
    cfg.keeptrials = 'yes';
    fltdata = ft_preprocessing(cfg,epoch_BL_Cross);
    for chan = 1:length(fltdata.label)-1 %the last channel is photodiode channel
        for ttt = 1:length(fltdata.trial)
            sig1(:,ttt) = fltdata.trial{ttt}(chan,:); %time*trl
            sig2(:,ttt) = fltdata.trial{ttt}(end,:);  %the photodiode channel
        end
        spec1 = nanmean(sig1.*conj(sig1),2);%time*1
        spec2 = nanmean(sig2.*conj(sig2),2);
        specX = abs(nanmean(sig1.*conj(sig2),2)).^2;
        coh.cohspctrm(chan,fff,:) = specX./(spec1.*spec2);%time*1
    end
end
% get coh for the sig occipital tagging sensors
cohcoh = nanmean(coh.cohspctrm(sigchanid,:,:),1);
TagCoh.Baseline_Coh(:,:,1,1) = squeeze(cohcoh);

save([PPath.ResPath 'TagCoh_' num2str(sid)],'TagCoh','-v7.3');
end






