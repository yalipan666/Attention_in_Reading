%%% grup analysis of Coh_RIFT
server = 0; % run analysis on the server or local computer
%%% set paths
[rootdir,cmap] = Get_Paths(server);
PPath.ResPath = [rootdir 'Results' filesep 'Coh' filesep]; %result path
%%% running setting
RunCond = 'WrdOn';
freqrange = 40:1:80;
SenSelectP = 0.01;
EpochType = {'PreTarg';'Targ';'PosTarg'}; % epochs of interest
loc2tar = [-1 0 1]; % loc2tar for {'PreTarg','Targ','PosTarg'}
%%% get the tag freq for pretarget and target
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
%%% get the exp information
conds = ExpInfo.CondID;
condnm = [{'all'},ExpInfo.CondName];
cond_col = find(strcmp(ExpInfo.EventHdr,'SentenceCondition'));
tag_col = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
fixdur_col = find(strcmp(ExpInfo.EventHdr,'fixation_duration'));
firstpass_col = find(strcmp(ExpInfo.EventHdr,'FirstPassFix'));
subjects = ExpInfo.subjects;
nsub = length(subjects);
%%% set output
TagCoh_all = [];
TagCoh_all.CondName = condnm;
TagCoh_all.EventHdr = ExpInfo.EventHdr;
TagCoh_all.subs = subjects;
TagCoh_all.PreTarg_SigTagSens = cell(nsub,2);
TagCoh_all.Targ_SigTagSens = cell(nsub,2);
TagCoh_all.PosTarg_SigTagSens = cell(nsub,2);
TagCoh_all.Coh_hdr = {'freq*time*sub*cond'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run locally get the whole structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% basic parameters
%%% combine data together into one big structure
TagCoh_all.freq = 45:75; %Hz, shorten the freq range to put two versions together
% adjust the freq range,unify freq version to[60 65] for pretarget and target
freqid_v1 = {6:36; 6:36; 6:36;}; %[45:75]Hz for freq version 1
freqid_v2 = {11:41; 1:31; 11:41;}; %[45:75]Hz
for ps = 1:nsub
    load([PPath.ResPath 'TagCoh_' num2str(ps)]);
    v = num2str(ExpInfo.FreqVersion(ps));
    for fd = 1:length(EpochType)
        eval(['TagCoh_all.' EpochType{fd} '_TrlIdx_equ(ps,:)= TagCoh.' EpochType{fd} '_TrlIdx_equ;']);
        eval(['TagCoh_all.' EpochType{fd} '_RT(ps,:)= TagCoh.' EpochType{fd} '_RT;']);
        eval(['TagCoh_all.' EpochType{fd} '_TrlInfo(ps,:)= TagCoh.' EpochType{fd} '_TrlInfo;']);
        eval(['TagCoh_all.' EpochType{fd} '_SigTagSens(ps,:)= TagCoh.' EpochType{fd} '_SigTagSens;']);
        %%% adjust the freq range of coherence to be able to plot two
        %%% version together (both into [60 65]HZ tagging
        eval(['fff = freqid_v' v '{fd};']);
        eval(['TagCoh_all.' EpochType{fd} '_Coh(:,:,ps,:)= TagCoh.' EpochType{fd} '_Coh(fff,:,1,:);']);
    end
end
%%% get time
TagCoh_all.time = TagCoh.time;
TagCoh_all.freq_halfwidth = TagCoh.freq_halfwidth;
clear TagCoh
TagCoh = TagCoh_all;
clear TagCoh_all
TagCoh.SenSelectP = [num2str(SenSelectP) '_one tail'];
for ep = 1:length(EpochType)
    eval(['nosigsub = find(strcmp(TagCoh.' EpochType{ep} '_SigTagSens(:,1),''No''));']);
    TagCoh.SigSubID{1,ep} = transpose(setdiff(1:nsub,nosigsub));
    disp([EpochType{ep} ': SentVersion and FreqVersion'])
    % get sentence version
    tabulate(ExpInfo.SentVersion(TagCoh.SigSubID{1,ep})) % V1_15sub, V2_14subs
    % get frequency version
    tabulate(ExpInfo.FreqVersion(TagCoh.SigSubID{1,ep})) % V1_15sub, V2_14subs
end
save([PPath.ResPath 'TagCoh.mat'],'TagCoh','-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. group plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting-coh
% % % plot two frequency versions combined
figtitle = 'Group_Coh';
tagfreq = [60 65]; %use freq version 1 as the standard
sigsubid = TagCoh.SigSubID{1,1}; % use the pre-target sig subjects for all epochs
whichfreq = [1 2 1]; % analyze which tag freq for {'PreTarg','Targ','PosTarg'}
% plot
h = figure('Name',figtitle,'color',[1 1 1]);
ncol = length(EpochType);
nrow = length(condnm)+1;
for ep = 1:ncol
    tf = whichfreq(ep);
    %%%% coh plot
    eval(['tmpdata = nanmean(TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,:),3);']);%%% freq*time*sub*cond
    tmpcc = [1 2 3];
    for cc = 1:nrow
        if ismember(cc,tmpcc)
            loc = ep + ncol*(cc-1);
            plotdata = squeeze(tmpdata(:,:,:,cc));
            subtitle = [EpochType{ep} 'et,' TagCoh.CondName{cc}];
            subplot(nrow,ncol,loc);
            pcolor(TagCoh.time, TagCoh.freq, plotdata); colorbar; axcopy;
            xlim([-0.3 0.3])
            if ep == 2 %tag at 65Hz
                ylim([55 75])
                set(gca,'YTick',55:5:75);
                set(gca,'YTickLabel',{'f2-10','f2-5',' f2 ','f2+5','f2+10'},'FontWeight','bold','FontSize',10);
            else
                ylim([50 70])
                set(gca,'YTick',50:5:70);
                set(gca,'YTickLabel',{'f1-10','f1-5',' f1 ','f1+5','f1+10'},'FontWeight','bold','FontSize',10);
            end
            if cc == 1
                caxis([0 0.015])
            else
                caxis([0 0.01])
            end
        else %get the contrast between conditions
            loc = ep + ncol*(cc-1);
            idx1 = 2; idx2 = 3;
            plotdata = squeeze(tmpdata(:,:,:,idx1)-tmpdata(:,:,:,idx2));
            subtitle = [EpochType{ep} 'et Coh-contrast'];
            subplot(nrow,ncol,loc);
            pcolor(TagCoh.time, TagCoh.freq, plotdata); colorbar; axcopy;
            caxis([-max(abs(caxis)) max(abs(caxis))])
            xlim([-0.3 0.3])
            if ep == 2 %tag at 65Hz
                ylim([55 75])
                set(gca,'YTick',55:5:75);
                set(gca,'YTickLabel',{'f2-10','f2-5',' f2 ','f2+5','f2+10'},'FontWeight','bold','FontSize',10);
            else
                ylim([50 70])
                set(gca,'YTick',50:5:70);
                set(gca,'YTickLabel',{'f1-10','f1-5',' f1 ','f1+5','f1+10'},'FontWeight','bold','FontSize',10);
            end
            caxis([-0.002 0.002])
        end
        colormap(cmap);shading interp
        hold on;
        plot([-1 1],[tagfreq(tf)-5 tagfreq(tf)-5],'-.k','LineWidth',2)
        plot([-1 1],[tagfreq(tf) tagfreq(tf)],'-.k','LineWidth',2)
        plot([-1 1],[tagfreq(tf)+5 tagfreq(tf)+5],'-.k','LineWidth',2)
        plot([-0.2 -0.2],[50 75],'-.k','LineWidth',2)
        plot([0 0],[50 75],'-.k','LineWidth',2)
        plot([0.2 0.2],[50 75],'-.k','LineWidth',2)
        set(gca,'XTick',-0.4:0.2:0.4);
        set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'},'FontWeight','bold','FontSize',10)
        title(subtitle,'FontWeight','bold','FontSize',14);
    end
end
ScSz = [0 0 1300 800];
set(gcf,'Position',ScSz);
saveas(h,[PPath.ResPath figtitle]);
set(h, 'renderer', 'painters')
saveas(h,['U:\writing\Full_Attention\' figtitle '.svg']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%============ permutation on the group level ============ %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	Group level-->  only shuffle conditions between subjects, keep time and
% 	freq! (cluster based!)
%%% get the data
stat_tagfreq = -3:3; % thefreq range used to average coherence(centre freq indexes as 0)
for ep = 1:length(EpochType)
    tf = whichfreq(ep);
    sigsubid = TagCoh.SigSubID{1,ep};
    % construct other stuff
    tmp = [];
    tmp.label = {'MEG1922'}; %pick a random label-> the actual data are the averaged coh over all tagsig_sensors
    tmp.freq = TagCoh.freq;
    tmp.time = TagCoh.time;
    tmp.dimord = 'chan_freq_time_subj';
    for c = 1:2 %condition loop
        eval(['powspctrm = TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,1+c);']);%%% freq*time*sub
        % add the chan dimention
        tmp.powspctrm = nan([1,size(powspctrm)]);
        tmp.powspctrm(1,:,:,:) = powspctrm;
        eval(['coh_' num2str(c) ' = tmp;']);
    end
    % do the statistic
    cfg = [];
    cfg.channel          = tmp.label;
    cfg.latency          = [0 0.2];
    cfg.frequency        = tagfreq(tf)+stat_tagfreq([1 end]);
    cfg.parameter        = 'powspctrm';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 5000;
    % construct the design matrix
    n_sub = length(sigsubid);
    design = zeros(2,2*n_sub);
    design(1,:) = [1:n_sub 1:n_sub];
    design(2,1:n_sub)        = 1;
    design(2,n_sub+1:2*n_sub) = 2;
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    stat = ft_freqstatistics(cfg, coh_1, coh_2);
    any(any(any(stat.mask)))
    % get the mask
    id_time = arrayfun(@(x) dsearchn(TagCoh.time',x),stat.time);
    id_freq = arrayfun(@(x) dsearchn(TagCoh.freq',x),stat.freq);
    mask_full = zeros(length(TagCoh.freq),length(TagCoh.time));
    if ep == 1 %%only found sig positive clusters
        id_mask = squeeze(stat.posclusterslabelmat==1);
    elseif ep == 2 %%only found sig negative clusters
        id_mask = squeeze(stat.negclusterslabelmat==1);
    else %% no sig clusters found
        id_mask = 0;
    end
    mask_full(id_freq,id_time)= id_mask;
    %%% save out
    stat.mask_full = mask_full;
    eval(['TagCoh.' EpochType{ep} '_stat = stat;']);
end
%%%% ==== stat for pre-pre-target  ====%%%%%
ep = 1;
tf = whichfreq(ep);
sigsubid = TagCoh.SigSubID{1,ep};
% construct other stuff
tmp = [];
tmp.label = {'MEG1922'}; %pick a random label-> the actual data are the averaged coh over all tagsig_sensors
tmp.freq = TagCoh.freq;
tmp.time = TagCoh.time;
tmp.dimord = 'chan_freq_time_subj';
for c = 1:2 %condition loop
    eval(['powspctrm = TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,1+c);']);%%% freq*time*sub
    % add the chan dimention
    tmp.powspctrm = nan([1,size(powspctrm)]);
    tmp.powspctrm(1,:,:,:) = powspctrm;
    eval(['coh_' num2str(c) ' = tmp;']);
end
% do the statistic
cfg = [];
cfg.channel          = tmp.label;
cfg.latency          = [-0.2 0];
cfg.frequency        = tagfreq(tf)+stat_tagfreq([1 end]);
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
% construct the design matrix
n_sub = length(sigsubid);
design = zeros(2,2*n_sub);
design(1,:) = [1:n_sub 1:n_sub];
design(2,1:n_sub)        = 1;
design(2,n_sub+1:2*n_sub) = 2;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;
stat = ft_freqstatistics(cfg, coh_1, coh_2);
any(any(any(stat.mask)))
% get the mask
id_time = arrayfun(@(x) dsearchn(TagCoh.time',x),stat.time);
id_freq = arrayfun(@(x) dsearchn(TagCoh.freq',x),stat.freq);
mask_full = zeros(length(TagCoh.freq),length(TagCoh.time));

id_mask = squeeze(stat.posclusterslabelmat==1);

mask_full(id_freq,id_time)= id_mask;
%%% save out
stat.mask_full = mask_full;
TagCoh.PrepreTarg_stat = stat;

%%% === save out the results ==== %%%
save([PPath.ResPath 'TagCoh.mat'],'TagCoh','-v7.3');

%% ====== plot the significant conditional contrast of coherence
%%% pre-target (show both the pre-target and the pre-pre-target
opq = cell(1,3); % mask for all epochs
opq{1} = TagCoh.PreTarg_stat.mask_full;
opq{2} = TagCoh.Targ_stat.mask_full;
opq{3} = TagCoh.PosTarg_stat.mask_full;
% mask out the non-sig area
opq{1}(opq{1}==0)=0.4;
opq{2}(opq{2}==0)=0.4;
opq{3}(opq{3}==0)=0.4;

figtitle = 'Coh_contrast_group';
h = figure('Name',figtitle,'color',[1 1 1]);
for ep = 1:length(EpochType)
    tf = whichfreq(ep);
    sigsubid = TagCoh.SigSubID{1,ep};
    eval(['coh_diff = TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,2)-TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,3);']); %low-high
    coh_diff = nanmean(coh_diff,3);
    % plot
    subplot(1,length(EpochType),ep)
    pcolor(TagCoh.time, TagCoh.freq, coh_diff.*opq{ep}); colorbar;
    caxis([-max(abs(caxis)) max(abs(caxis))])
    colormap(cmap);shading interp
    hold on;
    % add highlight outline
    [x,y] = meshgrid(TagCoh.time, TagCoh.freq);
    x = interp2(x, 2); % change to 4 for round corners
    y = interp2(y, 2); % change to 4 for round corners
    contourlines = opq{ep} ==1;
    contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
    dx = mean(diff(x(1, :))); % remove for round corners
    dy = mean(diff(y(:, 1))); % remove for round corners
    contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[0 0 0],'LineWidth',2);
    % add time and freq lines
    plot([-1 1],[tagfreq(tf)-5 tagfreq(tf)-5],'-.k','LineWidth',2)
    plot([-1 1],[tagfreq(tf) tagfreq(tf)],'-.k','LineWidth',2)
    plot([-1 1],[tagfreq(tf)+5 tagfreq(tf)+5],'-.k','LineWidth',2)
    if ep == 2 %tag at 65Hz
        ylim([55 75])
        set(gca,'YTick',55:5:75);
        set(gca,'YTickLabel',{'f2-10','f2-5',' f2 ','f2+5','f2+10'},'FontWeight','bold','FontSize',10);
    else
        ylim([50 70])
        set(gca,'YTick',50:5:70);
        set(gca,'YTickLabel',{'f1-10','f1-5',' f1 ','f1+5','f1+10'},'FontWeight','bold','FontSize',10);
    end
    plot([-0.2 -0.2],[50 75],'-.k','LineWidth',1)
    plot([0 0],[50 75],'-.k','LineWidth',1)
    plot([0.2 0.2],[50 75],'-.k','LineWidth',1)
    xlim([-0.3 0.3])
    set(gca,'XTick',-0.4:0.2:0.4);
    set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'},'FontWeight','bold','FontSize',10)
    title([EpochType{ep} 'et coh contrast'],'FontWeight','bold','FontSize',12);
    caxis([-0.002 0.002])
end
ScSz = [100 100 1400 300];
set(gcf,'Position',ScSz);
saveas(h,[PPath.ResPath figtitle]);
set(h, 'renderer', 'painters')
saveas(h,['U:\writing\Full_Attention\' figtitle '.svg']);


%%% === plot the pre-target epochs stat during both [-0.2 0]s and [0 0.2]s
opq{1} = double(TagCoh.PreTarg_stat.mask_full + TagCoh.PrepreTarg_stat.mask_full > 0);
opq{1}(opq{1}==0)=0.4;
figtitle = 'Coh_contrast_group_PrepreTarget';
h = figure('Name',figtitle,'color',[1 1 1]);
ep = 1;
tf = whichfreq(ep);
sigsubid = TagCoh.SigSubID{1,ep};
eval(['coh_diff = TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,2)-TagCoh.' EpochType{ep} '_Coh(:,:,sigsubid,3);']); %low-high
coh_diff = nanmean(coh_diff,3);
% plot
pcolor(TagCoh.time, TagCoh.freq, coh_diff.*opq{ep}); colorbar;
caxis([-max(abs(caxis)) max(abs(caxis))])
colormap(cmap);shading interp
hold on;
% add highlight outline
[x,y] = meshgrid(TagCoh.time, TagCoh.freq);
x = interp2(x, 2); % change to 4 for round corners
y = interp2(y, 2); % change to 4 for round corners
contourlines = opq{ep} ==1;
contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners
dx = mean(diff(x(1, :))); % remove for round corners
dy = mean(diff(y(:, 1))); % remove for round corners
contour(x+dx/2,y+dy/2,contourlines,1,'EdgeColor',[0 0 0],'LineWidth',2);
% add time and freq lines
plot([-1 1],[tagfreq(tf)-5 tagfreq(tf)-5],'-.k','LineWidth',2)
plot([-1 1],[tagfreq(tf) tagfreq(tf)],'-.k','LineWidth',2)
plot([-1 1],[tagfreq(tf)+5 tagfreq(tf)+5],'-.k','LineWidth',2)
ylim([50 70])
set(gca,'YTick',50:5:70);
set(gca,'YTickLabel',{'f1-10','f1-5',' f1 ','f1+5','f1+10'},'FontWeight','bold','FontSize',10);
plot([-0.2 -0.2],[50 75],'-.k','LineWidth',1)
plot([0 0],[50 75],'-.k','LineWidth',1)
plot([0.2 0.2],[50 75],'-.k','LineWidth',1)
xlim([-0.3 0.3])
set(gca,'XTick',-0.4:0.2:0.4);
set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'},'FontWeight','bold','FontSize',10)
title(figtitle,'FontWeight','bold','FontSize',12, 'Interpreter', 'none');
caxis([-0.002 0.002])
ScSz = [100 100 400 300];
set(gcf,'Position',ScSz);
saveas(h,[PPath.ResPath figtitle]);
set(h, 'renderer', 'painters')
saveas(h,['U:\writing\Full_Attention\' figtitle '.svg']);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Get total topograph of RFT %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get the occipital sensors
load([rootdir filesep 'Analyse_data' filesep TagCoh.subs{1} filesep 'epoch_BL_Cross.mat']) %%% just load random data that has all the labels
cfg = [];
erf = ft_timelockanalysis(cfg, epoch_BL_Cross);
for ep = 1:length(EpochType)
    %%% get labels
    figtitle = [EpochType{ep} ' RIFT sensors'];
    eval(['SigSen = TagCoh.' EpochType{ep} '_SigTagSens;']);
    tmp = zeros(size(erf.label));
    for bbb = 1:length(SigSen)
        tmpsen = SigSen{bbb};
        if ~strcmp(tmpsen,'No')
            for kkk = 1:length(tmpsen)
                idx = find(strcmp(erf.label,[tmpsen{kkk}(1:end-1) '1']));
                tmp(idx) = tmp(idx)+1;
            end
        end
    end
    erf.avg = repmat(tmp,1,size(erf.avg,2));
    h = figure('color', [1 1 1],'name',figtitle);
    cfg =[];
    cfg.layout = 'neuromag306mag.lay'; %'neuromag306cmb.lay';
    cfg.comment = ' ';
    ft_topoplotER(cfg,erf);
    colorbar('FontWeight','bold','FontSize',10);
    caxis([0 max(tmp)])
    colormap(cmap);
    title(figtitle,'FontWeight','bold','FontSize',16);
    text(0.65,0.3,'No. of sensors','FontWeight','bold','FontSize',14,'Rotation',270)
    hold on; plot([0 0],[-0.5 0.5],'-.w','LineWidth',2)
    saveas(h,[PPath.ResPath figtitle]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% coherence latency statistics ---- JackKnife %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get the cohenrece curve over the sig freq for epochs that are sig in
%%% the permutation test
%%% (also get the averaged cohernece over the sig freq and sig time in clusters, can be used for simple Ttest)
Epoch_Sig = {'PreTarg';'Targ'}; 
for ep = 1:length(Epoch_Sig)
    % get the freq range from the significant cluster from the gourp
    % permutation test
    eval(['mask_full = TagCoh.' Epoch_Sig{ep} '_stat.mask_full;']);
    % get the sig time_freq mask from TagCoh
    if ep == 3 || ep == 4 %postarget is not sig, mask_full = 0
        fid = dsearchn(TagCoh.freq',60)+stat_tagfreq;
        tid = dsearchn(TagCoh.time',0):dsearchn(TagCoh.time',200);
    else
        [fid,tid] = find(mask_full == 1);
        fid = unique(fid);
        tid = unique(tid);
        sigtim_range = [TagCoh.time(tid(1)) TagCoh.time(tid(end))];
        if ep == 1 %pre-target
            freq_rangeWRTtagfreq = [TagCoh.freq(fid(1)) TagCoh.freq(fid(end))]-60;
        elseif ep == 2 %target
            freq_rangeWRTtagfreq = [TagCoh.freq(fid(1)) TagCoh.freq(fid(end))]-65;
        end
        eval(['TagCoh.' Epoch_Sig{ep} '_sigtim_range = sigtim_range;']);
        eval(['TagCoh.' Epoch_Sig{ep} '_freq_rangeWRTtagfreq = freq_rangeWRTtagfreq;']);
    end
    sigsubid = TagCoh.SigSubID{1,ep};
    %%% get the coherence curve data over the sig freq
    eval(['tempdata = TagCoh.' Epoch_Sig{ep} '_Coh(fid,:,sigsubid,:);']);
    tempdata = squeeze(nanmean(tempdata,1)); %%% time*sub*cond
    eval(['TagCoh.' Epoch_Sig{ep} '_Coh_TimSubCond = tempdata;']);
    %%% get the averaged coherence oever the sig time-freq cluster
    data4Ttest = squeeze(nanmean(tempdata(tid,:,[2 3]),1));
    eval(['TagCoh.' Epoch_Sig{ep} '_stat.data4Ttest = data4Ttest;']);
end
%%% run jackknife and do statistic over the onset latency of the cohernece 
%%% only look at the time window between [0 rt]s
for ep = 1:length(Epoch_Sig)
% % % % %     %%% for plot
% % % % %     h = figure('Name',[Epoch_Sig{ep} '--Coh Curve during jackknife'],'color',[1 1 1],'Position',[100 100 1800 900]);
    sigsubid = TagCoh.SigSubID{1,ep};
    eval(['cohdata_all = TagCoh.' Epoch_Sig{ep} '_Coh_TimSubCond(:,:,[2 3]);']); %% time*sub*cond
    eval(['RT_all = TagCoh.' Epoch_Sig{ep} '_RT(sigsubid,[2 3])./1000;']); %% sub*cond
    HM_tim = []; %% when reach the half maximum amplitude (distant amplitude between max and min coh)  
    %%% doing jackknife loop
    allsigsub = 1:length(sigsubid);
    exclu_sub = [allsigsub 0]; %%[jackknife all]
    for sid = 1:length(exclu_sub)
        subid_jn = setdiff(allsigsub,exclu_sub(sid));
        %%% data after jackknife
        cohdata = cohdata_all(:,subid_jn,:);
        RT = RT_all(subid_jn,:);
        meanRFT = squeeze(nanmean(cohdata,2)); %%% time*cond
        seRFT = squeeze(nanstd(cohdata,0,2)./sqrt(size(cohdata,2)));
        meanRT = mean(RT);
        %%% get the time window (onset should be within the fixation duration time window)
        %%% note the time window doesn't affect the result, if choose [0 0.3]
        %%% will get the same results
        rt = min(meanRT);
        etp = find(TagCoh.time < rt); %%% rt window aligned with zero time point
        zero_tp = dsearchn(TagCoh.time',0);
        tprange = zero_tp:etp(end);
        timrange = TagCoh.time(tprange);
        for cc = 1:2 %%[Low High] condition loops
            %%%=== half amplitude latency
            tmptp = meanRFT(tprange,cc); %%% data in the timewindow: time*sub*cond
            [mintmp, minidx] = min(tmptp);
            %%% make sure that the maxtmp is later than the mintmp
            tmptp = tmptp(minidx:end);
            [maxtmp, maxidx] = max(tmptp);
            tmpidx = find(tmptp>((maxtmp-mintmp)/2+mintmp));
            if ~isempty(tmpidx)
                HM_tim(sid,cc) = timrange(tmpidx(1)+minidx);
            else
                HM_tim(sid,cc) = nan;
            end
        end
% % % % %         %%% plot each curve during the loop to double-check the onset
% % % % %         %%% latency estimation
% % % % %         subplot(5,7,sid)
% % % % %         a = shadedErrorBar(TagCoh.time,meanRFT(:,1),seRFT(:,1),{'color',colmat(1,:)},0.8);hold on;
% % % % %         b = shadedErrorBar(TagCoh.time,meanRFT(:,2),seRFT(:,2),{'color',colmat(2,:)},0.9);
% % % % %         xlim([0 0.4])
% % % % %         timeliney = get(gca,'ylim');
% % % % %         timelinex = repmat(HM_tim(sid,1),1,2);
% % % % %         plot(timelinex', timeliney','--','color',colmat(1,:),'LineWidth',1);
% % % % %         timelinex = repmat(HM_tim(sid,2),1,2);
% % % % %         plot(timelinex', timeliney','--','color',colmat(2,:),'LineWidth',1);
    end
    eval(['JkNf.' Epoch_Sig{ep} '_HM_tim = HM_tim ;']);
end
%%% using jacknife to do the statistics
tp_tpye = 'HM_tim';
N = length(sigsubid);
for ep = 1:length(Epoch_Sig)
    tj = [];
    eval(['tmptim = JkNf.' Epoch_Sig{ep} '_' tp_tpye ';']);
    Di = tmptim(1:end-1,1)-tmptim(1:end-1,2);
    J = nansum(Di)/N;
    SD = sqrt([(N-1)/N]*[nansum((Di-J).^2)]);
    D = tmptim(end,1) - tmptim(end,2);
    tj = D/SD;
    eval(['JkNf.'  Epoch_Sig{ep} '_' tp_tpye '_tvalue = tj;']);
end
% check if there's any sig
% lookup t-value table https://www.socscistatistics.com/pvalues/tdistribution.aspx
JkNf.PreTarg_HM_tim_pvalue = '0.147_one tail';
JkNf.Targ_HM_tim_pvalue = '0.243_one tail';
TagCoh.Jackknife_latency = JkNf;
save([PPath.ResPath 'TagCoh.mat'], 'TagCoh')

%% ========== plot the group level curve for Epoch_Sig
curve_xmax = 0.4; %%% the x-axis range
etp = find(TagCoh.time < curve_xmax); %%% rt window aligned with zero--saccadeonset
zero_tp = nearest(TagCoh.time,0);
tprange = zero_tp:etp(end);
timrange = TagCoh.time(tprange);
colmat = [0 114 189;217 83 25]./255;
figtitle = 'Coh Curve over TagSigCluster (not sig in latency)';
h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 1050 260]);
for ep = 1:length(Epoch_Sig)
    sigsubid = TagCoh.SigSubID{1,ep};
    eval(['cohdata_all = TagCoh.' Epoch_Sig{ep} '_Coh_TimSubCond(tprange,:,[2 3]);']); %% time*sub*cond
    meanRFT = squeeze(mean(cohdata_all,2));
    seRFT = squeeze(nanstd(cohdata_all,0,2)./sqrt(length(sigsubid)));
    eval(['HM = TagCoh.Jackknife_latency.' Epoch_Sig{ep} '_HM_tim(end,:,1);']);
    %%%% plot
    subplot(1,length(Epoch_Sig),ep)
    a = shadedErrorBar(timrange,meanRFT(:,1),seRFT(:,1),{'color',colmat(1,:)},0.8);hold on;
    b = shadedErrorBar(timrange,meanRFT(:,2),seRFT(:,2),{'color',colmat(2,:)},0.9);
    a.mainLine.LineWidth = 1;
    b.mainLine.LineWidth = 1;
    legendflex([a.mainLine,b.mainLine],{'Low';'High'},'anchor', {'ne','ne'}, 'buffer', [0 0],'Fontsize',7,'xscale',1,'box','off');
    %     a = plot(timrange,meanRFT(:,1),'color',colmat(1,:),'LineWidth',2);hold on;
    %     b = plot(timrange,meanRFT(:,2),'color',colmat(2,:),'LineWidth',2);
    %     legendflex([a,b],{'Low';'High'},'anchor', {'ne','ne'}, 'buffer', [0 0],'FontWeight','bold','Fontsize',8,'xscale',1,'box','off');
    set(gca,'FontSize',8,'FontWeight','bold');
    set(gca,'box','off')
    set(gca,'XTick',0:0.1:0.4);
    set(gca,'XTickLabel',{'FixOn','0.1','0.2','0.3','0.4'},'FontWeight','bold','FontSize',8)
    %     ylim([0.004 0.01])
    title([Epoch_Sig{ep} 'et at tagfreq #' num2str(ep)],'FontWeight','bold','FontSize',10)
    xlabel('Time (s)','FontWeight','bold','FontSize',10)
    ylabel('Coherence (r^2)','FontWeight','bold','FontSize',10)
    %%% plot the latency lines
    hold on;
    timeliney = get(gca,'ylim');
    timelinex = repmat(HM(1),1,2);
    plot(timelinex', timeliney','--','color',colmat(1,:),'LineWidth',1);
    timelinex = repmat(HM(2),1,2);
    plot(timelinex', timeliney','--','color',colmat(2,:),'LineWidth',1);
end
saveas(h,[PPath.ResPath figtitle]);
































 