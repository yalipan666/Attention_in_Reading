%%% grup analysis of Coh_RIFT
server = 0; % run analysis on the server or local computer
%%% set paths
[rootdir,cmap] = Get_Paths(server);
PPath.ResPath = [rootdir 'Results' filesep 'Coh_parallel' filesep]; %result path
%%% running setting
RunCond = 'WrdOn';
freqrange = 40:1:80;
SenSelectP = 0.01;
EpochType = {'PreTarg','Targ','PosTarg'};
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
freqid_v2 = {1:31; 11:41; 1:31;}; %[45:75]Hz
for ps = 1:nsub
    load([PPath.ResPath 'TagCoh_' num2str(ps)]);
    v = num2str(ExpInfo.FreqVersion(ps));
    for ep = 1:length(EpochType)
        eval(['TagCoh_all.' EpochType{ep} '_TrlIdx_equ(ps,:)= TagCoh.' EpochType{ep} '_TrlIdx_equ;']);
        eval(['TagCoh_all.' EpochType{ep} '_RT(ps,:)= TagCoh.' EpochType{ep} '_RT;']);
        eval(['TagCoh_all.' EpochType{ep} '_TrlInfo(ps,:)= TagCoh.' EpochType{ep} '_TrlInfo;']);
        eval(['TagCoh_all.' EpochType{ep} '_SigTagSens(ps,:)= TagCoh.' EpochType{ep} '_SigTagSens;']);
        %%% adjust the freq range of coherence to be able to plot two
        %%% version together (both into [60 65]HZ tagging
        eval(['fff = freqid_v' v '{ep};']);
        eval(['TagCoh_all.' EpochType{ep} '_Coh(:,:,ps,:)= TagCoh.' EpochType{ep} '_Coh(fff,:,1,:);']);
    end
    TagCoh_all.Baseline_Coh(ps,1) = mean(mean(TagCoh.Baseline_Coh(fff,:)));
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
    tabulate(ExpInfo.SentVersion(TagCoh.SigSubID{1,ep}))  
    % get frequency version
    tabulate(ExpInfo.FreqVersion(TagCoh.SigSubID{1,ep}))  
end
save([PPath.ResPath 'TagCoh_Parallel.mat'],'TagCoh','-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. group plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting-coh
% % % plot two frequency versions combined
figtitle = 'Group_Coh_Foveal';
tagfreq = [60 65]; %use freq version 1 as the standard
sigsubid = TagCoh.SigSubID{1,1}; % use the pre-target sig subjects for all epochs
% plot
whichfreq = [2 1 2]; 
h = figure('Name',figtitle,'color',[1 1 1]);
nrow = length(condnm)+1;
ncol = length(EpochType);
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
            if tf == 2 %tag at 65Hz
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
            if tf == 2 %tag at 65Hz
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
ScSz = [0 0 650 800];
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
    any(any(any(stat.mask))) %% no sig for pre-target, target and post-target
    %%% save out
    eval(['TagCoh.' EpochType{ep} '_stat = stat;']);
end
save([PPath.ResPath 'TagCoh_Parallel.mat'],'TagCoh','-v7.3');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Get total topograph of RFT %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get the occipital sensors
load([rootdir filesep 'Analyse_data' filesep TagCoh.subs{1} filesep 'epoch_BL_Cross.mat']) %%% just load random data that has all the labels
cfg = [];
erf = ft_timelockanalysis(cfg, epoch_BL_Cross);
for ep = 1:length(EpochType)
    %%% get labels
    figtitle = [EpochType{ep} 'Foveal RIFT sensors'];
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% plot coh curve for parafoveal + foveal %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the freq id, all in freq_version_v1
freqid_f1 = 13:19; %[57 63]Hz, target flicker at freq1 60Hz
freqid_f2 = 18:24; %[62 68]Hz, post-target flicker at freq2 65Hz
%%% get the foveal coherence
load('Z:\Full_Attention\Results\Coh_parallel\TagCoh_Parallel.mat')
TagCoh_fov = TagCoh;
%%% for pre-target epoch, F1 is para freq, F2 is fovea freq
PreTarg_F2 = squeeze(nanmean(TagCoh_fov.PreTarg_Coh(freqid_f2,:,TagCoh_fov.SigSubID{1,1},:),1));%time*sub*cond
Targ_F1 = squeeze(nanmean(TagCoh_fov.Targ_Coh(freqid_f1,:,TagCoh_fov.SigSubID{1,2},:),1));
PosTarg_F2 = squeeze(nanmean(TagCoh_fov.PosTarg_Coh(freqid_f2,:,TagCoh_fov.SigSubID{1,3},:),1));
% get the baseline coherence
Baseline = TagCoh_fov.Baseline_Coh(TagCoh_fov.SigSubID{1,2},1);
nnn = length(TagCoh_fov.time);
mean_bl = repmat(mean(Baseline),nnn,1);
se_bl = repmat(std(Baseline,0,1)./sqrt(length(Baseline)),nnn,1);

%%% get the parafoveal coherence
load('Z:\Full_Attention\Results\Coh\TagCoh.mat')
PreTarg_F1 = squeeze(nanmean(TagCoh.PreTarg_Coh(freqid_f1,:,TagCoh.SigSubID{1,1},:),1)); %time*sub*cond
Targ_F2 = squeeze(nanmean(TagCoh.Targ_Coh(freqid_f2,:,TagCoh.SigSubID{1,2},:),1));
PosTarg_F1 = squeeze(nanmean(TagCoh.PosTarg_Coh(freqid_f1,:,TagCoh.SigSubID{1,3},:),1));

%%% plot
colmat = [0 114 189;217 83 25]./255; %[blue-fovea, orange-para]
figtitle = 'Coh_Curve_Parallel';
h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 1200 600]);
nrow = length(condnm);
ncol = length(EpochType);
for cc = 1:nrow
    for ep = 1:ncol
        eval(['Coh_F1 =' EpochType{ep} '_F1(:,:,cc);']); %% time*sub
        eval(['Coh_F2 =' EpochType{ep} '_F2(:,:,cc);']); %% time*sub
        se_para = nanstd(Coh_F2,0,2)./sqrt(size(Coh_F2,2));
        se_fovea = nanstd(Coh_F1,0,2)./sqrt(size(Coh_F1,2));
        %%%% plot
        subplot(nrow,ncol,ep+ncol*(cc-1))
        a = shadedErrorBar(TagCoh.time,mean(Coh_F1,2),se_fovea,{'color',colmat(1,:)},0.8);hold on;
        b = shadedErrorBar(TagCoh.time,mean(Coh_F2,2),se_para,{'color',colmat(2,:)},0.9);
        c = shadedErrorBar(TagCoh.time,mean_bl,se_bl,{'color',[.5 .5 .5]},0.9);
        a.mainLine.LineWidth = 1;
        b.mainLine.LineWidth = 1;
        c.mainLine.LineStyle = '-.';
        legendflex([a.mainLine,b.mainLine],{'TagFreq_F1';'TagFreq_F2'},...
            'anchor', {'ne','ne'}, 'buffer', [0 0],'Fontsize',7,'xscale',1,'box','off','Interpreter', 'none');
        set(gca,'FontSize',7);
        set(gca,'box','off')
        xlim([-0.45 0.45])
        set(gca,'XTick',-0.4:0.2:0.4);
        set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'},'FontSize',7)
        title([EpochType{ep} 'et, ' condnm{cc}],'FontSize',7)
        xlabel('Time (s)','FontSize',7)
        ylabel('Coherence (r^2)','FontSize',7)
        %%% plot the latency lines
        timeliney = get(gca,'ylim');
        hold on;
        plot([0;0], timeliney','--','color',[0 0 0],'LineWidth',1);
        plot([0.2;0.2], timeliney','--','color',[0 0 0],'LineWidth',1);
    end
end 
saveas(h,[PPath.ResPath figtitle]);
set(gcf, 'renderer', 'painters')
saveas(h,['U:\writing\Full_Attention\' figtitle '.svg']);


%%% run t-test to compare coh during baseline VS Coh_F1 VS Coh_F2 during
%%% target fixations --> only select overlapped subjects between Coh_F2&Coh_F1
sub_inboth = intersect(TagCoh_fov.SigSubID{1,2},TagCoh.SigSubID{1,2});
freqid_f1 = 13:19; %[57 63]Hz, target flicker at freq1 60Hz
freqid_f2 = 18:24; %[62 68]Hz, post-target flicker at freq2 65Hz
tim_id = 501:700;%[0 0.2]s
cond_id = 1;
% get the overall foveal coherence 
Targ_fovea_all = squeeze(nanmean(nanmean(TagCoh_fov.Targ_Coh(freqid_f1,tim_id,sub_inboth,cond_id),1),2));
% get the baseline coherence
Baseline = TagCoh_fov.Baseline_Coh(sub_inboth,1);
% get the overall parafoveal coherence 
Targ_para_all = squeeze(nanmean(nanmean(TagCoh.Targ_Coh(freqid_f2,tim_id,sub_inboth,cond_id),1),2));
 
% t-test
Stat_Target.Targ_fovea_all = Targ_fovea_all;
Stat_Target.Targ_para_all = Targ_para_all;
Stat_Target.Baseline = Baseline;
[~,p,~,stats] = ttest(Targ_fovea_all, Baseline);stats.p = p
Stat_Target.foveaVSbaseline = stats;
[~,p,~,stats] = ttest(Targ_para_all, Baseline);stats.p = p
Stat_Target.paraVSbaseline = stats;
[~,p,~,stats] = ttest(Targ_para_all, Targ_fovea_all);stats.p = p
Stat_Target.paraVfovea = stats;
save('Z:\Full_Attention\Results\Coh_parallel\Stat_Target.mat','Stat_Target')

%%% plot
colmat = [0 114 189;217 83 25]./255; %[blue-Coh_F1, orange-Coh_F2]
figtitle = 'Ttest coh overall';
nsigsub = length(Targ_fovea_all);
group = [cellstr(repmat('Baseline',nsigsub,1));cellstr(repmat('Fovea',nsigsub,1));cellstr(repmat('Parafovea',nsigsub,1))];
grouporder={'Baseline','Fovea','Parafovea'};
h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 180 180]);
vdata = [Baseline;Targ_fovea_all;Targ_para_all];
vp = violinplot(vdata, group,'GroupOrder',grouporder);
vp(1).ViolinColor = [.5 .5 .5];
vp(2).ViolinColor = colmat(1,:);
vp(3).ViolinColor = colmat(2,:);
vp(1).ShowMean = 1; vp(2).ShowMean = 1;vp(3).ShowMean = 1;
vp(1,1).MedianPlot.Visible = 'off';
vp(1,2).MedianPlot.Visible = 'off';
vp(1,3).MedianPlot.Visible = 'off';
vp(1,1).MeanPlot.LineWidth = 1.5;
vp(1,2).MeanPlot.LineWidth = 1.5;
vp(1,3).MeanPlot.LineWidth = 1.5;
vp(1,1).BoxWidth = 0.01; vp(1,2).BoxWidth = 0.01;vp(1,3).BoxWidth = 0.01;
ylabel('Coherence (r^2)','FontSize',7,'FontWeight','normal','FontName','Arial');
set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
set(gca,'box','off','LineWidth',1)
xlim([0.5 3.5])
ylim([0 0.03])
%%% add stat
plot([1 2],[0.025 0.025],'k','LineWidth',1)
text(1.25,0.0255,'***','FontWeight','normal','FontSize',14,'FontName','Arial')
plot([2 3],[0.026 0.026],'k','LineWidth',1)
text(2.25,0.0275,'n.s.','FontWeight','normal','FontSize',7,'FontName','Arial')
plot([1 3],[0.028 0.028],'k','LineWidth',1)
text(2,0.0285,'***','FontWeight','normal','FontSize',14,'FontName','Arial')
title(figtitle,'FontWeight','normal','FontSize',7,'FontName','Arial');
set(gcf, 'renderer', 'painters')
saveas(h,[PPath.FigPath figtitle]);
saveas(h,[PPath.FigPath figtitle],'svg');

save([PPath.ResPath 'TagCoh_Parallel.mat'],'TagCoh','-v7.3');





 