% measure for reading speed based on the duration of the whole sentence,
% but not the fixation duration of each word (now we also include the
% saccadic durations between words)
% note: since the eyemovement data don't follow the Gaussian distribution,
% it's better to use Spearman correlation rather than Pearson correlation.

%% PPaths
PPath.Tag = 'Z:\Full_Attention\Results\Coh\';
PPath.EM = 'Z:\Full_Attention\Results\Behav\';
PPath.Data = 'Z:\Full_Attention\Analyse_data\';
PPath.PTB = 'Z:\Full_Attention\RawData\PTB_data\';
PPath.Result = 'Z:\Full_Attention\Results\Corr\';
Corr = [];
Corr.CondName = {'Low','High'};

%% ======= get the eye movements metrics for each participant ======= %%
load([PPath.EM 'BehaData']);
load([PPath.Tag 'TagCoh']);
% reading time for the whole sentence for all subjects
nsub = size(TagCoh.subs,1);
DurationPerwrd = nan(nsub,1); %averaged fixation duration of all words
for s = 1:length(TagCoh.subs)
    load([PPath.Data TagCoh.subs{s} '\EyeData.mat']); %EyeData
    load([PPath.Data TagCoh.subs{s} '\Event.mat']); %EyeData
    load([PPath.PTB TagCoh.subs{s} '.mat'],'Para'); %PTB data
    
    % get all fixations of each words
    tmp_allfixdu = cellfun(@(x) x(~isnan(x(:,7)),3),EyeData.TrlFixdata,'Uni',false); %%[allfix_duration]
    %     % remove outliers of fixations
    %     tmp_allfixdu = cellfun(@(x) x(x>=80 & x<= 1000),tmp_allfixdu,'Uni',false);
    % get sum of the fixation durations for each sentences
    tmp_allfixdu = cellfun(@(x) sum(x),tmp_allfixdu,'Uni',true);
    
    % word number in all sentences
    wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
    % divided by the number of words in each sentence
    tmp = tmp_allfixdu./wrdnum;
    DurationPerwrd(s,1) = mean(tmp);
end
BehaData.DurationPerwrd = DurationPerwrd;
save([PPath.EM 'BehaData'],'BehaData');
Corr.DurationPerwrd = BehaData.DurationPerwrd;
Corr.Speed = BehaData.Speed;


%% get the averaged Tagging coherence for each participant
EpochType = {'PreTarg';'Targ'};
figure;
for ep = 1:length(EpochType)
    eval(['Corr.' EpochType{ep} '_SigSub = TagCoh.SigSubID{1,ep};']);
    eval(['Corr.' EpochType{ep} ' = TagCoh.' EpochType{ep} '_stat.data4Ttest;']);
    eval(['cohtmp = Corr.' EpochType{ep} ';']);
    
    %%% =======  Corr of Tagging and EM ======= %%
    % Tag_pretarg diff, low-high
    tag_dif = cohtmp(:,1)-cohtmp(:,2);
    % Tag_pretarg_diff & reading speed
    rs = Corr.Speed(TagCoh.SigSubID{1,ep});
    [coef,pval] = corr([tag_dif rs],'type','Pearson');
    eval(['Corr.' EpochType{ep} '_ReadSped_CoefP = [coef(1,2),pval(1,2)];']);
    [coef(1,2),pval(1,2)] %% Pretarget: [0.1887    0.3180]; Target:[-0.4171    0.0157]
    
    %%% plot the correlation
    subplot(2,2,(ep-1)*2+1)
    scatter(tag_dif,rs,50,'filled')
    title(['Correlation for ' EpochType{ep} 'et'],'FontWeight','bold','FontSize',14);
    ylabel('Reading speed (#words/s)','FontWeight','bold','FontSize',12)
    xlabel([EpochType{ep} 'et coherence difference (r^2)'],'FontWeight','bold','FontSize',12)
    text(0,4.5,['[coef,pval] = [' num2str([coef(1,2),pval(1,2)]) ']'],'FontWeight','bold','FontSize',12);
    
    % averaged first fixation duration of target words
    if ep == 1
        ff = BehaData.FirstFix.Tag(TagCoh.SigSubID{1,ep},:);
        ff_avg = mean(ff,2);
        [coef,pval] = corr([tag_dif ff_avg],'type','Pearson');
        eval(['Corr.' EpochType{ep} '_FFavgTarg_CoefP = [coef(1,2),pval(1,2)];']);
    elseif ep == 2
        ff = BehaData.FirstFix.Pos(TagCoh.SigSubID{1,ep},:);
        ff_avg = mean(ff,2);
        [coef,pval] = corr([tag_dif ff_avg],'type','Pearson');
        eval(['Corr.' EpochType{ep} '_FFavgPosTarg_CoefP = [coef(1,2),pval(1,2)];']);
    end
    [coef(1,2),pval(1,2)] %% Pretarget: [-0.3826    0.0369]; Target:[0.3149    0.0743]
    
    %%% plot the correlation
    subplot(2,2,(ep-1)*2+2)
    scatter(tag_dif,ff_avg,50,'filled')
    title(['Correlation for ' EpochType{ep} 'et'],'FontWeight','bold','FontSize',14);
    ylabel('Fir fix dur of parafoveal words (ms)','FontWeight','bold','FontSize',12)
    xlabel([EpochType{ep} 'et coherence difference (r^2)'],'FontWeight','bold','FontSize',12)
    text(0,280,['[coef,pval] = [' num2str([coef(1,2),pval(1,2)]) ']'],'FontWeight','bold','FontSize',12);
end

%% correlation between attention flexibility & reading behaviour
PPath.Result = 'Z:\Full_Attention\Results\Corr\';
load([PPath.Result 'Corr.mat']);
%%%%%%% ============== individual reading speed ============%%%%%%
%%% coh_diff between pre_diff and targ_diff: need subjects shared in both
%%% pre-target sig and target sig coherence
[sig_both,sid_both_inpre,sid_both_intarg] = intersect(Corr.PreTarg_SigSub,Corr.Targ_SigSub); %only 27 subjects
tmp_pre = Corr.PreTarg(sid_both_inpre,1) - Corr.PreTarg(sid_both_inpre,2); %low-high
tmp_targ = Corr.Targ(sid_both_intarg,1) - Corr.Targ(sid_both_intarg,2); %low-high
diff_diff = tmp_pre - tmp_targ;% pretarget - target (i.e., positive valuse minue negative values)
Corr.SubID_inboth = sig_both;
Corr.Coh_diff_diff = diff_diff;
% correlations with individual reading speed
speed_both = Corr.Speed(sig_both,1);
[coef,pval] = corr([diff_diff speed_both],'type','Pearson'); [coef(1,2),pval(1,2)]
Corr.CohDiffDiff_ReadSped_CoedP = [coef(1,2),pval(1,2)]; %[0.2990    0.1298]
%%% do correlation seperately for pre-target and target
diff_pre = Corr.PreTarg(:,1) - Corr.PreTarg(:,2); %low-high
speed_pre = Corr.Speed(Corr.PreTarg_SigSub,1);
[coef,pval] = corr([diff_pre speed_pre],'type','Pearson'); [coef(1,2),pval(1,2)]
Corr.CohDiffPretarg_ReadSped_CoedP = [coef(1,2),pval(1,2)]; % [0.1887    0.3180]

diff_targ = Corr.Targ(:,1) - Corr.Targ(:,2); %low-high
speed_targ = Corr.Speed(Corr.Targ_SigSub,1);
[coef,pval] = corr([diff_targ speed_targ],'type','Pearson'); [coef(1,2),pval(1,2)]
Corr.CohDiffTarg_ReadSped_CoedP = [coef(1,2),pval(1,2)]; %[-0.4171    0.0157]

%%%%%%% ============== all first fixation durations ============%%%%%%
% correlation with the first fixation duration of all words in the
% sentence set
load('Z:\Full_Attention\Analyse_data\ExpInfo.mat');
nsub = length(ExpInfo.subjects);
FF = zeros(nsub,1);
for i = 1:nsub
    load(['Z:\Full_Attention\Analyse_data\' ExpInfo.subjects{i} '\Event.mat'])
    ff = Event.event_raw(Event.event_raw(:,8)==1,6);
    FF(i) = mean(ff);
end
[coef,pval] = corr(diff_diff, FF(Corr.SubID_inboth),'type','Pearson') %[-0.4678    0.0139]
Corr.CohDiffDiff_FirstFix_CoedP = [coef(1,2),pval(1,2)];
% save out
Corr.FirstFixation_allwords = FF;














