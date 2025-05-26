%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%========== stat on behavioral data ========== %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
server = 0; % run analysis on the server or local computer
%%% set paths
rootdir = Get_Paths(server);
PPath.ResPath = [rootdir 'Results' filesep 'Behav' filesep]; %result path
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
subjs = ExpInfo.subjects;
CondId = ExpInfo.CondID;

%% % gaze duration: all fixation on target before eye move out of taget,
%%% re-read is not included
BehaData = [];
BehaData.CondName = ExpInfo.CondName;
for sss = 1:length(subjs)
    load([rootdir 'Analyse_data' filesep subjs{sss} filesep 'Event.mat']); %Event
    if sss == 1
        % get the column id for some variables
        sentid = find(strcmp(Event.event_raw_header,'sentence_id'));
        wordloc = find(strcmp(Event.event_raw_header,'word_loc'));
        loc2targ = find(strcmp(Event.event_raw_header,'loc2targ'));
        sacdur = find(strcmp(Event.event_raw_header,'saccade2this_duration')); 
        megtrig = find(strcmp(Event.event_raw_header,'fixation_on_MEG')); 
        fixdur = find(strcmp(Event.event_raw_header,'fixation_duration'));
        firstPass = find(strcmp(Event.event_raw_header,'FirstPassFix'));%note!first pass more means the first fixation 
        cond = find(strcmp(Event.event_raw_header,'SentenceCondition'));
        pupil = find(strcmp(Event.event_raw_header,'PupilSize'));
        preorder = find(strcmp(Event.event_raw_header,'PreviousOrder'));
        %%% pre-config the gaze event
        GazeEvent.hdr = Event.event_raw_header([sentid wordloc loc2targ megtrig fixdur cond]);
    end
    
    %%%%%%%%%%%%%========== get gaze_event and save it out ========== %%%%%%%%%%%%%%
    %%% gaze: sum of the fixations before eyes leave a given word, the saccades duration ARE included here   
    %step1: find all the non-first fixations that in the same word and add
    %the saccade duration before it to the fixation
    tmp = Event.event_raw;
    tmp = tmp(tmp(:,fixdur)>=80,:); % get rid of the gaze that shorter than 80ms
    samewrdfix = tmp(:,preorder)==0; %index for the non-first-fixation before eyes leave this word
    tmp(samewrdfix,fixdur) = nansum([tmp(samewrdfix,fixdur) tmp(samewrdfix,sacdur)],2);
    id_firstfirx = find(tmp(:,firstPass)==1);
    n = length(id_firstfirx);
    GazeEvent.event = zeros(n,length(GazeEvent.hdr));
    for i = 1:n
        cur_i = id_firstfirx(i); %row id in the raw tmp array
        tmpgz = tmp(cur_i,fixdur);
        loopi = cur_i + 1;
        while loopi < n && tmp(loopi,preorder) == 0
            tmpgz = [tmpgz tmp(loopi,fixdur)];
            loopi = loopi + 1;
        end
        GazeEvent.event(i,:) = tmp(id_firstfirx(i),[sentid wordloc loc2targ megtrig fixdur cond]);
        GazeEvent.event(i,end-1) = nansum(tmpgz);
    end
    %%% get gaze data
    gz_pre1 = GazeEvent.event(:,3)==-1 & GazeEvent.event(:,6)==CondId(1);
    gz_pre2 = GazeEvent.event(:,3)==-1 & GazeEvent.event(:,6)==CondId(2);
    Gaze.Pre(sss,:) = [mean(GazeEvent.event(gz_pre1,5)) mean(GazeEvent.event(gz_pre2,5))];
    gz_tag1 = GazeEvent.event(:,3)==0 & GazeEvent.event(:,6)==CondId(1);
    gz_tag2 = GazeEvent.event(:,3)==0 & GazeEvent.event(:,6)==CondId(2);
    Gaze.Tag(sss,:) = [mean(GazeEvent.event(gz_tag1,5)) mean(GazeEvent.event(gz_tag2,5))];
    gz_pos1 = GazeEvent.event(:,3)==1 & GazeEvent.event(:,6)==CondId(1);
    gz_pos2 = GazeEvent.event(:,3)==1 & GazeEvent.event(:,6)==CondId(2);
    Gaze.Pos(sss,:) = [mean(GazeEvent.event(gz_pos1,5)) mean(GazeEvent.event(gz_pos2,5))];
    
    %%%%%%%%%%%%%========== get total gaze duration and regression probability========== %%%%%%%%%%%%%%
    %%% total gze: sum of all the fixations to a given word, the saccades duration IS included here
    togz_reg = zeros(n,4); %[loc2targ cond fixdur num_of_regression]
    for i = 1:n
        cur_i = id_firstfirx(i); %row id in the raw tmp array
        si = tmp(cur_i,sentid);
        wi = tmp(cur_i,loc2targ);
        cur_wrd = tmp(:,sentid)==si & tmp(:,loc2targ)==wi;
        togz_reg(i,1:3) = [wi tmp(cur_i,cond) sum(tmp(cur_wrd,fixdur))]; 
        %%% get the probablity of regressed fixations
        % We calculated the probability of making a regression into a word 
        % as the proportion of trials on which there was at least one regression 
        % from a later part of the sentence back to that word. 
        togz_reg(i,4) = any(tmp(cur_wrd,preorder) > 0);
    end
    % get total gaze durations in each condition
    tg_pre1 = togz_reg(:,1)==-1 & togz_reg(:,2)==CondId(1);
    tg_pre2 = togz_reg(:,1)==-1 & togz_reg(:,2)==CondId(2);
    TotalGaze.Pre(sss,:) = [mean(togz_reg(tg_pre1,3)) mean(togz_reg(tg_pre2,3))];
    tg_tag1 = togz_reg(:,1)==0 & togz_reg(:,2)==CondId(1);
    tg_tag2 = togz_reg(:,1)==0 & togz_reg(:,2)==CondId(2);
    TotalGaze.Tag(sss,:) = [mean(togz_reg(tg_tag1,3)) mean(togz_reg(tg_tag2,3))];
    tg_pos1 = togz_reg(:,1)==1 & togz_reg(:,2)==CondId(1);
    tg_pos2 = togz_reg(:,1)==1 & togz_reg(:,2)==CondId(2);
    TotalGaze.Pos(sss,:) = [mean(togz_reg(tg_pos1,3)) mean(togz_reg(tg_pos2,3))];
    % get regression probability in each condition
    RegsProb.Pre(sss,:) = [mean(togz_reg(tg_pre1,4)) mean(togz_reg(tg_pre2,4))].*100;
    RegsProb.Tag(sss,:) = [mean(togz_reg(tg_tag1,4)) mean(togz_reg(tg_tag2,4))].*100;
    RegsProb.Pos(sss,:) = [mean(togz_reg(tg_pos1,4)) mean(togz_reg(tg_pos2,4))].*100;

    %%%%%%%%%%%%%========== get first fixation ========== %%%%%%%%%%%%%%
    %%% remove outliers, only for the first fixation
    tmp = Event.event_raw;
    tmp = tmp(tmp(:,fixdur)>=80 & tmp(:,fixdur)<=1000,:);
    %%% get duration & pupil size: [cond1 cond2]
    % index for first fixation
    ff = tmp(:,firstPass)==1;
    pre1 = tmp(:,loc2targ)==-1 & tmp(:,cond)==CondId(1);
    pre1_ff = tmp(pre1 & ff,[fixdur pupil]);
    pre2 = tmp(:,loc2targ)==-1 & tmp(:,cond)==CondId(2);
    pre2_ff = tmp(pre2 & ff,[fixdur pupil]);
    tag1 = tmp(:,loc2targ)==0 & tmp(:,cond)==CondId(1);
    tag1_ff = tmp(tag1 & ff,[fixdur pupil]);
    tag2 = tmp(:,loc2targ)==0 & tmp(:,cond)==CondId(2);
    tag2_ff = tmp(tag2 & ff,[fixdur pupil]);
    pos1 = tmp(:,loc2targ)==1 & tmp(:,cond)==CondId(1);
    pos1_ff = tmp(pos1 & ff,[fixdur pupil]);
    pos2 = tmp(:,loc2targ)==1 & tmp(:,cond)==CondId(2);
    pos2_ff = tmp(pos2 & ff,[fixdur pupil]);
    %%% duration of the first fixation 
    FirstFix.Pre(sss,:) = [mean(pre1_ff(:,1))   mean(pre2_ff(:,1))];
    FirstFix.Tag(sss,:) = [mean(tag1_ff(:,1))   mean(tag2_ff(:,1))];
    FirstFix.Pos(sss,:) = [mean(pos1_ff(:,1))   mean(pos2_ff(:,1))];
    %%% pupil size of the first fixation
    PupilSize.Pre(sss,:) = [mean(pre1_ff(:,2))   mean(pre2_ff(:,2))];
    PupilSize.Tag(sss,:) = [mean(tag1_ff(:,2))   mean(tag2_ff(:,2))];
    PupilSize.Pos(sss,:) = [mean(pos1_ff(:,2))   mean(pos2_ff(:,2))];
end

%% simple paired t-tests
% durations of first fixations
[~,p,~,stat] = ttest(FirstFix.Pre(:,1),FirstFix.Pre(:,2));
stat.p = p
FirstFix.Pre_stat = stat;
[~,p,~,stat] = ttest(FirstFix.Tag(:,1),FirstFix.Tag(:,2));
stat.p = p
FirstFix.Tag_stat = stat;
[~,p,~,stat] = ttest(FirstFix.Pos(:,1),FirstFix.Pos(:,2));
stat.p = p
FirstFix.Pos_stat = stat;
% durations of gze fixations
[~,p,~,stat] = ttest(Gaze.Pre(:,1),Gaze.Pre(:,2));
stat.p = p
Gaze.Pre_stat = stat;
[~,p,~,stat] = ttest(Gaze.Tag(:,1),Gaze.Tag(:,2));
stat.p = p
Gaze.Tag_stat = stat;
[~,p,~,stat] = ttest(Gaze.Pos(:,1),Gaze.Pos(:,2));
stat.p = p
Gaze.Pos_stat = stat;
% durations of total gze 
[~,p,~,stat] = ttest(TotalGaze.Pre(:,1),TotalGaze.Pre(:,2));
stat.p = p
TotalGaze.Pre_stat = stat;
[~,p,~,stat] = ttest(TotalGaze.Tag(:,1),TotalGaze.Tag(:,2));
stat.p = p
TotalGaze.Tag_stat = stat;
[~,p,~,stat] = ttest(TotalGaze.Pos(:,1),TotalGaze.Pos(:,2));
stat.p = p
TotalGaze.Pos_stat = stat;
% pupil size of first fixation
[~,p,~,stat] = ttest(PupilSize.Pre(:,1),PupilSize.Pre(:,2));
stat.p = p
PupilSize.Pre_stat = stat;
[~,p,~,stat] = ttest(PupilSize.Tag(:,1),PupilSize.Tag(:,2));
stat.p = p
PupilSize.Tag_stat = stat;
[~,p,~,stat] = ttest(PupilSize.Pos(:,1),PupilSize.Pos(:,2));
stat.p = p
PupilSize.Pos_stat = stat;

% probability of regression  
[~,p,~,stat] = ttest(RegsProb.Pre(:,1),RegsProb.Pre(:,2));
stat.p = p
RegsProb.Pre_stat = stat;
[~,p,~,stat] = ttest(RegsProb.Tag(:,1),RegsProb.Tag(:,2));
stat.p = p
RegsProb.Tag_stat = stat;
[~,p,~,stat] = ttest(RegsProb.Pos(:,1),RegsProb.Pos(:,2));
stat.p = p
RegsProb.Pos_stat = stat;



%% get the individual reading speed
%%%  estimating the sentence duration by getting the first and final fixations of this sentence
PPath.Data = 'Z:\Full_Attention\Analyse_data\';
PPath.PTB = 'Z:\Full_Attention\RawData\PTB_data\';
nsub = length(subjs);
Speed = zeros(nsub,1);
for sss = 1:nsub
    % number of words in each sentence
    load([PPath.PTB subjs{sss} '.mat'],'Para'); %PTB data
    wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
    % duration of sentence presentation/reading
    load([rootdir 'Analyse_data' filesep subjs{sss} filesep 'EyeData.mat']);
    %%% unify the number of sentences in both eye and PTB data
    wrdnum = wrdnum(EyeData.TrlId);
    sent_dur = nan(size(wrdnum));
    for i = 1:length(sent_dur)
        tmp = EyeData.TrlFixdata{i,1};
        valid_id = find(~isnan(tmp(:,7)));
        if isempty(valid_id) % no eye movement data in this sentence
             sent_dur(i,1) = nan;
        else
            sent_dur(i,1) = (tmp(valid_id(end),2) - tmp(valid_id(1),1))./1000; %fixation_off - fixation_on; unit in second
        end
    end
    % average across sentences to get individual reading speed
    Speed(sss,1) = nanmean(wrdnum./sent_dur);
end

%% save out
BehaData.FirstFix = FirstFix;
BehaData.Gaze = Gaze;
BehaData.TotalGaze = TotalGaze;
BehaData.PupilSize = PupilSize;
BehaData.RegsProb = RegsProb;
BehaData.Speed = Speed;
save([PPath.ResPath 'BehaData.mat'],'BehaData');

%% %%%%%%%%%%%%%=============== plotting =================%%%%%%%%%%%%%%%%
FigNam = {'FirstFix';'Gaze';'TotalGaze';'RegsProb'};
colmat = [188 188 0;0 188 188]./255;
EpochType = {'Pre','Tag','Pos'};
grouporder={'Low','High'};
for ff = 1:length(FigNam)
    eval(['tmp = BehaData.' FigNam{ff} ';']);
    if ff == 4
        unit = ' (%)';
    else
        unit = ' (ms)';
    end
    figtitle = [FigNam{ff} '_Ttest'];
    nsigsub = size(tmp.Pre,1);
    group = [cellstr(repmat('Low',nsigsub,1)); cellstr(repmat('High',nsigsub,1))];
    figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 700 280]);
    subtitles = {'Pre-target','Target','Post-target'};
    for mmm = 1:length(EpochType)
        eval(['vdata = [tmp.' EpochType{mmm} '(:,1); tmp.' EpochType{mmm} '(:,2)];']);
        h = subplot(1,length(EpochType),mmm);
        vp = violinplot(vdata, group,'GroupOrder',grouporder);
        vp(1).ViolinColor = colmat(1,:);
        vp(2).ViolinColor = colmat(2,:);
        vp(1).ShowMean = 1; vp(2).ShowMean = 1;
        vp(1,1).MeanPlot.LineWidth = 3;
        vp(1,2).MeanPlot.LineWidth = 3;
        vp(1,1).MedianPlot.Visible = 'off';
        vp(1,2).MedianPlot.Visible = 'off';
        vp(1,2).BoxWidth = 0.01; vp(1,2).BoxWidth = 0.01;
        ylabel([FigNam{ff} unit],'FontSize',10,'FontWeight','normal','FontName','Arial');
        xlabel('Target lexical frequency','FontSize',10,'FontWeight','normal');
        set(gca,'FontSize',10,'FontWeight','normal','FontName','Arial');
        set(gca,'box','off','LineWidth',1)
        %%% plot the line linking each subject
        hold on;
        x1 = vp(1,1).ScatterPlot.XData;
        y1 = vp(1,1).ScatterPlot.YData;
        x2 = vp(1,2).ScatterPlot.XData;
        y2 = vp(1,2).ScatterPlot.YData;
        plot([x1; x2],[y1; y2],'Color',[.9 .9 .9],'linewidth',0.5)
        %%% add subtitles and stats
        eval(['p = tmp.' EpochType{mmm} '_stat.p;']);
        p = round(p,3);
        title([subtitles{mmm} ',p=' num2str(p)],'FontSize',10,'FontWeight','normal','FontName','Arial');
        xlim([0.5 2.5])
    end
    saveas(h,[PPath.ResPath figtitle]);
    %%% save as svg
    set(gcf, 'renderer', 'painters')
    saveas(gcf,[PPath.ResPath figtitle],'svg');
end


%% =========== acc for yes-or-no question ===========%%%%%%%%
PTBpath = [rootdir 'RawData' filesep 'PTB_data' filesep];
ptb_file = ExpInfo.PTBFiles;
acc = [];
for sss = 1:length(ptb_file)
   load([PTBpath ptb_file{sss}]);
   tmp = nanmean(Result.CORR);
   acc = [acc; tmp];
end
meanacc = mean(acc);
stdacc = std(acc);
BehaData.YesNoQ_AvgStd = [meanacc stdacc]; %%[0.941,0.054]
save([PPath.ResPath 'BehaData.mat'],'BehaData');


