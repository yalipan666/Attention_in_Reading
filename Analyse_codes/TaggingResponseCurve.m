%%% try to plot the responsive curve of RIFT for each subject at each word
%%% locations
%%% coherence was averaged over 200ms after the fixation onset of each word
%%% or 200ms after the baseline onset as a 'chance level'
% Note: since the skip rate is different for words in different location,
% trial number is different across words. However, the averaging time
% window is the same for each words
% note: analysis timewindow are first-pass durations

function TaggingResponseCurve(sid)
%%% set paths
server = 1; % run analysis on the server or local computer
rootdir = Get_Paths(server); %get the root dir and the colormap of cbrewer(RdBu)
PPath.ResPath = [rootdir 'Results' filesep 'CohCurve' filesep]; %result path
if (~exist(PPath.ResPath,'dir'))
    mkdir(PPath.ResPath);
end
%%% setting
freq_halfwidth = 5; % for hilbert filter
stat_tagfreq = -3:3;% the freq range used to average coherence(centre freq indexes as 0)
%%% get the tag freq for pretarget and target
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
if ExpInfo.FreqVersion(sid) == 1
    tagfreq_all = [60 65];
else
    tagfreq_all = [65 60];
end
pds = {'MISC004','MISC005'}; % photodiodes for the first and second tagging freq (i.e., flicker target and flicker postartget)
EpochType = {'PreTarg';'Targ'};

%%% get the exp information
load([rootdir 'Results' filesep 'Coh' filesep 'TagCoh.mat']);

%%% set the outputs
TagCurve.row_hdr = {'1stTagFreq_PretargetSigSen','2ndTagFreq_TargetSigSen'};
TagCurve.column_hdr = {'baseline','n-4','n-3','n-2','n-1','n','n+1','n+2','n+3','n+4'};
Loc2Target = [nan -4 -3 -2 -1 0 1 2 3 4];
n = length(TagCurve.column_hdr);
TagCurve.coh = zeros(2,n); %[tagfreq1; tagfreq2] 1st row is the 1st tagging freq
TagCurve.trlnum = zeros(2,n);
TagCurve.fixdur_avg = zeros(2,n);
TagCurve.fixdur_std = zeros(2,n);

%%%%%%%%%%%%%%%%%========= loop over epochs ========%%%%%%%%%%%%%%%
for tfi = 1:length(tagfreq_all)
    tagsub = TagCoh.SigSubID{1,tfi};
    sub_id = tagsub(sid); %the real subject id in the whole datasets
    subname = ExpInfo.subjects{sub_id};
    datapath = [rootdir 'Analyse_data' filesep subname filesep];
    eval(['sigchan = TagCoh.' EpochType{tfi} '_SigTagSens{sub_id,2};']);
    %%% get coloum index
    tag_col = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
    fixdur_col = find(strcmp(ExpInfo.EventHdr,'fixation_duration'));
    firstpass_col = find(strcmp(ExpInfo.EventHdr,'FirstPassFix'));
    trig_col = find(strcmp(ExpInfo.EventHdr,'fixation_on_MEG'));
    
    %%% get the freq range
    tagrange = tagfreq_all(tfi)+stat_tagfreq;
    
    %%% run the coherence analysis, loop over words
    for i = 1:n
        clear epoch
        loc = Loc2Target(i);
        if isnan(loc)
            load([datapath 'epoch_BL_Cross.mat']);
            epoch = epoch_BL_Cross;
            clear epoch_BL_Cross
        else
            % loadin data
            if i == 2
                load([datapath 'data_icaclean.mat']); %data
                load([datapath 'Event.mat']);
                load([datapath 'hdr.mat']);
            end
            % get epochs for the word [-0.5 0.5]s
            AnaTW = 1000; %ms
            PPara.pretrig = -AnaTW/2;
            PPara.posttrig = AnaTW/2;
            tw = [80 AnaTW];
            PPara.SR = 1000; %sampling rate
            PPara.timezero = 'fixation_on_MEG';
            PPara.badsens = ExpInfo.BadSensor{sub_id};
            validduration = Event.event_raw(:,fixdur_col)>= tw(1) & Event.event_raw(:,fixdur_col)<= tw(2);
            thisword = Event.event_raw(:,tag_col) == loc;
            firstpass = Event.event_raw(:,firstpass_col) == 1;
            PPara.event_all = Event.event_raw(validduration & thisword & firstpass,:);
            % no need to reject trials
            PPara.rejectvisual = 0;
            epoch = Get_Epoch(hdr,data,PPara,trig_col);
        end
        
        %%%  normalize pd and remove trials that without photodies signal
        pdi = find(strcmp(epoch.label,pds{tfi})); % get pd label index
        rmtrl = [];
        for ttt = 1:length(epoch.trial)
            %%% remove trials that have no pd-004 signal
            if max(epoch.trial{ttt}(pdi,:)) < 0.005
                rmtrl = [rmtrl; ttt];
            end
            epoch.trial{ttt}((pdi),:) = zscore(epoch.trial{ttt}((pdi),:),0,2); % zscore across each trial
        end
        cfg = [];
        cfg.trials = setdiff(1:length(epoch.trial),rmtrl);
        epoch = ft_selectdata(cfg, epoch);
        
        %%% run cohernece
        %%%% just to get a dummy coherence struct
        cfg            = [];
        cfg.output     = 'fourier';
        cfg.channel    = {'MEGGRAD',pds{tfi}}; % the last chan is always the photodiode
        cfg.method     = 'mtmconvol';
        cfg.taper      = 'hanning';
        cfg.foi        = tagrange;
        cfg.toi        = -0.5:0.5:0.5;
        cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
        cfg.keeptrials = 'yes';
        cfg.pad        = 'nextpow2';
        fourier = ft_freqanalysis(cfg,epoch);
        %%% get coherence spctrm
        cfg            = [];
        cfg.method     = 'coh';
        cfg.channelcmb = {'MEGGRAD',pds{tfi}}; %% all channell combinations calculated together
        coh            = ft_connectivityanalysis(cfg, fourier); %% the raw coherence
        coh.label      = fourier.label(1:length(coh.labelcmb));
        coh.time       = epoch.time{1};
        
        %%% empty the cohspctrm; do the real coherence using hilbert complex
        coh.cohspctrm = zeros(length(coh.label),length(coh.freq),length(coh.time));
        for fff = 1:length(coh.freq)
            %%%% get hilbert filert data
            cfg = [];
            cfg.channel    = {'MEGGRAD',pds{tfi}};
            cfg.bpfilter   = 'yes';
            cfg.bpfreq     = [tagrange(fff)-freq_halfwidth tagrange(fff)+freq_halfwidth];
            cfg.hilbert    = 'complex';
            cfg.keeptrials = 'yes';
            fltdata = ft_preprocessing(cfg,epoch);
            for chan = 1:length(fltdata.label)-1
                for ttt = 1:length(fltdata.trial)
                    sig1(:,ttt) = fltdata.trial{ttt}(chan,:); %time*trl
                    sig2(:,ttt) = fltdata.trial{ttt}(end,:);
                end
                spec1 = nanmean(sig1.*conj(sig1),2);%time*1
                spec2 = nanmean(sig2.*conj(sig2),2);
                specX = abs(nanmean(sig1.*conj(sig2),2)).^2;
                coh.cohspctrm(chan,fff,:) = specX./(spec1.*spec2);%time*1
            end
        end
        %%% finally get the sig occipital tagging sensors
        mini_fixdur = 0.2;
        tp = [dsearchn(epoch.time{1}',0):dsearchn(epoch.time{1}',mini_fixdur)];
        cohcoh = squeeze(coh.cohspctrm);%% chan*freq*time
        cohcoh = nanmean(nanmean(nanmean(coh.cohspctrm(sigchan,:,tp),1),2),3); %averaging over chan,freq,and time
        TagCurve.coh(tfi,i) = cohcoh;
        TagCurve.trlnum(tfi,i) = size(epoch.trialinfo,1);
        TagCurve.fixdur_avg(tfi,i) = nanmean(epoch.trialinfo(:,fixdur_col));
        TagCurve.fixdur_std(tfi,i) = nanstd(epoch.trialinfo(:,fixdur_col));
    end
end
save([PPath.ResPath 'TagCurve_' num2str(sid)],'TagCurve')


%% %%%%%%%%%%%%%%%% combine and plot the curve locally %%%%%%%%%%%%%%%%%%
server = 0; % run analysis on the server or local computer
rootdir = Get_Paths(server); %get the root dir and the colormap of cbrewer(RdBu)
PPath.ResPath = [rootdir 'Results' filesep 'CohCurve' filesep]; %result path
TagCurve_all.row_hdr = {'1stTagFreq_PretargetSigSen','2ndTagFreq_TargetSigSen'};
TagCurve_all.column_hdr = {'baseline','n-4','n-3','n-2','n-1','n','n+1','n+2','n+3','n+4'};
n = length(TagCurve_all.column_hdr);
nsub = 29;
TagCurve_all.coh = zeros(2,n,nsub);
TagCurve_all.trlnum = zeros(2,n,nsub);
TagCurve_all.fixdur_avg = zeros(2,n,nsub);
TagCurve_all.fixdur_std = zeros(2,n,nsub);
for sid = 1:nsub
    load([PPath.ResPath 'TagCurve_' num2str(sid)])
    TagCurve_all.coh(:,:,sid) = TagCurve.coh;
    TagCurve_all.trlnum(:,:,sid) = TagCurve.trlnum;
    TagCurve_all.fixdur_avg(:,:,sid) = TagCurve.fixdur_avg;
    TagCurve_all.fixdur_std(:,:,sid) = TagCurve.fixdur_std;
end
clear TagCurve
TagCurve = TagCurve_all;
clear TagCurve_all
save([PPath.ResPath 'TagCurve'],'TagCurve')

%%% plot, bar with error
fignam = 'Tagging response curve';
h = figure('Name',fignam,'color',[1 1 1],'Position',[1 1 800 360]);
for tfi = 1:length(tagfreq_all)
    err = std(TagCurve.coh(tfi,:,:),0,3)/sqrt(nsub);
    meancoh = mean(TagCurve.coh(tfi,:,:),3);
    subplot(1,2,tfi);
    errorbar(meancoh(2:end),err(2:end),'color','k','LineWidth',1)
    xticklabels();
    hold on;
    plot([1 n-1],[meancoh(1) meancoh(1)],'k-.')
    set(gca,'box','off','LineWidth',1)
    set(gca,'XTickLabel',TagCurve.column_hdr(2:end),'FontSize',10,'FontWeight','normal','FontName','Arial')
    title(['tagging location: ' EpochType{tfi+1} 'et']);
    %%% add p value
    for i = 2:n
        [~,p,~,stats] = ttest(TagCurve.coh(tfi,1,:), TagCurve.coh(tfi,i,:));
        text(i-0.6,meancoh(i),num2str(round(1000*p)/1000));
    end
    xlim([0.5 n-1])
    hold on;
end
saveas(h,'TagCurve.fig')




