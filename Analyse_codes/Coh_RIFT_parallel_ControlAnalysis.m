%%% 2025-05-09  add control analysis to address the concern that ""Perhaps
%%% this simultaneous coherence for N and N+1 is simply because on some
%%% trials attention was on N and on other trials it was on N+1; the resulting
%%% figure shows simultaneous coherence for both?"

function Coh_RIFT_parallel_ControlAnalysis
%%% set paths
server = 0; % run analysis on the server or local computer
[rootdir] = Get_Paths(server); %get the root dir and the colormap of cbrewer(RdBu)
PPath.ResPath = [rootdir 'Results' filesep 'Coh_parallel' filesep]; %result path
if (~exist(PPath.ResPath,'dir'))
    mkdir(PPath.ResPath);
end
%%% setting
freq_halfwidth = 5; % for hilbert filter
SenSelectP = 0.01; % sig threshold for selection tagging sensors
EpochType = {'Targ'}; % epochs of interest
loc2tar = [0];

%%% load in the tagsigsensors
load([PPath.ResPath 'Sensors4ControlAnalysis.mat']); %[foveal parafoveal]

%%% get the exp information
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
tag_col = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
subjects = ExpInfo.subjects;

%%% set output
Results = [];
Results.hdr = {'{sub}','trls*[fovea parafoveal noise/baseline]'};

%% %%%% ============= subjects run! ===============%%%%%%%%%
s = 0; %index for the subjects in analysis
for sid = 1:length(subjects)
    
    %only run for participants with tagsigsensors for both foveal and parafoveal tagging
    if  any(Sensors4ControlAnalysis{sid,1}==0) || any(Sensors4ControlAnalysis{sid,2}==0)
        continue
    end
    s = s+1;
    
    fprintf(['***** analyzing: s' num2str(sid) '**** \n\n']);
    sub = subjects{sid};
    
    %%% get the tag freq for pretarget and target
    if ExpInfo.FreqVersion(sid) == 1
        tagfreq = [60 65]; %%[foveal parafoveal] tagging responses
    else
        tagfreq = [65 60];
    end
    
    %%%%%%%% 1. loading data
    PPath.SaveData = [rootdir 'Analyse_data' filesep sub filesep];
    load([PPath.SaveData 'epoch_WrdOn']); % epoch
    
    % get the time_index of timewin [0 0.2]
    tim_0 = find(epoch.time{1,1} == 0);
    tim_200 = find(epoch.time{1,1} == 0.2);
    tim_id = tim_0:tim_200;
    Pow = [];
    for ep = 1:2 %loop over foveal and parafoveal
        SigTagSens = Sensors4ControlAnalysis{sid, ep};
        % select the epoch data
        flickid     = epoch.trialinfo(:,tag_col)==loc2tar;
        cfg         = [];
        cfg.trials  = find(flickid);
        cfg.channel = epoch.label(SigTagSens);
        data        = ft_selectdata(cfg, epoch);
        
        %%%%%%%%%%% get the power %%%%%%%%%%%
        % get hilbert filert data
        cfg = [];
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = [tagfreq(ep)-freq_halfwidth tagfreq(ep)+freq_halfwidth];
        cfg.keeptrials = 'yes';
        cfg.hilbert    = 'abs';
        amplitude = ft_preprocessing(cfg,data);  % chan*time*trls
        
        % get the averaged power over sigsensors and across the timewin [0 0.2]
        amp_avg = cellfun(@(x) nanmean(nanmean(x(:,tim_id))),amplitude.trial,'Uni',false);
        Pow(:,ep) = cell2mat(amp_avg).^2; %trls*3 --[fovea parafovea baselline/noise]
    end
    % get the baseline/noise level power at 55Hz
    freq_noise = 55;
    cfg = [];
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [freq_noise-1 freq_noise+1]; %[54 56]
    cfg.keeptrials = 'yes';
    cfg.hilbert    = 'abs';
    amplitude = ft_preprocessing(cfg,data);  % chan*time*trls
    % get the averaged power over sigsensors and across the timewin [0 0.2]
    amp_avg = cellfun(@(x) nanmean(nanmean(x(:,tim_id))),amplitude.trial,'Uni',false);
    Pow(:,3) = cell2mat(amp_avg).^2; %trls*3 --[fovea parafovea baselline/noise]
    
    %%% save it out
    Results.Pow{s,1} = Pow;
end

%%% get the statistics
% get the divergence index for power at foveal and parafoveal:
% [fovea-para]./[fovea+para]--> near zero means similar values between
% fovea and para
clear epoch data


%%% plot all the individual power-differences
colmat = [128 128 128;0 114 189;217 83 25]./255; %[grey:foveal-para, blue:Fovea-noise, orange:Para-noise]
figtitle = 'Single-trl PowDiff_subject plot ';
n = length(Results.Pow);
id_1 = [1 1 2]; id_2 = [2 3 3];
grouporder={'Fovea Vs. Para','Fovea Vs. Noise','Para Vs. Noise'};
figure('Name',figtitle,'color',[1 1 1],'Position',[0 0 660 320]);
for p = 1:3 
    subplot(1,3,p);
    c_map = colmat(p,:);
    for s = 1:n
        pow_all = Results.Pow{s,1};
        ntrl = size(pow_all,1);
        x = ones(ntrl,1).*s;
        y = [pow_all(:,id_1(p))-pow_all(:,id_2(p))]./[pow_all(:,id_1(p))+pow_all(:,id_2(p))];
        scatter(x,y,5,c_map,'filled')
        ylim([-1.2 1.2])
        set(gca,'FontSize',10,'FontWeight','normal','FontName','Arial');
        if p == 1
            xlabel('Participant ID','FontSize',10,'FontWeight','normal','FontName','Arial');
            ylabel('Normalized power difference (n.a.)','FontSize',10,'FontWeight','normal','FontName','Arial');
        end
        title(grouporder{p},'FontWeight','normal','FontSize',12,'FontName','Arial');
        plot([1 n],[0 0],'r--')
        hold on;
    end
end
set(gcf, 'renderer', 'painters')
saveas(gcf,figtitle,'fig');
saveas(gcf,figtitle,'svg');


%%% histgram
fovea_para_all = [];
fovea_BL_all = [];
para_BL_all = [];
pow_mini_BL = [];%min(fovea,para)-noise
for s = 1:n
    y = Results.Pow{s,1};
    n_y = size(y,1);
    sub_id = ones(n_y,1).*s;
    fovea_para_all = [fovea_para_all; [sub_id  (y(:,1)-y(:,2))./(y(:,1)+y(:,2))]];
    fovea_BL_all = [fovea_BL_all; [sub_id  (y(:,1)-y(:,3))./(y(:,1)+y(:,3))]];
    para_BL_all = [para_BL_all; [sub_id  (y(:,2)-y(:,3))./(y(:,2)+y(:,3))]];
    
    mini_BL = nan(n_y,1);
    for t = 1:n_y
        mini_BL(t,:) = min(y(t,[1 2]))-y(t,3);
    end
    pow_mini_BL = [pow_mini_BL; [sub_id  mini_BL]];
end
%%% plot the histogram
figtitle = 'Single-trl PowDiff_Hist';
figure('Name',figtitle,'color',[1 1 1],'Position',[0 0 450 200]);
subplot(1,3,1);
histogram(fovea_para_all(:,2),1000,'Normalization','probability','EdgeColor',colmat(1,:,:))
hold on; plot([0 0],[0 0.02],'r--')
set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
xlabel('Normalized power difference (n.a.)','FontSize',7,'FontWeight','normal','FontName','Arial');
ylabel('Histogram probability','FontSize',7,'FontWeight','normal','FontName','Arial');
xlim([-1.1 1.1]);ylim([0 0.016])

subplot(1,3,2);
histogram(fovea_BL_all(:,2),1000,'Normalization','probability','EdgeColor',colmat(2,:,:))
hold on; plot([0 0],[0 0.02],'r--')
set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
xlim([-1.1 1.1]);ylim([0 0.016])

subplot(1,3,3);
histogram(para_BL_all(:,2),1000,'Normalization','probability','EdgeColor',colmat(3,:,:))
hold on; plot([0 0],[0 0.02],'r--')
set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
xlim([-1.1 1.1]);ylim([0 0.016])

%%% save the figure out
set(gcf, 'renderer', 'painters')
saveas(gcf,figtitle,'fig');
saveas(gcf,figtitle,'svg');


%%% save the data out into excel for R (LMM)
filename = 'SingleTrl_PowDiff.xlsx';
range = 'A2';
sheet = "fovea_para";
writematrix(fovea_para_all, filename,'Sheet', sheet, 'Range', range);
sheet = "fovea_noise";
writematrix(fovea_BL_all, filename,'Sheet', sheet, 'Range', range);
sheet = "para_noise";
writematrix(para_BL_all, filename,'Sheet', sheet, 'Range', range);
sheet = "minipow_noise";
writematrix(pow_mini_BL, filename,'Sheet', sheet, 'Range', range);




