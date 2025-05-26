% Get eye movement data from EyeLink and plot the eyemovement trace
% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise

function EyeData = Get_EyeData(eyefile, WordLocMat,Trigger,TrackedEye)
%%% input paras:
% eyefile: .asc file of the eyelink file
% WordLocMat: matrix of xy coordinates of words [trial*words*4]
% Trigger: struct of triggers, SentOn and SentOff are needed here 
% TrackedEye: the eye that was tracked by eye-tracker, 'L' or 'R'

tic
%%% Get important event details from eyelink file 
if isempty(eyefile)
    error('No eyelink files found..');
end
%Open eyelink ASC file
fid = fopen(eyefile);
% reading file line by line -- loop
BlockEnd = 0;
TriggerMat = [];%% getting triggers
TrlFixdata = cell(size(WordLocMat,1),1);
AllEyeEvents = {};
eidx = 0; %% index for all events
SentenceScr = 0; %% whether a sentence is on the screen or not
TrlId = []; %% index of trial
while ~BlockEnd
    strline = fgetl(fid);
    isevent = 0;
    %%% if this line is a sentence start
    if contains(strline, 'Sentence')
        ti = sscanf(strline, 'MSG %f Sentence_ %f');
        ti = ti(2); %% index of the trial
        isevent = 1;
        eidx = eidx + 1;
        txt = 'Sentence';
        vals = ti;
        display(['sentence-' num2str(ti)]);
        TrlId = [TrlId; ti];
    end
    %%% if this line is a fixation start
    if contains(strline, 'SFIX') 
        vals = sscanf(strline, ['SFIX ' TrackedEye ' %f']); %%vals-[StartTime]
        isevent = 1;
        eidx = eidx + 1;
        txt = 'SFIX';
    end
    %%% if this line is a fixation end
    if contains(strline, 'EFIX') && exist('ti','var')
        vals = sscanf(strline, ['EFIX ' TrackedEye ' %f %f %f %f %f %f']); %%vals-[StartTime EndTime Duration AveragedXPosition AveragedYPosition AveragedPupilSize]
        if SentenceScr == 1
            TrlFixdata{ti,1} = [TrlFixdata{ti,1}; vals'];
        end
        isevent = 1;
        eidx = eidx + 1;
        txt = 'EFIX';
    end
     %%% if this line is a saccade start
    if contains(strline, 'SSACC')
        vals = sscanf(strline, ['SSACC ' TrackedEye ' %f']); %%vals-[STime]
        isevent = 1;
        eidx = eidx + 1;
        txt = 'SSACC';
    end
    %%% if this line is a saccade end
    if contains(strline, 'ESACC')
        vals = sscanf(strline, ['ESACC ' TrackedEye ' %f %f %f %f %f %f %f %f %f']); %%vals-[STime ETime Dur StartXPosit StartYPosit EndXPosit EndYPosit Amplitude(degree) PeakVelocity(deg/sec)]
        isevent = 1;
        eidx = eidx + 1;
        txt = 'ESACC';
    end
    %%% if this line is a blink start
    if contains(strline, 'SBLINK')
        vals = sscanf(strline, ['SBLINK ' TrackedEye ' %f']); %%vals-[STime]
        isevent = 1;
        eidx = eidx + 1;
        txt = 'SBLINK';
    end
    %%% if this line is a blink end
    if contains(strline, 'EBLINK')
        vals = sscanf(strline, ['EBLINK ' TrackedEye ' %f %f %f']); %%vals-[STime ETime Dur]
        isevent = 1;
        eidx = eidx + 1;
        txt = 'EBLINK';
    end
    %%% if this line is a trigger
    if contains(strline, 'Trigger')
        vals = sscanf(strline, 'MSG %f Trigger_ %f');
        isevent = 1;
        eidx = eidx + 1;
        txt = 'Trigger';
        TriggerMat = [TriggerMat;[vals(2) vals(1)]];
        if vals(2) == Trigger.SentOn
            SentenceScr = 1;
        end
        if vals(2) == Trigger.SentOff
            SentenceScr = 0;
        end
    end
    if isevent 
        AllEyeEvents{eidx,1} = txt;
        AllEyeEvents{eidx,2} = vals;
    end
    if contains(strline, 'end of block') 
        BlockEnd = 1;
    end
end
fclose(fid);
EyeData.TriggerMat_header = {'Trigger','TimeStamps'};
EyeData.TriggerMat = TriggerMat;
EyeData.AllEyeEvents = AllEyeEvents;

%% get the fixation of each word from eye link
for tt = 1:size(WordLocMat,1) %% trial index
    fixdata = TrlFixdata{tt,1};
    if isempty(fixdata) %%% no dixation data during this trial, accidentally end the sentence before really read it
        fixdata = nan(1,6);
    end
    allwordloc = squeeze(WordLocMat(tt,:,:));
    %%% exclude the non words coordinates
    allwordloc(~allwordloc) = nan; 
    %%% get the space between words
    xspace = allwordloc(2,1) - allwordloc(1,3);
    %%% define the word x-coordinates: the space BEFORE the word and the
    %%% x-length of the word!
    xposlim = [allwordloc(:,1)-xspace allwordloc(:,3)];
    %%% extend the y coordinates by +-200 pixel (or other random size--y coordinate doesn't matter here)
    yposlim = [allwordloc(:,2)-200 allwordloc(:,4)+200];
    %%% get word index for a given eye-fixation
    wrdid = [];
    for ff = 1:size(fixdata,1)
        x = fixdata(ff,4);
        y = fixdata(ff,5);
        %%% use x y coordinates to define the word id of fixation 
        tempid = find(x>xposlim(:,1) & x<xposlim(:,2) & y>yposlim(:,1) & y<yposlim(:,2));
        if isempty(tempid)
            wrdid(ff,1) = nan;
        else
            wrdid(ff,1) = tempid;
        end
    end
    TrlFixdata{tt,1}(:,7) = wrdid;
end
%% remove empty cells from TrlFixdata (no fixations for a given trial)
TrlFixdata = TrlFixdata(TrlId,:); %more than one fixation in this trial
EyeData.TrlFixdata_header = {'StartTime', 'EndTime', 'Duration', 'AveragedXPosition', 'AveragedYPosition', 'AveragedPupilSize','FixWrdId'};
EyeData.TrlFixdata = TrlFixdata;
EyeData.TrlId = TrlId;

%% get the saccade data, pure saccade event WITHOUT EYE BLINK
SSacc_id = find(strcmp(AllEyeEvents(:,1),'SSACC'));
ESacc = strcmp(AllEyeEvents(:,1),'ESACC');
ESacc_suppose = zeros(size(ESacc));
ESacc_suppose(SSacc_id+1) = 1;
if length(ESacc_suppose)>length(ESacc)
    ESacc_suppose = ESacc_suppose(1:end-1);
end
Sacc_enevt_id = logical(ESacc.*ESacc_suppose); %WITHOUT EYE BLINK
tmp = AllEyeEvents(Sacc_enevt_id,2);
datalength = cellfun(@length,tmp);
nonsaccade = datalength~=9;
tmp(nonsaccade) = [];
tmp = cellfun(@(x) x',tmp,'UniformOutput', false);
SaccadeData = cell2mat(tmp);
SaccadeData_header = {'STime','ETime','Dur','StartXPosit','StartYPosit','EndXPosit','EndYPosit','Amplitude(degree)','PeakVelocity(deg/sec)'};
EyeData.SaccadeData_header = SaccadeData_header;
EyeData.SaccadeData = SaccadeData;
 
     
% %% Step-6: plot the eye movements trace
% t = 60; %%% trial id, randomly chosen
% %%% using eye data from EyeLink
% Edata = EyeData.TrlFixdata{t,1};
% WrdLocs = 2.*squeeze(Result.WordLocation(t,:,:));
% p = figure('Name',File.Eye,'color',[1 1 1]);
% subplot(2,1,1)
% title(['Eyemovement Trace of Trl-' num2str(t) '-EyeLink Data']);
% xlim([1 1920])
% ylim([1 1080]) %resolution of eye-link screen
% set(gca,'xlim',[1 1920],'ylim',[1 1080], 'YDir','reverse')
% % plot word locations
% x = [WrdLocs(:,1),WrdLocs(:,3),WrdLocs(:,3),WrdLocs(:,1),WrdLocs(:,1)]';
% y = [WrdLocs(:,2),WrdLocs(:,2),WrdLocs(:,4),WrdLocs(:,4),WrdLocs(:,2)]';
% hold on;
% line(x,y);
% plot eye movements dynamically
% cvalue = colormap(jet(size(Edata,1)));
% for i = 1:size(Edata,1)
%     hold on;
%     scatter(Edata(i,4),Edata(i,5), 20, cvalue(i,:),'filled');
%     text(Edata(i,4),Edata(i,5),num2str(Edata(i,7)));
%     pause(Edata(i,3)/1000)
% end
% 
% %%% using eye data from PTB
% Edata = Result.EYEdata{t,1};
% WrdLocs = 2.*squeeze(Result.WordLocation(t,:,:));
% subplot(2,1,2)
% title(['Eyemovement Trace of Trl-' num2str(t) '-PTB Data']);
% xlim([1 1920])
% ylim([1 1080])
% set(gca,'xlim',[1 1920],'ylim',[1 1080], 'YDir','reverse')
% % plot word locations
% x = [WrdLocs(:,1),WrdLocs(:,3),WrdLocs(:,3),WrdLocs(:,1),WrdLocs(:,1)]';
% y = [WrdLocs(:,2),WrdLocs(:,2),WrdLocs(:,4),WrdLocs(:,4),WrdLocs(:,2)]';
% hold on;
% line(x,y);
% % plot eye movements dynamically
% cvalue = colormap(jet(size(Edata,1)));
% for i = 1:size(Edata,1)
%     hold on;
%     scatter(Edata(i,1),Edata(i,2), 5, cvalue(i,:));
%     pause(0.001);
% end
% saveas(p,[PPath.SaveData 'EyeDataPlot'])

toc
end

