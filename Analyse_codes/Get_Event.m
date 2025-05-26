% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise

function Event = Get_Event(EyeData,CondMat,TargLoc,Trigger_MEG,Trigger,EventHdr)
%%% get triggers from both meg and eyelink
Event.Trigger_MEG = Trigger_MEG;
Event.SentOnTrig_header = {'TriggerValue','EyeLink','MEG'};
eyeidx = EyeData.TriggerMat(:,1) == Trigger.SentOn;
trig_eye = EyeData.TriggerMat(eyeidx,:);
megidx = Trigger_MEG(:,1) == Trigger.SentOn;
trig_meg = Trigger_MEG(megidx,2);
sentid = EyeData.TrlId; % sentence id in the data
n_meg = size(trig_meg,1);
n_eye = size(trig_eye,1);
if n_eye>n_meg %late start of meg recording
    warning('trial number in eye data does not match the brain data, use brain data as a reference to cut eye data')
    id_new = n_eye-n_meg+1:n_eye;
    trig_eye = trig_eye(id_new,:);
    sentid = sentid(id_new);
end
Event.SentOnTrig = [trig_eye trig_meg]; %[trigger eye meg]
% correct eyelink triggers aligned with sentence onset of MEG triggers
eyeneedminus = trig_eye(:,2)-trig_meg; % EyeLink-MEG

%% % get the whole event of all fixations --- raw event
% NextOrder:location of next fixated word compared to current word;
% FirstPassFix: fixate in the 1st pass time
event_raw_header = EventHdr;
event_raw = []; % total fixaiton duration event
TargLoc = TargLoc(sentid,:); %[tarloc1 tarloc2]
CondMat = CondMat(sentid,:);

% get saccade metric
saccade_eye = EyeData.SaccadeData(:,[2 3]); %[saccadeoff_eye saccade_duration]

for t = 1:length(sentid)
    trldata_tmp = EyeData.TrlFixdata{t}; %['StartTime', 'EndTime', 'Duration', 'AveragedXPosition', 'AveragedYPosition', 'AveragedPupilSize','FixWrdId']
    % select valid trldata---fixation on a real word
    valword = ~isnan(trldata_tmp(:,7));
    % get the information of the trial event
    trldata = trldata_tmp(valword,[1 3 7]); %['StartTime','Duration','FixWrdId']
    pupilsize = trldata_tmp(valword,6);
    evtn = size(trldata,1);
    cur_sent = repmat(sentid(t),evtn,1);%mat of current sentence id
    word_loc = trldata(:,3);%word location mat of all fixations
    % find the distance between a given word and the nearest target
    loc2targ = [];
    for www = 1:size(word_loc,1) %loop over all fixations
        cur_wl = word_loc(www); % word location of the current fixation
        tmptargloc = []; %try to figure out to which targe word the current word is closest (Degno's sentence set contain two targets in each sentence)
        for p = 1:size(TargLoc,2)
            tmptargloc = [tmptargloc abs(cur_wl-TargLoc(t,p))];
        end
        [~,near_targ] = min(tmptargloc);
        loc2targ(www,1) = cur_wl-TargLoc(t,near_targ);
    end
    fixon_eye = trldata(:,1);
    sac2this_dur = [];
    for pp = 1:evtn
        r = find(saccade_eye(:,1) == fixon_eye(pp)-1);
        if length(r)==1
            sac2this_dur(pp,1) = saccade_eye(r,2);
        else
            sac2this_dur(pp,1) = nan;
        end
    end
    fix_on_MEG = fixon_eye-repmat(eyeneedminus(t),evtn,1);
    fix_dur = trldata(:,2);
    trl_event = [cur_sent word_loc loc2targ sac2this_dur fix_on_MEG fix_dur];
    
    %%% read direction
    word_loc_plus1 = [word_loc(2:end); nan];
    NextOrder = word_loc_plus1 - word_loc;
    word_loc_minus1 = [nan; word_loc(1:(end-1));];
    PreOrder = word_loc_minus1 - word_loc;
    
    %%% first pass fixation
    %%% all first fixation in the first pass, no re-read after reading the biggest (last) word, could miss some words!
    uniqueword = unique(word_loc,'stable'); %%% all first fixation
    [~,biggest_fixword_id] = max(uniqueword);
    Firstidx = [];
    for pp = 1:biggest_fixword_id
        tmp = find(word_loc == uniqueword(pp));
        Firstidx = [Firstidx;tmp(1)];
    end
    FirstFix = zeros(evtn,1);
    FirstFix(Firstidx,1) = 1;
    %get the whole event
    event_raw = [event_raw;[trl_event NextOrder FirstFix PreOrder repmat(CondMat(t,1),evtn,1) pupilsize]];
end
Event.event_raw_header = event_raw_header;
Event.event_raw = event_raw;


