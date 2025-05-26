% copy from Lexical/Analyse_codes
% 20210720 pre-processing data on bluebear: get raw megdata, eyedata, event,
% and all ica components

function S2_PreProcessing(sid) %id of sub
% number of .fif files in SV task
%%% set paths
server = 1;
rootdir = Get_Paths(server);
PPath.RawMEG = [rootdir 'RawData' filesep 'MEG_data' filesep];
PPath.RawPTB = [rootdir 'RawData' filesep 'PTB_data' filesep];
PPath.RawEye = [rootdir 'RawData' filesep 'EyeLink_data' filesep];

%%% settings for pre-processing

MPara.detrend    = 'yes';
MPara.SR         = 1000;
MPara.pretrig    = 0; % for randomly epoching data during ica
MPara.bpfilter   = 'yes';
MPara.bpfreq     = [0.5 100];
% nouth-filter the line noise around 50Hz
MPara.bsfilter = 'yes';
MPara.bsfreq   = [49 51];
%%% get file names
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
subjects = ExpInfo.subjects;
PTBFiles = ExpInfo.PTBFiles;
EyeFiles = ExpInfo.EyeFiles;
sub = subjects{sid};
TrackedEye = ExpInfo.TrackEye{sid};% get tracked eye of eye-link 
MPara.sub = sub;
MPara.badsens = ExpInfo.BadSensor{sid};
File.PTB = PTBFiles{sid};
File.Eye = EyeFiles{sid};
PPath.SaveData = [rootdir 'Analyse_data' filesep sub filesep];
if ~exist(PPath.SaveData,'dir')
    mkdir(PPath.SaveData);
end

%%% display processing sub
disp(['*** PreProcessing: sub--' sub]);

%% get raw MEG data
if ~exist([PPath.SaveData 'data.mat'],'file')
    tmpf = [sub filesep sub(3:8) filesep sub(10:end)];
    File.MEG{1,1} = [tmpf '.fif'];
    % loop to find the append data
    f = 1; % index of append data
    while 1
        tmpmeg = [tmpf '-' num2str(f) '.fif'];
        if exist([PPath.RawMEG tmpmeg],'file')
            File.MEG{1,f+1} = tmpmeg;
        else
            break
        end
        f = f + 1;
    end
    [hdr,data,Trigger_MEG] = Get_MEGData(PPath,File);
    save([PPath.SaveData 'hdr'],'hdr');
    save([PPath.SaveData 'Trigger_MEG'],'Trigger_MEG','-v7.3');
    save([PPath.SaveData 'data'],'data','-v7.3');
else
    load([PPath.SaveData 'hdr']);
    load([PPath.SaveData 'data']);
    load([PPath.SaveData 'Trigger_MEG']);
end

%% remove artefacts with ICA on server
%%% get the data based on the first and alst trigger (get rid of the
%%% cHPI at the begining and end of data collection) with 1 sec extra
MPara.start = Trigger_MEG(1,2)-1000;
MPara.end = Trigger_MEG(end,2)+1000;

%%% band-pass filter data
% get standard fieldtrip strcture
epoch          = [];
epoch.label    = hdr.label;
epoch.trial{1} = data(:,MPara.start:MPara.end);
epoch.time{1}  = (MPara.start:MPara.end)/MPara.SR;
% get bad channels
keepchan = {'MEG'};
for cc = 1:length(MPara.badsens)
    keepchan{cc+1} = ['-MEG' MPara.badsens{cc}];
end
% preprocessing
cfg          = [];
cfg.channel  = keepchan;  %% only valid MEG channel
cfg.detrend  = MPara.detrend;     % Remove slow drifts
cfg.bpfilter = MPara.bpfilter;
cfg.bpfreq   = MPara.bpfreq;
cfg.bsfilter = MPara.bsfilter;
cfg.bsfreq   = MPara.bsfreq;
epoch4ICA    = ft_preprocessing(cfg, epoch);
% save the start and end timepoints
epoch4ICA.trlinfo = [MPara.start MPara.end];% [start end]timepoints
clear epoch
% run ICA
comp = ICA4rawdata(epoch4ICA);
save([PPath.SaveData 'ica.mat'],'epoch4ICA','comp','-v7.3')

%% get eyemovement metrics from eyelink
%first converting data from .edf to .asc (C:/toolbox/SR Research/edfconverter/)
if ~exist([PPath.SaveData 'Event.mat'],'file')
    if ~exist([PPath.SaveData 'EyeData.mat'],'file')
        eyefile = [PPath.RawEye File.Eye '.asc'];
        load([PPath.RawPTB File.PTB],'Result');
        WordLocMat = 2.*Result.WordLocation;
        EyeData = Get_EyeData(eyefile, WordLocMat,ExpInfo.Trigger,TrackedEye);
        save([PPath.SaveData 'EyeData'],'EyeData');
    else
        load([PPath.SaveData 'EyeData']);
    end
    %%% get event
    load([PPath.RawPTB File.PTB],'Para'); % Para
    Event = Get_Event(EyeData,Para.CondMat,Para.FlkWords,Trigger_MEG,ExpInfo.Trigger,ExpInfo.EventHdr);
    save([PPath.SaveData 'Event'],'Event');
end

% end
disp(['*** PreProcessing done! --- ' sub]);
end




