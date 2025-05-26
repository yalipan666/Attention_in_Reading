% % % copy from Lexical/Analyse_codes
% % % 20210719 remove artefcts based on ica components, get all kinds of epochs based on
% % % event

function S3_Get_all_epoches(server,sid)
%%% set paths
rootdir = Get_Paths(server);

%%% basic settingup
RunCond = 'WrdOn'; %%'WrdOn';%%% epoch aligned with RunCond
PPara.filename = RunCond;
AnaTW = 1000;
DoBaseline = 1; % get epochs for baseline period
ica_rmcomps = 1; % remove ica components
%%% get file names
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
epoch_select = 'abs(PPara.event_all(:,3))<2'; %% select: n-1,n,n+1 [pretart,targ,postarg]
conds = ExpInfo.CondID;% condition mat for diff tasks
subjects = ExpInfo.subjects;
sub = subjects{sid};
PPara.sub = sub;
PPara.badsens = ExpInfo.BadSensor{sid};
PPath.SaveData = [rootdir 'Analyse_data' filesep sub filesep];

%% ========== artefacts removal based on ICA components after preprocessing
if ica_rmcomps
    load([PPath.SaveData 'ica']); %'epoch4ICA','comp'
    
    %%% plot the components for visual inspection
    figure
    cfg           = [];
    cfg.layout    = 'neuromag306mag.lay'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    cfg.colormap  = 'jet';
    cfg.component = 1:30;       % specify the component(s) that should be plotted
    ft_topoplotIC(cfg, comp)
   
    %%% identifying components
    cfg            = [];
    cfg.channel    = 1:15;
    cfg.continuous = 'no';
    cfg.viewmode   = 'component';
    cfg.layout     = 'neuromag306mag.lay';
    ft_databrowser(cfg, comp);
    colormap jet;
    
    %%% rejecting components back to the un-down-sampling data
    cfg = [];
    cfg.component = input('component(s) ID that needed to be removed []: ');
    epoch4ICA = ft_rejectcomponent(cfg, comp, epoch4ICA);
    clear comp
    
    %%% put the ica clean epoch data back to the raw data
    load([PPath.SaveData 'data']); %raw meg data
    load([PPath.SaveData 'hdr']); %meg hdr
    senid = cellfun(@(x) find(strcmp({x},hdr.label)),epoch4ICA.label);
    tp = epoch4ICA.trlinfo; % time points for the start and end
    data(senid,tp(1):tp(2)) = epoch4ICA.trial{1};
    clear epoch4ICA
    save([PPath.SaveData 'data_icaclean.mat'], 'data','-v7.3')
    
    % delete ica.mat
    delete([PPath.SaveData 'ica.mat']); %'epoch4ICA','comp'
else
    load([PPath.SaveData 'data_icaclean.mat']); % data
    load([PPath.SaveData 'hdr']); %meg hdr
end

%% ================ epoching ================== %%
load([PPath.SaveData 'Event'])
%%% get trialinfo index
PPara.SR = 1000; %sampling rate
PPara.timezero = 'fixation_on_MEG';
trig_col = find(strcmp(Event.event_raw_header,PPara.timezero));
fixdur_col = find(strcmp(Event.event_raw_header,'fixation_duration'));
firstpass_col = find(strcmp(Event.event_raw_header,'FirstPassFix'));

%%%% get the baseline epochs
if DoBaseline
    TrigOn = Event.Trigger_MEG(Event.Trigger_MEG(:,1) == ExpInfo.Trigger.Fix,2);
    epoch_leng = AnaTW;
    event_all = nan(length(TrigOn),size(Event.event_raw,2));
    event_all(:,trig_col) = TrigOn;
    event_all(:,fixdur_col) = epoch_leng.*ones(size(event_all,1),1);
    PPara.event_all = event_all;
    PPara.pretrig = 0;
    PPara.posttrig = PPara.pretrig + epoch_leng; %epoch end timepoint, aligned with marker
    epoch_BL_Cross = Get_Epoch(hdr,data,PPara,trig_col);
    save([PPath.SaveData 'epoch_BL_Cross'],'epoch_BL_Cross','-v7.3');
end

%%% get the epochs for sentence reading
PPara.pretrig = -AnaTW/2;
PPara.posttrig = AnaTW/2;
tw = [80 1000]; %valid duration of fixations
validduration = Event.event_raw(:,fixdur_col)>= tw(1) & Event.event_raw(:,fixdur_col)<= tw(2);
% get first fixations
first =  Event.event_raw(:,firstpass_col)==1;
PPara.event_all = Event.event_raw(validduration&first,:);
eval(['PPara.event_all = PPara.event_all(' epoch_select ',:);']);
%%%% change the fixation onset trigger to fixation offset
if contains(RunCond,'Off')
    PPara.event_all(:,trig_col) = PPara.event_all(:,trig_col) + PPara.event_all(:,fixdur_col);
end
epoch = Get_Epoch(hdr,data,PPara,trig_col);

%%% save out 
save([PPath.SaveData 'epoch_' PPara.filename],'epoch','-v7.3');
disp(['*** epoching done! --- ' sub]);
clear data

end



