% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise

function [hdr,data,Trigger_MEG] = Get_MEGData(PPath,File)
nfif = length(File.MEG);
data = [];
Trigger_MEG = [];
meg_len = [];
for p = 1:nfif
    cfg.dataset = [PPath.RawMEG File.MEG{p}];
    %%% read in multiple MEG datasets
    data_tmp = ft_read_data(cfg.dataset);
    event = ft_read_event(cfg.dataset);
    trig_tmp = [[event(strcmp('Trigger',{event.type})).value]' [event(strcmp('Trigger',{event.type})).sample]'];
    % concatenating multiple datasets
    add_last_len = sum(meg_len); % adding the data_length from the last .fif file
    meg_len = [meg_len size(data_tmp,2)];
    data = [data data_tmp];
    Trigger_MEG = [Trigger_MEG; [trig_tmp(:,1) trig_tmp(:,2)+add_last_len]];
end
hdr =  ft_read_header([PPath.RawMEG File.MEG{1}]);

