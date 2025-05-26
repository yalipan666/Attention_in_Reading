% copy from Lexical/Analyse_codes
% 20210719 clean functions to make scripts more concise

function comp = ICA4rawdata(epoch4ICA)
%%% ICA downsampling 
tic
cfg            = [];
cfg.resamplefs = 200;
cfg.detrend    = 'no';
data_train = ft_resampledata(cfg, epoch4ICA);
toc
%%% ICA decomposition
tic
cfg                 = [];
cfg.method          = 'runica';
cfg.runica.maxsteps = 100;
comp                = ft_componentanalysis(cfg,data_train);
toc
end

