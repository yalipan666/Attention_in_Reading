%%% set paths
function [rootdir,cmap] = Get_Paths(server)
if server
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20220208/
    ft_defaults
    rootdir = '/rds/projects/2018/jenseno-reading/Full_Attention/';
else
    addpath Z:\fieldtrip-20220208\
    ft_defaults
    rootdir = 'Z:\Full_Attention\';
end
addpath(genpath([rootdir 'Analyse_codes' filesep]))
cmap = colormap(cbrewer('div','RdBu',32));
cmap = colormap(flipud(cmap)); close;
end