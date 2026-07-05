% =========================================================================
% Basic waveform reading and plotting; XA, Yu (2026/07/01)
% =========================================================================
% >> Make sure enter the main interface
% CPATH='\ReadSAC';
CPATH=pwd;
addpath([CPATH,'\func']);
addpath([CPATH,'\ExampleData']);

% >> Prepare the Path and files
PATH=[CPATH,'\','ExampleData'];
FILE={
    'pws1.sac'
    'pws2.sac'
    'pws3.sac'
};
fullfiles=fullfile(PATH, FILE);
%% 
% >> Obtain results and plot
out=Fun_PlotWaves(fullfiles, ...
    'PeriodBand', [1 30], ...
    'FilterOrder', 4, ...
    'Normalize', true, ...
    'TimeXLim', [-150 150], ...
    'LineWidth',0.7);

rmpath([CPATH,'\func']);
rmpath([CPATH,'\ExampleData']);
