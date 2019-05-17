%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Script written by: Kirsten Hilger & Olaf Sporns (2018)
%%
%% Project: Dynamic Modularity and IQ  
%% Subject: Compute connectivity matrices for all 
%%          windows of specific size that slides (with a specifc shift) 
%%          through the matrix of nodes x time poins (
%%          preprocessed raw BOLD signal data). Windows are tapered with 
%%          Gaussian kernel and triangle (see Fukushima et al., 2017) to 
%%          increase sensitivity.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath(genpath(pwd)) % add this folder and allfolders below to the path

% Load matlab structure containing subject x scanning sequences (3) on the
% highest level and nodes x time series for each subject and for each 
% scanning sequence on the lower level. Structure is located in 00_Data.

load timeseries-gs.Yeo2011.mm316_281_renamed_KH.mat

S = 281; % number of subjects

% reorder the time series by subject ID and put them into one array
ids = char(subjects.id);
ids_num = str2num(ids(:,end-4:end));
[aa bb] = sort(ids_num,'ascend');

% build empty data frame of size: nodes x time points x subjects
ts = zeros(114,885,281); 

% loop over all subjects and fill the data frame
for s=1:S
    y = subjects(bb(s)).mx645;
    yt = size(y,2);
    ts(:,1:size(y,2),s) = y;
end;

W = 156;  % window length
dt = 10;  % window offset (shift in time points)
run_length = 885;  % length of entire run (in time points)
offsets = 1:dt:885-W+1;   % offsets (starting points for the windows)
sigma = 21; % tapering parameter

%%%% create windows

% returns FCall, i.e., 4D matrix with: nodes x nodes x windows x subjects
% FCall contains the connectivity matrices of all subjects for all windows

FCall = zeros(114,114,length(offsets),S);

for s=1:S
    ts_temp = squeeze(ts(:,:,s));
    FCall = taper_gaussian(ts,W,sigma,dt); % calls the tapering function 
end;

% save FCall under specified name +  window length (W), window
% offset (dt), length of the entire run, and the starting point (offsets)
save FCall70_tapered FCall W dt run_length offsets

% nice function for plotting connectivity matrices in a video
for i=1:73; imagesc(squeeze(fcall(:,:,i,1)),[-0.8 0.8]); pause(0.1); drawnow; end;


