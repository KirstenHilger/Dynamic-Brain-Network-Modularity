%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Script written by: Kirsten Hilger (help by Joshua Faskowitz, 2018)
%%
%% Project: Dynamic Modularity and IQ  
%% Subject: Get connectivity matrices for all participants for static 
%%          network construction, i.e., edges based on time courses 
%%          averaged across all time points of preprocessed fMRI BOLD data.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Script needs the function get_slide_window_conMats.m located in 
% same folder or below. 

clear all
close all

addpath(genpath(pwd)) % add this folder and folders below

% Load matlab structure containing subject x scanning sequences (3) on the
% highest level and nodes x time series for each subject and for each 
% scanning sequence on the lower level. Structure is located in 00_Data.

loadedData = load('timeseries-gs.Yeo2011.mm316_281_renamed_KH.mat') ; 
                 
numberOfSubject = size(loadedData.subjects,2) ; % size = all subjects


for idx = 1:numberOfSubject % loop over all subjects

    disp(idx) % displays for which subjects we do it actually
    
    tmpData = loadedData.subjects(idx).mx645 ; % consider only the column 
    % mx645 in the lower-level matlab structure loaded above.  It as is the 
    % high-resolution fMRI data we used in the actual study. 
    tmpName = loadedData.subjects(idx).id ; % name (index) of subject
    
    [~,~,conMat] = get_slide_window_conMats(tmpData,156,39) ; % calls 
    % function that computed actually the matrices
    
    % When you want to save not the connectivity matrices but other things 
    % that the function returns (for example the conMATs for each window 
    % additionally to the static connectivity matrix) you just have
    % to replace the first tilde with a name and change the name below for
    % saving all stuff. We used this function (get_slide_window_conMats)
    % here only for getting the static connectivity matrices. Because 
    % the function used for the dynamic matrices
    % B_FCmatrices_Dynamic_sliding_windows_tapered.m does also the
    % window tapering. This function here provides untapered windows.
    
    % But generally, get_slide_window_conMat.m can also be 
    % used for computing dynamic matrices (matrices for each window). 
    
    % You have to provide as input for the function: Input time series data
    % (loeaded above), WindowSize (in TP), WindowShift (in TP)

    outName = strcat(tmpName,'_FC_Mat_EnhNKI_TR645_Window_staticFC.mat') ;
    save(outName,'conMat') % saves resulting matrices for every subjects
    
    % for further analyses, we recommend to combine these matrices (2D) 
    % into one single matrix file (3D) called 'FCstatic_all.mat'. 
    % This woul best fit the other scripts for further analyzing
    % modularity etc.
    
end

% note: to display data as graphical map: imagesc(tmpOut(:,:,nW)
% replace nW with the number of the window you want to see.

