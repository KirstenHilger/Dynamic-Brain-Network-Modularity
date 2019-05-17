%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Function written by: Kirsten Hilger (help by Joshua Faskowitz, 2018)
%%
%% Project: Dynamic Modulariyt and IQ  
%% Subject: Get connectivity matrices for all participants for static 
%%          network construction, i.e., edges based on time courses 
%%          averaged across all time points of preprocessed fMRI BOLD data.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ tvMats , windowsUsed, conMat ] = get_slide_window_conMats(inMat,windowSz,shiftSz)

% function to compute static or time-varying connectivity

% Input arguments: inMat = input Matrix (for each sub a matrix with nodes x
% raw Bold signal (num nodes x num TPs), as variable in a matlab structure); 
% windoSz = Windowsize in TP; shiftSz; WindowShiftSize in TP.

% Returns: tvMats = timevarying Matrices; the connectivity matrices for each window; 
% windowUsed = Size of the window finally used; conMats = static matrices
% connectivity matrix, i.e., 1 FC map for whole duration of scan, Pearson
% Correlation across all time points.

% These things the function returns won't be saved untill not explicitely
% metioned in the script that calles this function!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
   error('need at least three arguments') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nNode , nTp ] = size(inMat) ; % specifies the matrix for each subject

% first comupte the normalized connectivity matrix for whole data (conMat)
conMat = corr(inMat') ; % Input Matrix (nodes x time points)
% clear the diagonal
conMat(1:nNode+1:end) = NaN ; % replaces diagonal with NaNs

%%%%%%%%%%%%%%%%%%%%%%% specify window parameter %%%%%%%%%%%%%%%%%%%%%%%

if windowSz > nTp % if WindowSize < number total timepoints --> Error!
   error('window size bigger than number timepoints') 
end

% calculate how many tvMats we will calculate
winInds = (0:shiftSz:nTp)' + (1:windowSz) ;
% calculate rows that exceed the total nTp (length of whole data in TPs)
% and sets these ti 0.
valid_winInds = winInds(sum(winInds > nTp,2) == 0,:) ;

% num windows, and preallocate the output matrix
nWin = size(valid_winInds,1) ;
tvMats = zeros([nNode nNode nWin]) ; % creates tvMats as empty matrix

% loop through each window
for idx = 1:nWin

    disp(idx) % show in which window we are actually
    
    % the data to work on 
    tmpMat = inMat(:,valid_winInds(idx,:)) ;
    
    % get the correlation
    tmpCon = corr(tmpMat') ;
    tmpCon(1:nNode+1:end) = NaN ;
    tvMats(:,:,idx) = tmpCon ;
    
end

windowsUsed = [ valid_winInds(:,1) valid_winInds(:,end) ] ; % returns also 
% shape of the final window used for all analyses. After the loop, because
% it is the same for all windows (and later on also for all Subs).
