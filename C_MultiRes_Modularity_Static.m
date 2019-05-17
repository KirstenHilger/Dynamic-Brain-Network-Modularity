%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Script written by: Kirsten Hilger & Olaf Sporns (2018)
%%         
%% Project: Dynamic Modularity and IQ   
%% Subject: Determines for each subject's static connectivity matrix the  
%%          "best" modular partition on different spatial resolutions.   
%%          The script provides also average module size and average number 
%%          of modules. All metrics are computed across a full range of 
%%          different spatial reslutions (gamma) and based on the Louvain
%%          algorithm (Blondel et al., 2008).
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath(genpath(pwd)) % add this folder and all folder below to the path

% load 3D matrix containing all subjects static connectivity matrices 
load FCstatic_all   % (size nodes x nodes x subjects)
FC = squeeze(FC);

ind = setdiff(1:114,24); % exclude one node (index 24) in all matrices as
% one subject had missing values here
R = 100; % runs of Louvain, recommend at least 100
N = 113; % number of nodes
S = 281; % number of subjects

tau = 0.1; % for the consensus step over the Louvain Runs (R)
reps = 10; % for the consensus step over the Louvain Runs (R)


%specify range of spatial resolutions (here 0.1 - 6) and steps (here 0.1)
gam = [0.1:0.1:6]; G = length(gam); 
tic % timing

%% use parallel to make computation more efficient (computationally intense)
parfor s=1:S % loop over subjects)
    
    disp(num2str(s)); % display actual subject working on 
            
        FCi = squeeze(FC(:,:,s)); % squeeze FCi, i.e., subject-specific 
        % static connectivity matrix
        FCi = (FCi+FCi')./2; % use only half of the matrix
        FCi = FCi(ind,ind); % indices
        FCi(isnan(FCi)) = 0; % nans = 0
        FCi = fisherZTransform(FCi); % fisher Z transform all correlations
        
        for g=1:G % loop over range of gamma
            
            gamma = gam(g); % resolution parameter gamma
            
            qmax = -10;
            ciall = zeros(N,R);
            
            for r=1:R % loop over runs (of Louvain)
                % perform community detection on each FCi and each gamma
                [ci q] = community_louvain(FCi,gamma,[],'negative_asym');
                if(q>qmax)
                    qmax = q;
                    CI(:,s,g) = ci; % community vector based on partition 
                    % that maximizes Q
                    Q(s,g) = q;
                end;
                ciall(:,r) = ci;
                
            end;
            CI2(:,s,g) = consensus_und(agreement(ciall),tau,reps); %builds 
            % consensus over all Loauvain runs.
            
        end;
    
end;

toc

% note, CI ist the community vector based on the (out of 100 runs Louvain)
% partition where Q is maximal. CI2 is the community vector based on an
% agreement partition (over all 100 partitions). However, both serve nearly 
% identical results. We suggest to rely on CI for simplicity.

%% compute numbers of modules and avg. module size
for s=1:S % loop over subjects
    
    for g=1:G % loop over gamma
        
        ua = unique(CI(:,s,g));
        ub = histc(CI(:,s,g),ua);
        nM(s,g) = sum(ub>1); % number of modules
        mM(s,g) = mean(ub(ub>1)); % size of modules
        
    end;
    
end;

%% saves the modular partitions with all computed parameters into one
% matlab structure named "mod_data_static.mat".

save mod_data_static CI CI2 Q nM mM ind R N S G gam tau reps

