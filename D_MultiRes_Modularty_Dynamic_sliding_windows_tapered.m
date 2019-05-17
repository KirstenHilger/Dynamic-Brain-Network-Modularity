%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Script written by: Kirsten Hilger & Olaf Sporns (2018)
%%         
%% Project: Dynamic Modularity and IQ   
%% Subject: Determines for each subject and each connectivity matrix 
%%          (based on sliding tapered windows) the "best" modular
%%          bain network partition on different spatial resolutions.   
%%          The script provides also average module size and average number 
%%          of modules. All metrics are computed across a full range of 
%%          different spatial reslutions (gamma) and based on the Louvain
%%          algorithm (Blondel et al., 2008).
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath(genpath(pwd)) % add this folder and all folder below to the path

% load 4D matrix containing all subjects dynamic connectivity matrices (one
% matrix for each window). 4D matrix: nodes x nodes x windows x subjects
load FCall70_tapered 

ind = setdiff(1:114,24);

R = 100;   % Runs of louvain
N = 113;   % Number of nodes
S = 281;   % Number of subjects
T = 70;    % Number of time windows

%%%%%%%%%% change this for lower temporal resolution %%%%%%%%%%%%%%%%%%%%%%
T = 70;                         % number of time-varying Windows  
win = 1:1:T;                    % change to e.g. to 1:4:73 to get 19 Wins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = 0.1;   % for the consensus step over the Louvain Runs (R)
reps = 10;   % for the consensus step over the Louvain Runs (R)

%specify range of spatial resolutions (here 0.1 - 6) and steps (here 0.1)
gam = [0.1:0.1:6]; G = length(gam);  

% extract FC windows
FC = squeeze(FCall(:,:,win,:));  clear FCall

tic % timing

%% use parallel to make computation more efficient (computationally intense
parfor s=1:S   % loop over subjects
    
    disp(num2str(s))
    
    for t=1:T   % loop over windows
        
        FCi = squeeze(FC(:,:,t,s));  % squeeze FCi, i.e., subject-specific
        % connectivity matrix
        FCi = (FCi+FCi')./2; % use only half of the matrix
        FCi = FCi(ind,ind); % indices
        FCi(isnan(FCi)) = 0; % nans = 0
        FCi = fisherZTransform(FCi); % fisher Z transform all correlations
        
        for g=1:G   % loop over range of gamma
            
            gamma = gam(g); % resolution parameter gamma
            
            qmax = -10;
            ciall = zeros(N,R);
            
            for r=1:R   % loop over runs (of Louvain)
                % perform community detection on each FCi and each gamma
                [ci q] = community_louvain(FCi,gamma,[],'negative_asym');
                if(q>qmax)
                    qmax = q;
                    CI(:,t,s,g) = ci; % community vector based on partition
                    % that maximizes Q
                    Q(t,s,g) = q;
                end;
                ciall(:,r) = ci;
                
            end;
            CI2(:,t,s,g) = consensus_und(agreement(ciall),tau,reps);%builds 
            % consensus over all Loauvain runs.
            
        end;
                
    end;
    
end;

toc

% note, CI ist the community vector based on the (out of 100 runs Louvain)
% partition where Q is maximal. CI2 is the community vector based on an
% agreement partition (over all 100 partitions). However, both serve nearly 
% identical results. We suggest to rely on CI for simplicity.



%% 
% compute numbers of modules and avg. module size. Done for each sub and 
% each gamma
for s=1:S % loop over subjects
    
    for t=1:T % loop over windows
        
        for g=1:G % loop over gamma
            
            ua = unique(CI2(:,t,s,g));
            ub = histc(CI2(:,t,s,g),ua);
            nM(t,s,g) = sum(ub>1);  % number of modules
            mM(t,s,g) = mean(ub(ub>1)); % size of modules
            
        end;
        
    end;
    
end;

toc


%% 
% compute variability (operationalized as mean inverse mutual information
% criteria (vi) and mean mutual information criteria (mi) between 
% modular partitions over time (over windows). Done for each subjects and 
% each gamma in two different ways, i.e., based on adjacent time windows
% and based on all time windows.

for s=1:S % loop over subjects
    s
    for g=1:G % loop over gamma
        ci = squeeze(CI(:,:,s,g));
        for t=1:T-1
            [vi(t) mi(t)] = partition_distance(ci(:,t),ci(:,t+1));
        end;
        vidiff(s,g) = mean(vi); % the mean vi between adjacent time slices
        midiff(s,g) = mean(mi);
        cnt = 1;
        for t1=1:T-1
            for t2=t1+1:T
                [vi2(cnt) mi2(cnt)] = partition_distance(ci(:,t1),ci(:,t2));
                cnt = cnt+1;
            end;
        end;
        vidiff2(s,g) = mean(vi2); % the mean vi between all T time slices
        midiff2(s,g) = mean(mi2);
    end;
end;

%% 
% saves the modular partitions with all computed parameters into one
% matrix named "mod_data_dynamic_70Wins_tapered_makoto.mat".

save mod_data_dynamic_70Wins_tapered_makoto CI CI2 Q nM mM ind R N S T G gam tau reps vidiff vidiff2

