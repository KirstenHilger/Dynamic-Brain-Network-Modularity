%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Script written by: Kirsten Hilger & Olaf Sporns (2018)
%%         
%% Project: Dynamic Modularity and IQ   
%% Subject: Script a) determines the relation between static modularity
%%          and individual differences (e.g., in intelligence) 
%%          for all spatial resolution levels of gamma, b) finds optimal 
%%          spatial resolution level (gamma) for each subject 
%%          (via highest agreement between static NW and Yeo 17-NW  
%%          partition), c) computes variability (SD) in modularity over 
%%          time, d) counts number of extreme modularity states in two
%%          different ways, e) serves a node-specific measure of temporal
%%          brain network stability  via temporal coclassification,
%%          and f) calculates associations between these
%%          characteristics and an individual difference measures (here 
%%          intelligence) by controlling for effects of age, sex, hand  
%%          and mean framewise displacement (FD).  
%%          
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath(genpath(pwd))  % add this folder and all folder below to the path

%%
%%%%%%%%%%%%%%%%% ID-Analyses on static connectivity %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load results of community detection for static networks (community
% assignement vector CI)
load mod_data_static 

% load behavioral data file (here measure of interest, i.e., IQ = column
% 5, and variables of no interest are listed in columns 2, 3, 4, 11 (age,
% sex, hand, mean framewise displacement)
load Subjects_N281_Behav_Dyn.mat  

% compute association (spearman correlation, rho) between IQ and Q by 
% controlling for age, sex, hand, mean FD in STATIC networks for all 
% resolution levels (gamma)
[rho_Q p_Q] = partialcorr(unnamed(:,5),(Q),unnamed(:,[2:4 11]),'type','spearman')

% compute association (spearman correlation, rho) between IQ and numbers of 
% modules by controlling for age, sex, hand, mean FD in STATIC networks 
% for all resolution levels (gamma)
[rho_Mn p_Mn] = partialcorr(unnamed(:,5),(nM),unnamed(:,[2:4 11]),'type','spearman')

gam = [0.1:0.1:6]; G = length(gam); % specify gamma


%%
%%%%%%%%%%%% finding subject-specific optimal resolution levels %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load yeo 17 network standard parcellation (Yeo et al., 2011)
load Yeo17_NWidx_ex24.mat
P1 = Y17(:,2); P1n = max(P1);
P2 = Y17(:,1); P2n = max(P2); % 17 NW partition, because this is also the
% partition where our 114 node atlas comes from.

% compute similarity between CI (community identity) vector & Yeo partition
% for each subject and each gamma
for s=1:S
    disp(num2str(s));
    for g=1:G
        vi(s,g) = partition_distance(squeeze(CI(:,s,g)),P2);
    end;
end;

% detect minima of the difference between Yeo17 and ind. partition,
% note that the min here is not only the minimum over all 60 gammas, the
% minima detection bases now on the consensus of the gaussian-filtered 
% gamma(x-axis)-vi(y-axis) distribution with varying filters (1-5). 

filts = [1:1:5];
for s=1:S
    for f=1:length(filts)
        ivv = imgaussfilt(vi(s,:),filts(f));
        [yp xp] = findpeaks(-ivv,'MinPeakWidth',3);
        [myp mxp] = max(yp);
        if isempty(yp)
            opt_gam_Val2_temp(s,f) = NaN;
        else
            opt_gam_Val2_temp(s,f) = xp(mxp);
        end;
    end;
    opt_gam_Val2(s) = round(nanmedian(opt_gam_Val2_temp(s,:)));
end;

% saves vector of subject-specific optimal gammas (saved are indices not
% the gamma values)
save opt_gam_Val_Final opt_gam_Val2 


%%
%%%%%%%%%%%%%%%%% ID-Analyses on dynamic connectivity %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load results of community detection for dynamic networks (community
% assignement vector CI)
load ('mod_data_dynamic_70_tapered_Makoto')
load FCall70_tapered % load dynamic connectivity matrices
load Subjects_N281_Behav_Dyn.mat  % load behavioural data

g = opt_gam_Val2; % specify vector of subject-specific optimal gammas 
%(computed in the section above)
 
Q_optGam = zeros(281,70); % build empty data frame (subjects x windows) 

% calculate vector of modularity values for each subject corresponding to
% the optimal gamma
for s=1:S
    Q_optGam(s,:) = Q([1:1:70],s,g(s));
end    

Q_optGamT = Q_optGam'; % transpose vector

% compute association (spearman correlation, rho) between IQ and mean/std 
% of Q by controlling for age, sex, hand, mean FD in DYNAMIC networks
[rho_Q1 p_Q1] = partialcorr(unnamed(:,5),(mean(Q_optGamT))',unnamed(:,[2:4 11]),'type','spearman')
[rho_Q2 p_Q2] = partialcorr(unnamed(:,5),(std(Q_optGamT))',unnamed(:,[2:4 11]),'type','spearman')

% compute association (spearman correlation, rho) between IQ and mean/std 
% of "number of modules"  by controlling for age, sex, hand, mean FD in 
% DYNAMIC networks
[rho_nM1 p_nM1] = partialcorr(unnamed(:,5),(mean(nM_optGamT))',unnamed(:,[2:4 11]),'type','spearman')
[rho_Q2 p_Q2] = partialcorr(unnamed(:,5),(std(Q_optGamT))',unnamed(:,[2:4 11]),'type','spearman')


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identifies windows of low/high modularity episodes (states of over- or 
% disconnectivity) in two different ways a) subject-specific, i.e.,
% relative to individual mean Q over time (mean Q +/- 50% of mean Q) and b)
% general, i.e., relative to the group-averaged Q over time (mean Q +/-
% 50% of mean Q).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Way (a) subject specific:

%   provides count of instances (wins) where Q exceeds the
%   individual-specific upper threshold, i.e., mean Q + (e.g. 50% of mean Q)
for s=1:S
    qp = mean(Q_optGamT(:,s)) + (mean(Q_optGamT(:,s)) * 0.50);
    aboveVals = Q_optGamT(:,s)>=qp;
    SumQgp(s) = sum(aboveVals);
end 

%   provides count of instances (wins) where Q exceeds the
%   individual-specific lower threshold, i.e., mean Q - (e.g. 50% of mean Q)
for s=1:S
    qp = mean(Q_optGamT(:,s)) - (mean(Q_optGamT(:,s)) * 0.50);
    belowVals = Q_optGamT(:,s)<=qp;
    SumBQgp(s) = sum(belowVals);
end 

% relation between count of occurances and ID measure, here: intelligence
[rho_Qp2 p_Qp2] = partialcorr(unnamed(:,5), SumQgp(:),unnamed(:,[2:4 11]),'type','spearman')
[rho_Qp3 p_Qp3] = partialcorr(unnamed(:,5), SumBQgp(:),unnamed(:,[2:4 11]),'type','spearman')


%%%% Way (b) general (group-averaged):

for g=1:G
    disp(num2str(g)) % shows at which of 60 gam we are
    Qg = squeeze(Q(:,:,g)); % gamma-dependent Q-Threshold is determined as 
    %mean Q over all subjects and all windows at this specific gamma
    
    % high Q episodes
    prch(g) = prctile(Qg(:),75);
    ff = Qg>=prch(g);
    [x y] = find(ff);
    FCtemp = zeros(113);
    for i=1:length(x)
        FCtemp = squeeze(FCall(ind,ind,x(i),y(i))) + FCtemp;
    end;
    FCtemp = FCtemp./length(x);
    FCh(:,:,g) = FCtemp;
    
    % low Q episodes
    prcl(g) = prctile(Qg(:),25);
    ff = Qg<=prcl(g);
    [x y] = find(ff);
    FCtemp = zeros(113);
    for i=1:length(x)
        FCtemp = squeeze(FCall(ind,ind,x(i),y(i))) + FCtemp;
    end;
    FCtemp = FCtemp./length(x);
    FCl(:,:,g) = FCtemp;
    
    % all (mean over all windows & subjects)
    FCm(:,:,g) = squeeze(mean(mean(squeeze(FCall(ind,ind,:,:)),4),3));
end; % we got gamma-sepcifc states
gam_index = opt_gam_Val2;

% relation with intelligence (now we focus on individual optimal gamma)
for s=1:S
    nMg(:,s) = squeeze(nM(:,s,gam_index(s))); % number of modules
    Qg(:,s) = squeeze(Q(:,s,gam_index(s))); % modularity
    numhg(s) = sum(Qg(:,s)>=prch(gam_index(s)));  % number of high-Q windows
    numlg(s) = sum(Qg(:,s)<=prcl(gam_index(s)));  % number of low-Q windows
end;
    
[rho_G1_25 p_G1] = partialcorr(unnamed(:,5),numhg',unnamed(:,[2:4 11]),'type','spearman')
[rho_G2_25 p_G2] = partialcorr(unnamed(:,5),numlg',unnamed(:,[2:4 11]),'type','spearman')

% plotting group-mean time-averaged mean functional connectivity matrix at e.g., gamma 35 
imagesc(FCh(:,:,35)) % plot mean FC at e.g., gamma = 35 for the high-modularty states (average across all subjects and across all high-Q windows)
imagesc(FCm(:,:,35)) % plot mean FC at e.g., gamma = 35 for the modularty states that are not characterized of high/low modularity (average across all subjects and across all high-Q windows)
imagesc(FCl(:,:,35)) % plot mean FC at e.g., gamma = 35 for the low-modularty states (average across all subjects and across all low-Q windows)



%%
%%%%%%%% Variability in the modular partition itself (VI) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:S
    vidiff_ID(s,:) = vidiff(s,g(s)); % vidiff vector, i.e., mean inverse
    % mutual information criteria between different partitions over time
    % based on adjacent time windows
end    

% relation with ID-measure, e.g., intelligence
[rho_vi1 p_vi1] = partialcorr(unnamed(:,5),vidiff_ID,unnamed(:,[2:4 11]),'type','spearman')

for s=1:S
    vidiff2_ID(s,:) = vidiff2(s,g(s)); % vidiff vector, i.e., mean inverse
    % mutual information criteria between different partitions over time
    % based on all time windows
end   

% relation with ID-measure, e.g., intelligence
[rho_vi2 p_vi2] = partialcorr(unnamed(:,5),vidiff2_ID,unnamed(:,[2:4 11]),'type','spearman')




%%
%%%%%%%%% size of biggest module (measure of global integration) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find largest module for each sub in each win at each gam
clear M 

for s=1:S
    for t=1:T
        for g=1:G
            ci = squeeze(CI(:,t,s,g));
            ui = unique(ci);
            hi = histc(ci,ui);
            M(t,s,g) = mean(hi);
        end;
    end;
end;

gam_index = opt_gam_Val2;

% takes largest module of partition on optimal gamma 
for s=1:S
    Mg(:,s) = squeeze(M(:,s,gam_index(s)));
end;

% relation with intelligence
[rho_M1 p_M1] = partialcorr(unnamed(:,5),mean(Mg)',unnamed(:,[2:4 11]),'type','spearman')
[rho_M2 p_M2] = partialcorr(unnamed(:,5),std(Mg)',unnamed(:,[2:4 11]),'type','spearman')




%%
% node- or network-specific stability (coclassification) and Intelligence %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get agreement and compute block-wise averages

g = opt_gam_Val2; % vector of optimal gammas

Ag = zeros(N,N,S); 
Agm1 = zeros(P1n,P1n,S); % vector for aggregation within 7 network partition
Agm2 = zeros(P2n,P2n,S); % vector for aggregation within 17 network partition
%Agm3 = zeros(P3n,P3n,S);

for s=1:S % loop over subjects
    Ag(:,:,s) = agreement(CI(:,:,s,g(s)))./T;           % compute time-
    % average coclassification values at subject-specific gamma 
    % agreement = number of times any two nodes were assigned to the same 
    % module. Because we divide through T, values are standardized.  
    AgT4(:,:,s) = agreement(CI(:,50,s,g(s))); 
    % aggregate CC vals network wise in 3 network partitions:
    Agm1(:,:,s) = module_density(Ag(:,:,s),P1,1);       % Yeo 17
    Agm2(:,:,s) = module_density(Ag(:,:,s),P2,1);       % Yeo 7
    %Agm3(:,:,s) = module_density(Ag(:,:,s),P3,1);      % TPN/TNN
end;

% plotting of CC Values for one subject: 
imagesc(Ag(:,:,1))

 % plotting CC Values between timpe point 4 for one subjects
 imagesc(AgT4(:,:,1)) 

%%%%%% relation to IQ controlled for age, sex, handedness and mean FD

for n1=1:P1n % 7NW
    for n2=1:P1n
        [rho_RSN1(n1,n2) p_RSN1(n1,n2)] = partialcorr(unnamed(:,5),squeeze(Agm1(n1,n2,:)),unnamed(:,[2:4 11]),'type','spearman');
    end;
end;

rho_RSN1(:,:) % 7x7 matrix of correlation values between network-specific 
% stability and IQ association
p_RSN1(:,:) % corresponding p values

for n1=1:P2n % 17NW
    for n2=1:P2n
        [rho_RSN2(n1,n2) p_RSN2(n1,n2)] = partialcorr(unnamed(:,5),squeeze(Agm2(n1,n2,:)),unnamed(:,[2:4 11]),'type','spearman');
    end;
end;

rho_RSN2(:,:) % 17x17 matrix of correlation values between network-specific 
% stability and IQ association
p_RSN2(:,:) % corresponding p values

for n1=1:N % 113 Nodes (each node)
    for n2=1:N
        [rho_Nodes(n1,n2) p_Nodes(n1,n2)] = partialcorr(unnamed(:,5),squeeze(Ag(n1,n2,:)),unnamed(:,[2:4 11]),'type','spearman');
    end;
end;

for n1=1:P3n % TPN vs. TNN
    for n2=1:P3n
        [rho_RSN3(n1,n2) p_RSN3(n1,n2)] = partialcorr(unnamed(:,5),squeeze(Agm3(n1,n2,:)),unnamed(:,[2:4 11]),'type','spearman');
    end;
end;

rho_RSN3(:,:) % 2x2 matrix
p_RSN3(:,:) % p values




%%
%%%%%%%%%%%%%% Plooting high/low (1/3) IQ people %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low IQ people

k = load("Subjects_N281_Behav_Dyn.mat") % load beahavioral data

[a, IQ_indices_sorted] = sort(k.unnamed(:,5)) % sort in accordance to IQ

g = opt_gam_Val2; % specify vector of optimal gammas

Ag = zeros(N,N,94); % replace all but the upper third with zeros

for s=1:94 % now only 94 values of lower N = 94 IQ peaple
    new_s = IQ_indices_sorted(s);
    Ag(:,:,s) = agreement(CI(:,:,new_s,g(new_s)))./T;   % compute CC vals
end;

Ag_lowIQ = mean(Ag,3)

imagesc(Ag_lowIQ) % show mean coclassification matrix of low-IQ people


% high IQ people

k = load("Subjects_N281_Behav_Dyn.mat")

[a, IQ_indices_sorted] = sort(k.unnamed(:,5),'descend') % sort descend

g = opt_gam_Val2; % specify vector of optimal gammas

Ag = zeros(N,N,94); % replace all but the upper third with zeros


for s=1:94 % now only 94 values of high N = 94 IQ peaple
    new_s = IQ_indices_sorted(s);
    Ag(:,:,s) = agreement(CI(:,:,new_s,g(new_s)))./T;    % compute CC vals
end;

Ag_highIQ = mean(Ag,3)

imagesc(Ag_highIQ) % show mean coclassification matrix of high-IQ people

Ag_Diff = Ag_lowIQ - Ag_highIQ % calculate high - low-IQ difference 
%coclassification matrix

imagesc(Ag_Diff) % show mean coclassification matrix of high-IQ peaple 
%minus low-IQ people




%%
%%%%%%%%%%%%%% Plooting high/low (10%) IQ people %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low IQ people

k = load("Subjects_N281_Behav_Dyn.mat")

[a, IQ_indices_sorted] = sort(k.unnamed(:,5))

g = opt_gam_Val2; 

Ag = zeros(N,N,28); 


for s=1:28
    new_s = IQ_indices_sorted(s);
    Ag(:,:,s) = agreement(CI(:,:,new_s,g(new_s)))./T;           % compute CC vals
end;

Ag_lowIQ = mean(Ag,3)

imagesc(Ag_lowIQ)

imagesc(Ag_lowIQ(11:17,11:17))

imagesc(Ag_lowIQ(67:73,67:73))

% high IQ people

k = load("Subjects_N281_Behav_Dyn.mat")

[a, IQ_indices_sorted] = sort(k.unnamed(:,5),'descend')

g = opt_gam_Val2; 

Ag = zeros(N,N,28); 


for s=1:28
    new_s = IQ_indices_sorted(s);
    Ag(:,:,s) = agreement(CI(:,:,new_s,g(new_s)))./T;           % compute CC vals
end;

Ag_highIQ = mean(Ag,3)

imagesc(Ag_highIQ)

imagesc(Ag_highIQ(11:17,11:17))

imagesc(Ag_highIQ(67:73,67:73))


Ag_Diff = Ag_lowIQ - Ag_highIQ

imagesc(Ag_Diff)



















    
