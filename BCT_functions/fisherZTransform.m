
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coded by: Olaf Sporns  
%%% Date: 28.09.2018
%%% Platform: MATLAB
%%% Purpose: Definition of function for Fisher-Z tranformation 
	% (here: of all correlations in Adj matrix)
%%% Note: Script has to be added when running Louvain_modules.m		 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z=fisherZTransform(C) 	% Input is the correlation matrix C
	
% fisherZTransform Compute Fisher z-transform of correlation matrix
%
% Input Arguments
%__________________________________________________________________________
%
%   C -- Correlation matrix
%
% Output Arguments
%__________________________________________________________________________
%
%   Z -- Z-score matrix (arctanh(C)) with zeroed diagonal

Z=atanh(C);    					% Output is Z, the transformed matrix
for i=1:length(Z)
    Z(i,i)=0;
end
end
