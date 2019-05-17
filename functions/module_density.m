function [Mden Mden_in Mden_bw Mden_in_lo Mden_bw_hi] = module_density(M,ci,flag)

% set flag=1 if main diagonal is to be discounted (assumes modules on main
% diagonal are square)

mnum = max(ci);
Mden = zeros(mnum);

for i=1:mnum
    iind = logical(ci==i);
    for j=1:mnum
        jind = logical(ci==j);
        Mij = M(iind,jind);
        if ((i==j)&&(flag==0))
            Mij = M(iind,jind);
        end;
        if ((i==j)&&(flag==1))
            Mij = M(iind,jind);
            mask = find(ones(size(Mij)).*~eye(size(Mij,1)));
            Mij = Mij(mask);
        end;
        Mden(i,j) = nanmean(Mij(:));
    end;
end;

Mden(isnan(Mden)) = 0;

m1mask = eye(mnum); ff1 = find(m1mask);
m2mask = ones(mnum).*~eye(mnum); ff2 = find(m2mask);
% mean density of within/between module blocks
% note: modules composed of a single node may have density 0
Mden_in = nanmean(Mden(ff1));
Mden_bw = nanmean(Mden(ff2));
% lowest within vs highest between
Mden_in_lo = min(nonzeros(Mden(ff1)));
Mden_bw_hi = max(Mden(ff2));


