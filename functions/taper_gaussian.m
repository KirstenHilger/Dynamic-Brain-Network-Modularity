function trFC = taper_gaussian(bold,width,sigma,Nstep)
% Compute time-resolve FC using tapered Gaussian window
%
% [Input]
%  bold: BOLD time series [nodes x time samples x subjects]
%  width: Width of rectangle to be convolved with Gaussian kernel in TR (= 66 in [1])
%  sigma: Sigma of Gaussian kernel in TR (= 9 in [1])
%  Nstep: Number of steps (skips) in TR (= 3 in [1])
%
% [Output]
%  trFC: Time-resolved FC [nodes x nodes x windows x subjects]
%
% [1] Fukushima M, Betzel RF, He Y, de Reus MA, van den Heuvel MP, Zuo XN, 
%     Sporns O, Fluctuations between high- and low-modularity topology in 
%     time-resolved functional connectivity, NeuroImage, 180, 406-416, 2018
% 
% September 3, 2018: M.Fukushima

% Number of nodes, time samples, subjects
[Nv,Nt,Nsub] = size(bold);

% Normalize BOLD
nbold = zeros(size(bold));
for ss = 1:Nsub
  for vv = 1:Nv
    nbold(vv,:,ss) = zscore(bold(vv,:,ss));
  end
end

%%% Compute Gaussian taperd window
f = @(m,s) exp(-0.5*m.^2/s)/sqrt(2*pi*s);
window_rect = [zeros(1,sigma-1) ones(1,width) zeros(1,sigma-1)];
arb = 3e3; % An arbitrary large number
window_temp = conv(window_rect,f(-arb:arb,sigma^2));
x = 1:(width+2*sigma-2); y = window_temp(x+arb);
% Resample weights of window to make window(central_sample)=max(window); may be unnecessary
xx = (1+0.5):((width+2*sigma-2)-0.5); yy = spline(x,y,xx);
window_taper = yy/sum(yy);
% Window length
Ntw = length(window_taper);

%%% Windowing
% Number of windows
Nw = length(1:Nstep:(Nt-Ntw+1));
trFC = zeros(Nv,Nv,Nw,Nsub);
for ss = 1:Nsub
  ww = 1;
  for tt = 1:Nstep:(Nt-Ntw+1)
    % Weighted correlation
    trFC(:,:,ww,ss) = weightedcorrs(nbold(:,tt:tt+Ntw-1,ss)',window_taper');
    ww = ww+1;
  end
end