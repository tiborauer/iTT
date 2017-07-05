function [Z, XYZ, th] = spm_uc_iTT(Z,XYZ,p,iter)
% computes thresholds using two-theshold (TT) correction and apply them
% FORMAT [Z, tau, XYZ] = spm_uc_TT(Z,p,iter,XYZ);
% 
% Required input:
%   Z    - Statistics {unfiltered} 
%   XYZ  - location of voxels {voxel coords}
%   p    - upper probability threshold {1-p}
%   iter - determines whether iterative optimization is to be employed
%
% Provided output:
%   Z   - Statistics {filtered with TT}
%   tau - upper threshold in measure of the statistical parameter
%   XYZ - location of voxels {voxel coords} {filtered with TT}
%_______________________________________________________________________
% Copyright (C) 2011 Biomedizinische NMR Forschungs GmbH

% Tibor Auer
% $Id: spm_uc_TT.m 2011-07-05 $ 

%-calculating thresholds
%--------------------------------------------------------------------------
th = th_comp(Z, XYZ, p, iter);

%-applying thresholds
%--------------------------------------------------------------------------
[Z XYZ] = th_app(Z,XYZ,th);

end

%==========================================================================
%-Computation of two-threshold
%==========================================================================

function [th stat] = th_comp(Z0, XYZ, lt, iter)
% computes the data-driven thresholds
% FORMAT th = comp(Z, ut, iter);
%
% Required inputs:
%   Z       - Statistics {unfiltered}
%   XYZ     - location of voxels {voxel coords}
%   lt      - upper probability threshold {p} (e.g. 0.05)
%   iter    - determines whether iterative optimization is to be employed
%
% Fields of output th:
%   .tal    - lower threshold of activation in the measure of the parameter
%   .tau    - upper threshold of activation in the measure of the parameter
%   .tdl    - lower threshold of deactivation in the measure of the parameter
%   .tdu    - upper threshold of deactivation in the measure of the parameter
%

%-intitializing
%--------------------------------------------------------------------------
iTT = IniFile(fullfile(spm('dir'),'toolbox','iTT','config.ini'));
try
    ut = 1-(iTT.thresholds.ut/2); % retrieving upper probability threshold
    SHOW = iTT.config.show;
catch % iTT configuration file is not available
    ut = 0.99995;
    SHOW = true;    
end
SKIP = 0.25;
lt = 1 - (lt/2);
stat.cuts = [0.2 0.8];
main = false;
if check_licence({'fmincon'}), disp('Optimisation Toolbox is not available!'); end
iter = iter && ~check_licence({'fmincon'});

%-histogram
%--------------------------------------------------------------------------
Z0 = round(Z0*1000)/1000; Z0(isnan(Z0)) = 0;
Z = remove_peak(Z0);

nbin = ceil(sqrt(numel(Z)));
go = true;
while go
    [f0, xc] = hist(Z,nbin);
    if any(abs(diff(f0)) > max(f0)/3)
        nbin = ceil(nbin*0.9);
    else
        go = false;
    end
end
% [f0, xc] = hist(Z,-20:0.02:20);
f0 = f0'; xc = xc';

%-deteriming the place of the main histogram
%--------------------------------------------------------------------------
    q25=round(numel(f0)*SKIP); q75=round(numel(f0)*(1-SKIP));
    pf = max(f0(q25:q75)); indpf=round(mean(find(f0(q25:q75)==pf)))+q25-1;

if main    
    STEP=round(numel(f0)/10);
    d=(mean(f0(indpf:indpf+STEP))-pf)/abs(mean(f0(indpf:indpf+STEP))-pf);
    
    ind=indpf; max0 = pf;
    if d > 0
        while max(f0(ind:ind+STEP)) > max0
            max0 = max(f0(ind:ind+STEP));
            ind = ind + STEP;
        end
        pf = max(f0(ind-STEP:ind));
        indpf=round(mean(find(f0(ind-STEP:ind)==pf)))+ind-STEP-1;
    else
        while max(f0(ind-STEP:ind)) > max0
            max0 = max(f0(ind-STEP:ind));
            ind = ind - STEP;
            if (ind <= STEP), STEP=ind-1; end
        end
        STEP=round(numel(f0)/10);
        pf = max(f0(ind:ind+STEP));
        indpf=round(mean(find(f0(ind:ind+STEP)==pf)))+ind-1;
    end
    
    [fh xh] = hist(f0(1:indpf),numel(f0(1:indpf))); fh = fh(1:floor(numel(fh)/2)); indh_l = max(xh(fh==max(fh)));
    [fh xh] = hist(f0(indpf:end),numel(f0(indpf:end))); fh = fh(1:floor(numel(fh)/2)); indh_r = max(xh(fh==max(fh)));
    i_l = []; i_r = [];
    for i = 0:indpf-1
        if isempty(i_l) && (f0(indpf-i) <= indh_l), i_l = indpf-i ; end
        if isempty(i_r) && (f0(indpf+i) <= indh_r), i_r = indpf+i ; end
    end
    if isempty(i_l), i_l = 1; end
    if isempty(i_r), i_r = 2*indpf; end
    
    main_f00 = f0(i_l:i_r);
    main_xc0 = xc(i_l:i_r);
    
    for r = 1:10
        main_xc = interp(main_xc0,r);
        main_xc = [min(main_xc0); main_xc(isin(main_xc,min(main_xc0),max(main_xc0))); max(main_xc0)];
        main_f0 = pchip(main_xc0,main_f00,main_xc);
        ff = hist(main_f0,0:pf*0.1:pf); ff = ff(2:end);
        OK = true;
        for i = 1:numel(ff)-1
            OK =  OK && (sum(ff(i:i+1)) >= 3);
        end
        if OK, break, end
    end
else
    main_xc = xc;
    main_f0 = f0;
end
% plot(main_xc0, main_f00), hold on, plot(main_xc, main_f0,'r');

%-Determinig cuts
%--------------------------------------------------------------------------
if iter
    stat.cuts = fmincon(@(x)calc_fit(x,main_xc,main_f0,xc,f0,lt),stat.cuts,[],[],[],[],[0.1 0.5],[0.5 0.9],@cut_con,optimset('Algorithm','interior-point'));
end

%-Gaussian-fit, calculating the cumulative distribution form the
%Gaussian curve and determining threshold
%--------------------------------------------------------------------------
fres = Gauss_fit(main_xc(main_xc<=0.5), main_f0(main_xc<=0.5), stat.cuts);
G = feval(fres, xc);
inty = cumsum(G);
th.tal = xc(find(inty>=(lt*max(inty)),1));
th.tdl = xc(find(inty>=((1-lt)*max(inty)),1));
med = find(xc>mean([th.tal th.tdl]),1);
i1 = find(f0<(stat.cuts(2)*max(main_f0))); i1 = i1(i1>med); i1 = i1(1);
i2 = find(xc==th.tal);
i1 = round((i1+i2)/2);
Tol = std(f0(i1:i2)-G(i1:i2));
md_sig = abs(max(f0(i2:end)-G(i2:end)));
if iter
    if md_sig <= Tol
        STEP = 0.001;
        for i = max(Z):-STEP:th.tal
            if max(abs(calc_ut(i,Z0,XYZ,lt,fres,xc,0))) > Tol, break; end
        end
        th.tau = i + STEP;
    else
        [th.tau fval exitflag output] = fminbnd(@(x)mean(abs(calc_ut(x,Z0,XYZ,lt,fres,xc,0))),th.tal,max(Z));
    end
    median=mean([th.tal th.tdl]);
    th.tdu = median-(th.tau-median);
else
    th.tau = xc(find(inty>=(ut*max(inty)),1));
    th.tdu = xc(find(inty>=((1-ut)*max(inty)),1));
end

%-Calculating sub-histogram
%--------------------------------------------------------------------------
Z = th_app(Z0,XYZ,th); Z0(Z>0) = 0;
Z = Z0(Z0~=0);
[f1, xc] = hist(Z,xc); f1 = f1'; xc = xc'; i1=find(xc==th.tal);
G = feval(fres, xc); eG = G(i1:end); eF0 = f0(i1:end); eF1 = f1(i1:end);
P = curve_sub(eF0,eF1); GP = curve_sub(eF0,eG);

stat.P = P; stat.GP = GP;
stat.FP = curve_sub(eG,eF1)-curve_sub(eG,eF0);
stat.TP = (P - stat.FP);

stat.FPR = stat.FP/sum(G);
stat.TPR = stat.TP/GP;

%-popping up a figure containing the histogram tha Gaussian-fit and the
%thresholds calculated
%--------------------------------------------------------------------------
if SHOW
    f = figure(100);
    set(f, 'Name', 'Distribution of the given statistical parameter');
    set(f, 'NumberTitle', 'off');
    clf;

    xlim([min(xc) max(xc)]);
    hold on;
    h = plot(xc,f0,'w.');
    set(h,'MarkerEdgeColor','b');
    set(h, 'MarkerSize',5);
    
    h = plot(xc,f1,'w.');
    set(h,'MarkerEdgeColor','g');
    set(h, 'MarkerSize',5);
    set(h, 'LineWidth',2);
    
    h = plot(fres,'r');
    set(h, 'LineWidth',2);
    
    plot(xc,ones(1,numel(xc))*pf*stat.cuts(1),'-k');
    plot(xc,ones(1,numel(xc))*pf*stat.cuts(2),'-k');
    set(gca,'ytick',pf*stat.cuts);
    set(gca,'yticklabel', num2str(stat.cuts'*100,3));
    
    text = {'original data' 'data  after removal' 'fitting' '"80%" line' '"30%" line'};
    legend(text);
    
    ths = [th.tdu, th.tdl, th.tal, th.tau];
    yl = ylim*1.5;
    plot(ones(1,yl(2)+1)*th.tdu,yl(1):yl(2),'-k');
    plot(ones(1,yl(2)+1)*th.tau,yl(1):yl(2),'-k');
    plot(ones(1,yl(2)+1)*th.tdl,yl(1):yl(2),'-k');
    plot(ones(1,yl(2)+1)*th.tal,yl(1):yl(2),'-k');
    ylim(ylim/1.5);
    set(gca,'xtick',ths);
    set(gca,'xticklabel', num2str(ths',3));
    hold off;
end

end


function fres = Gauss_fit(x, y, cuts)
% 0.3 and 0.8 determine the central part of the distribution
%----------------------------------------------------------------------
pf = max(y); ind = [];
for i = 1:numel(y)
    if ~isin(y(i),pf*cuts(1),pf*cuts(2))
        ind(end+1) = i;
    end
end
outliers = excludedata(x,y,'indices',ind);
%-Gaussian-fit
%--------------------------------------------------------------------------
ftype = fittype('gauss1');
opt = fitoptions(ftype);
% opt = fitoptions(opt,'StartPoint',[max(y),5,max(y)],'Exclude',outliers);
%    'Lower', [max(f0)/2,-10,0], 'Upper', [max(f0)*2,10,max(f0)]
try
    fres = fit(x,y,ftype,opt);
catch
    opt = fitoptions(opt,'StartPoint',[max(y),5,1],'Exclude',outliers);
    fres = fit(x,y,ftype,opt);
end
end


%==========================================================================
%-Iterative functions
%==========================================================================

function diff = calc_fit(cuts,main_xc,main_f0,xc,f0,lt)
pf = max(main_f0);
fres = Gauss_fit(main_xc, main_f0, cuts);
inty = cumsum(feval(fres, xc));
tal = xc(find(inty>=(lt*max(inty)),1));

diff = abs(f0(xc==tal)-pf*cuts(1))+abs(feval(fres,tal)-pf*cuts(1));
end

function [c ceq] = cut_con(x)
TOL=0.2;
c = x(1)+TOL-x(2);
ceq = [];
end

function diff = calc_ut(ut,Z0,XYZ,lt,fres,xc,Tol)
y0 = feval(fres,xc);
inty = cumsum(y0);
th.tal = xc(find(inty>=(lt*max(inty)),1));
th.tdl = xc(find(inty>=((1-lt)*max(inty)),1));
median=mean([th.tal th.tdl]);
th.tau = ut;
th.tdu = median-(ut-median);
Z = th_app(Z0,XYZ,th);

Z0(Z~=0) = 0;
Z = Z0(Z0~=0);
[f0, xc] = hist(Z,xc); f0 = f0'; xc = xc';
i1=find(xc==th.tal);
diff = y0(i1:end)+Tol-f0(i1:end);
end

%==========================================================================
%-Applying two-threshold
%==========================================================================

function [Z XYZ]= th_app(Z,XYZ,th)
%-filtering out any voxel with Z below lower threshold
%--------------------------------------------------------------------------
Q      = find(Z > th.tal);
Z      = Z(:,Q);
XYZ    = XYZ(:,Q); 
%-finding clusters contaning voxel with Z greater then upper threshold
%--------------------------------------------------------------------------
A     = spm_clusters(XYZ);
Q     = [];
for i = 1:max(A)
    j = find(A == i);
    if any(Z(j) > th.tau); Q = [Q j]; end
end
Z     = Z(:,Q);
XYZ   = XYZ(:,Q);
end

%==========================================================================
%-Tool
%==========================================================================

function out = remove_peak(in)
% removes values with outlyingly high frequency
% FORMAT out = remove_peak(in)
%
% Required input:
%   in  - data with maximum precision of 0.0001

d = round(in(:)*1000);
mincc = min(d);
maxcc = max(d);
xc = double((mincc:maxcc)');
R = maxcc-mincc;
f = zeros(R+1,1);
for i = 1:length(d)
    f(d(i)-mincc+1) = f(d(i)-mincc+1) + 1;
end
ind = f~=0; f = f(ind); xc = xc(ind);

go = true;

while go
    th = mean(f)+3*std(f);
    for i = 1:numel(xc)
        if f(i) > th
            d(d==xc(i)) = NaN;
            f(i) = NaN;
            xc(i) = NaN;
        end
    end
    if any(isnan(d))
        d = d(~isnan(d));
        f = f(~isnan(f));
        xc = xc(~isnan(xc));
    else
        go = false;
    end
end

out = d/1000;

end

function res = isin(n, nl, nh)
% returns whether a given test value is in a given range
% FORMAT res = isin(n, nl, nh);
%
% Required input:
%   n   - test value
%   nl  - lower end of the range
%   nh  - higher end of the range

for i = 1:numel(n)
    res(i) = (n(i) >= nl) && (n(i) <= nh);
end
end

function cd = curve_sub(c1,c2)
i1d = c1 > c2;
ad = c1 - c2;
cd = sum(ad(i1d));
end