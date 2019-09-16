function [rsq2, pvals,resids,slope,intercept,R2Total] = ...
    thresholdplotter(thresholdval, xdata, ydata, lstyle, suppressplot)

% Use stats_regress2line to get threshold

hi = (xdata >= thresholdval);
lo = (xdata < thresholdval);
resids = [];
slope = [];
intercept = []; 

%Below threshold line
if sum(lo) ~= 0
    %lower equation
    lm1 = fitlm(xdata(lo), ydata(lo));
    resids = [resids; lm1.Residuals{:,1}];
    intercept = [intercept; lm1.Coefficients.Estimate(1)];
    slope = [slope; lm1.Coefficients.Estimate(2)];
    
    
    rsq2.lo = lm1.Rsquared.Ordinary;
    lmsum = anova(lm1,'summary');
    pvals.lo = lmsum{2,5}; 
    
%     p_lo = polyfit(xdata(lo), ydata(lo),1);
%     yfit_lo = polyval(p_lo,xdata(lo));
%     rsq2.lo = rsquare(ydata(lo), yfit_lo); 

    if suppressplot == 0 
        hold on
        m1 = lm1.Coefficients{2,1};
        b1 = lm1.Coefficients{1,1};
        t(1) = plot([min(xdata(lo)),max(xdata(lo))], ...
            [m1*min(xdata(lo))+b1, m1*max(xdata(lo))+b1],...
            'linewidth',1.5,'color','r', 'linestyle', lstyle);
    else
        t(1) = 0; 
    end
else 
    rsq2.lo = NaN;
end

%Above threshold line
if sum(hi) ~= 0
    %upper equation
    lm2 = fitlm(xdata(hi), ydata(hi));
    resids = [resids; lm2.Residuals{:,1}];
    intercept = [intercept; lm2.Coefficients.Estimate(1)];
    slope = [slope; lm2.Coefficients.Estimate(2)];
    
    rsq2.hi = lm2.Rsquared.Ordinary;
    lmsum = anova(lm2,'summary');
    pvals.hi = lmsum{2,5};
    
%     p_hi = polyfit(xdata(hi), ydata(hi),1);
%     yfit_hi = polyval(p_hi,xdata(hi));
%     rsq2.hi = rsquare(ydata(hi), yfit_hi); 
    
    if suppressplot == 0
        hold on
        m2 = lm2.Coefficients{2,1};
        b2 = lm2.Coefficients{1,1};
        t(2) = plot([min(xdata(hi)),max(xdata(hi))], ...
            [m2*min(xdata(hi))+b2, m2*max(xdata(hi))+b2],...
            'linewidth',1.5,'color','b', 'linestyle', lstyle);
    else
        t(2) = 0; 
    end
else
    rsq2.hi = NaN;
end

R2Total = (lm1.SSR+lm2.SSR)/(lm1.SST+lm2.SST);

