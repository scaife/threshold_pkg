function results = thresholdAnalysisPlotGenerator(...
    Y, X, tLabel)

f = figure; 
monteCarloThresh = 10000;
monteCarloSims = 1e3;
yLimMax = ceil(max(Y)/10)*10;

idx = all([~isnan(X),~isnan(Y)],2);
x = X(idx); 
y = Y(idx);
[results.thresh, results.r2_total] = stats_regress2lines_m2(x,y,3,0,...
    monteCarloThresh,...
    monteCarloSims);
[results.r2,results.pval,results.resids,results.slope,results.intercept, results.r2_total2] ...
    = thresholdplotter(results.thresh(1),x,y,'-',0);

ylim([0 yLimMax])

figsettings(f,'','',tLabel, 14)

drawnow