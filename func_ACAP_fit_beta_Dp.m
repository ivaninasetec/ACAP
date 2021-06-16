function [eqout] = func_ACAP_fit_beta_Dp(DATA_FIT_I,StrEquation,optbeta)
DATA_FIT_I.LOGD10REL = log(DATA_FIT_I.D10./DATA_FIT_I.D50);
DATA_FIT_I.LOGD20REL = log(DATA_FIT_I.D20./DATA_FIT_I.D50);
DATA_FIT_I.LOGD30REL = log(DATA_FIT_I.D30./DATA_FIT_I.D50);
DATA_FIT_I.LOGD40REL = log(DATA_FIT_I.D40./DATA_FIT_I.D50);
DATA_FIT_I.LOGD60REL = log(DATA_FIT_I.D60./DATA_FIT_I.D50);
DATA_FIT_I.LOGD80REL = log(DATA_FIT_I.D80./DATA_FIT_I.D50);
DATA_FIT_I.LOGDPILOGD50 = DATA_FIT_I.LOGDPI.*DATA_FIT_I.LOGD50;
DATA_FIT_I.LOGDPIBETA_FIX_EMP = DATA_FIT_I.LOGDPI.*DATA_FIT_I.BETA_FIX_EMP;
DATA_FIT_I.LOGD50BETA_FIX_EMP = DATA_FIT_I.LOGD50.*DATA_FIT_I.BETA_FIX_EMP;
DATA_FIT_I.LOGDPIBETA_FIX_FIT = DATA_FIT_I.LOGDPI.*DATA_FIT_I.BETA_FIX_FIT;
DATA_FIT_I.LOGD50BETA_FIX_FIT = DATA_FIT_I.LOGD50.*DATA_FIT_I.BETA_FIX_FIT;

DATA_FIT_I.CLAYB = DATA_FIT_I.CLAY./(1-DATA_FIT_I.GRAVEL);
DATA_FIT_I.SILTB = DATA_FIT_I.SILT./(1-DATA_FIT_I.GRAVEL);
DATA_FIT_I.SANDB = DATA_FIT_I.SAND./(1-DATA_FIT_I.GRAVEL);

fitF0 = fitlm(DATA_FIT_I,StrEquation,'RobustOpts','off');

scatter(fitF0.Variables.BETA_FIX_FIT,fitF0.Fitted,'.k');
hold on
title(strcat('\beta_i (R2=',num2str(fitF0.Rsquared.Ordinary),' RMSE:',num2str(fitF0.RMSE),')'));
minbeta = min([fitF0.Variables.BETA_FIX_FIT',fitF0.Fitted']);
maxbeta = max([fitF0.Variables.BETA_FIX_FIT',fitF0.Fitted']);
xlim([minbeta, maxbeta]);
ylim([minbeta, maxbeta]);
xlabel('\beta_{measured}');
ylabel('\beta_{fitted}');
plot([minbeta,maxbeta],[minbeta,maxbeta],'k','LineWidth',2);
grid on
hold off

    EquationStr = strcat(optbeta,' =', num2str(fitF0.Coefficients.Estimate(1),'%+f'));

for i=1:length(fitF0.Coefficients.Estimate)-1
    EquationStr = strcat(EquationStr,num2str(fitF0.Coefficients.Estimate(1+i),'%+f'),'.*',fitF0.Coefficients.Properties.RowNames{1+i});
end
    eqout= EquationStr;
end

