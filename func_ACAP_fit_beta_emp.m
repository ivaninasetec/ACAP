%-------------------------------------------------------------------------------
%MIT License
% 
% Copyright (c) 2021 Ivan Campos-Guereta Diez
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% - The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
% - Campos-Guereta et al, 2020 paper must be quoted whith the following data when 
%   refering to this code or when the code is applied: 
%   'CAMPOS-GUERETA, I., DAWSON, A. & THOM, N. 2021. An alternative continuous form
%    of Arya and Paris model to predict the soil water retention curve of a soil. 
%    Advances in Water Resources, 154, 103968.'
%    https://doi.org/10.1016/j.advwatres.2021.103968
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%-------------------------------------------------------------------------------

function func_ACAP_fit_beta_emp(F0_70, StrEquation, optbeta)
F0_70.LOGD10REL = log(F0_70.D10./F0_70.D50);
F0_70.LOGD20REL = log(F0_70.D20./F0_70.D50);
F0_70.LOGD30REL = log(F0_70.D30./F0_70.D50);
F0_70.LOGD40REL = log(F0_70.D40./F0_70.D50);
F0_70.LOGD60REL = log(F0_70.D60./F0_70.D50);
F0_70.LOGD80REL = log(F0_70.D80./F0_70.D50);

F0_70.CLAYB = F0_70.CLAY./(1-F0_70.GRAVEL);
F0_70.SILTB = F0_70.SILT./(1-F0_70.GRAVEL);
F0_70.SANDB = F0_70.SAND./(1-F0_70.GRAVEL);

fitF0 = fitlm(F0_70,StrEquation,'RobustOpts','off');

scatter(fitF0.Variables.BETA_FIX_FIT,fitF0.Fitted,'.k');
hold on
title(strcat('\beta_{fix}: measured vs fitted (R2=',num2str(fitF0.Rsquared.Ordinary),' RMSE:',num2str(fitF0.RMSE),')'));
minbeta = min([fitF0.Variables.BETA_FIX_FIT',fitF0.Fitted']);
maxbeta = max([fitF0.Variables.BETA_FIX_FIT',fitF0.Fitted']);
xlim([minbeta, maxbeta]);
ylim([minbeta, maxbeta]);
ylabel('\beta_{fitted}');
xlabel('\beta_{measured}');
plot([minbeta,maxbeta],[minbeta,maxbeta],'k','LineWidth',2);
grid on
hold off

EquationStr = strcat(optbeta,' =', num2str(fitF0.Coefficients.Estimate(1),'%+e'));
for i=1:length(fitF0.Coefficients.Estimate)-1
    EquationStr = strcat(EquationStr,num2str(fitF0.Coefficients.Estimate(1+i),'%+f'),'.*',fitF0.Coefficients.Properties.RowNames{1+i});
end
    disp(EquationStr);
end
