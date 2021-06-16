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

function func_PLOT_predicted_measured_Sei(SEIPREDICTED,SEIMEASURED,NAME)
RESIDUAL = (SEIPREDICTED-SEIMEASURED);
RESIDUALABS = abs(RESIDUAL);
% RESIDUALPOS = RESIDUAL(RESIDUAL>0);
% SEIPOS = SEIMEASURED(RESIDUAL>0);
% NLENGTHPOS = length(RESIDUALPOS);
% NLENGTH = length(RESIDUAL);
% RESIDUALNEG = RESIDUAL(RESIDUAL<0);
% SEINEG = SEIMEASURED(RESIDUAL<0);

% NLENGTHNEG = length(RESIDUALNEG);
% RESIDUALPOSREL = RESIDUALPOS./(1-SEIPOS);
% RESIDUALNEGREL = -RESIDUALNEG./SEINEG;

[R2,RMSE,ME] = func_math.Rsquared_rmse(SEIMEASURED,SEIPREDICTED);

% fitpos = fitdist(RESIDUALPOSREL', 'half normal');
% fitneg = fitdist(RESIDUALNEGREL', 'half normal');
fitall = fitdist(RESIDUALABS', 'half normal');

% bound95POS = icdf(fitpos,0.95,0);
% bound95NEG = icdf(fitneg,0.95,0);
bound95ALL = icdf(fitall,0.90,0);

seipredict_x = linspace(min(SEIMEASURED),max(SEIMEASURED),100);
% seipredict_y_sup = seipredict_x+bound95POS.*(1-seipredict_x);
% seipredict_y_inf = seipredict_x-bound95NEG.*seipredict_x;
seipredict_y_sup = seipredict_x+bound95ALL;
seipredict_y_inf = seipredict_x-bound95ALL;

%% START FIGURE
figure1 = figure(...
    'Name',NAME,...
    'Color',[1,1,1],...
    'Units','centimeters',... % centimeters, pixels, points, characters
    'Position',[0,0,10,10],...%[left bottom width height]
    'PaperSize',[10,10]); %[w h] For printing purposes
axes1 = axes(...
    'Parent',figure1,...
    'Units','normalized',... %normalized (max 1,1), inches centimeters points pixels
    'OuterPosition',[0.01,0,1.03,1.0]);%Accounting for text, titles...[left bottom width height]
%     'Position',[0.1,0.1,0.9,0.9]...left bottom width height...

hold(axes1,'on');
s1=5;
c1 = [0,0,0];
scatter(SEIMEASURED,SEIPREDICTED,s1,c1,...
    'MarkerFaceAlpha',0.5,...
    'MarkerEdgeAlpha',0.5,...
    'MarkerEdgeColor',[0 0 0]);

% plot(seipredict_x,seipredict_y_sup,'LineStyle','--','LineWidth',2,'Color',[0 0 0]);
% plot(seipredict_x,seipredict_y_inf,'LineStyle','--','LineWidth',2,'Color',[0 0 0]);
plot([0,1],[0,1],'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

title({NAME,'Saturation: measured vs predicted'});
xlabel('measured');
ylabel('predicted');
xlim([0,1]);
ylim([0,1]);
box(axes1,'on');
grid(axes1,'on');

set(axes1,...
    'DataAspectRatio',[1 1 1],...
    'FontSize',8,...
    'GridAlpha',1,...
    'GridColor',[0 0 0],...
    'GridLineStyle',':',...
    'PlotBoxAspectRatio',[1 1 1],...
    'XColor',[0 0 0],...
    'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
    'XTickLabel',{'0%','','','','','50%','','','','','100%'},...
    'YColor',[0 0 0],...
    'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
    'YTickLabel',{'0%','','','','','50%','','','','','100%'});

% Create textbox
annotation(figure1,'textbox',...
    [0.65 0.15 0.25 0.1],... %[x y w h]
    'String',{strcat('R2   = ',num2str(round(R2,4))),strcat('RMSE = ',num2str(round(RMSE,4))),strcat('ME   = ',num2str(round(ME,4)))},...
    'LineStyle','none',...
    'Color',[0 0 0],...
    'FontWeight','bold',...
    'FontSize',8,...
    'FontName','Courier',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);

hold off

end

