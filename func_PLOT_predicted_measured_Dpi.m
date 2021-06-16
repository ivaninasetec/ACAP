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

function func_PLOT_predicted_measured(DPIPREDICTED,DPIMEASURED,NAME)

[R2,RMSE,ME] = func_math.Rsquared_rmse(log(DPIMEASURED),log(DPIPREDICTED));

dpipredict_x = linspace(min(DPIMEASURED),max(DPIMEASURED),100);

fitmeas = fitdist(DPIMEASURED(DPIMEASURED>1E-15)', 'Loglogistic');
fitpredic = fitdist(DPIPREDICTED(DPIPREDICTED>1E-15)', 'Loglogistic');

mindpi = min(icdf(fitmeas,0.01),icdf(fitpredic,0.01));
maxdpi = max(icdf(fitmeas,0.99),icdf(fitpredic,0.99));

%% START FIGURE
figure1 = figure(...
    'Name',NAME,...
    'Color',[1,1,1],...
    'Units','centimeters',... % centimeters, pixels, points, characters
    'Position',[0,0,10,10],...%[left bottom width height]
    'PaperSize',[10,10]);
axes1 = axes(...
    'Parent',figure1,...
    'Units','normalized',... %normalized (max 1,1), inches centimeters points pixels
    'OuterPosition',[0.01,0,1.03,1.0]);%Accounting for text, titles...[left bottom width height]
%     'Position',[0.1,0.1,0.9,0.9]...left bottom width height...

hold(axes1,'on');
s1=3;
c1 = [0,0,0];
scatter(DPIMEASURED,DPIPREDICTED,s1,c1,...
    'MarkerFaceAlpha',0.5,...
    'MarkerEdgeAlpha',0.5,...
    'MarkerEdgeColor',[0 0 0]);

plot([min(DPIMEASURED),max(DPIMEASURED)],[min(DPIMEASURED),max(DPIMEASURED)],'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

title({NAME,'Pore size: measured vs predicted'});
xlabel('measured (m)');
ylabel('predicted (m)');
xlim([mindpi,maxdpi]);
ylim([mindpi,maxdpi]);
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
    'YColor',[0 0 0],...
    'XScale','log',...
    'YScale','log',...
    'XMinorTick','off',...
    'YMinorTick','off');%,...    
%     'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
%     'XTickLabel',{'0%','','','','','50%','','','','','100%'},...

%     'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
%     'YTickLabel',{'0%','','','','','50%','','','','','100%'});

% Create textbox
annotation(figure1,'textbox',...
    [0.55 0.15 0.35 0.1],... %[x y w h]
    'String',{strcat('R2(log)   = ',num2str(round(R2,4))),strcat('RMSE(log) = ',num2str(round(RMSE,4))),strcat('ME(log)   = ',num2str(round(ME,4)))},...
    'LineStyle','none',...
    'Color',[0 0 0],...
    'FontWeight','bold',...
    'FontSize',8,...
    'FontName','Courier',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);

hold off

end
