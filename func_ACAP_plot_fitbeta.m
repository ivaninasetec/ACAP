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

function func_ACAP_plot_fitbeta(DATA_FIT_I,COD,NAME)
F1 = DATA_FIT_I(DATA_FIT_I.COD==COD,:);
%DG0=F1.D50./0.001;
DG0=F1.DG0;
numerator = sqrt(3)./(sqrt(2).*F1.EVOID.^(3/2))  .*  pi./4  .* (F1.DPI./DG0).^3;
denominator = pi./6  .* (F1.DGI./DG0).^3;

lognumerator = log(numerator);
logdenominator = log(denominator);

lmfit   = fitlm(logdenominator,lognumerator,'Intercept',false);
outbeta = lmfit.Coefficients{1,1};
RMSE    = lmfit.RMSE;
R2      = lmfit.Rsquared.Ordinary;

minylog = min([lognumerator',outbeta.*logdenominator']);
maxylog = max([lognumerator',outbeta.*logdenominator']);
miny = exp(minylog);
maxy = exp(maxylog);
minx = exp(minylog./outbeta);
maxx = exp(maxylog./outbeta);

close all;
figure1 = figure(...
    'Name',NAME,...
    'Color',[1,1,1],...
    'Units','centimeters',... % centimeters, pixels, points, characters
    'Position',[0,0,8.5,8.5],...%[left bottom width height]
    'PaperSize',[8.5,8.5]); %[w h] For printing purposes
axes1 = axes(...
    'Parent',figure1,...
    'Units','normalized',... %normalized (max 1,1), inches centimeters points pixels
    'OuterPosition',[0.01,0,1.03,1.0]);%Accounting for text, titles...[left bottom width height]
%     'Position',[0.1,0.1,0.9,0.9]...left bottom width height...

hold(axes1,'on');
s1=10;
c1 = [0,0,0];

scatter(denominator,numerator,s1,c1,...
    'MarkerFaceAlpha',1.0,...
    'MarkerEdgeAlpha',1.0,...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[0 0 0]);
plot([minx,maxx],[miny,maxy],'LineStyle','-','LineWidth',2,'Color',[0 0 0]);

%title({['UNSODA COD: ',int2str(COD)],NAME});

xlabel('$$\frac{\pi }{6}{\left( {\frac{{{D_{p,i}}}}{{{D_{g0}}}}} \right)^3}$$', 'Interpreter', 'Latex')
ylabel('$$\frac{{\sqrt 3 }}{{\sqrt 2 \cdot{e^{3/2}}}}\frac{\pi }{4}{\left( {\frac{{{D_{p,i}}}}{{{D_{g0}}}}} \right)^3}$$', 'Interpreter', 'Latex')



% xlabel('\pi/6·(D_g');
% ylabel('predicted');
xlim([minx,maxx]);
ylim([miny,maxy]);
box(axes1,'on');
%grid(axes1,'on');

set(axes1,...
    'FontSize',7,...
    'XColor',[0 0 0],...
    'XScale','log',...
    'YScale','log',...
    'YColor',[0 0 0]);
%     'GridAlpha',1,...
%     'GridColor',[0 0 0],...
%     'GridLineStyle',':',...
%     'DataAspectRatio',[1 1 1],...
%     'PlotBoxAspectRatio',[1 1 1],...
%     'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
%     'XTickLabel',{'0%','','','','','50%','','','','','100%'},...
% 'YColor',[0 0 0],...
%     'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
%     'YTickLabel',{'0%','','','','','50%','','','','','100%'});

% Create textbox
annotation(figure1,'textbox',...
    [0.55 0.20 0.35 0.1],... %[x y w h]
    'String',{strcat('\beta_{fix.fit}   =',num2str(round(outbeta,4))), strcat('R2(log)   = ',num2str(round(R2,4))),strcat('RMSE(log) = ',num2str(round(RMSE,4)))},...
    'LineStyle','none',...
    'Color',[0 0 0],...
    'FontWeight','bold',...
    'FontSize',7,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);

hold off


end

