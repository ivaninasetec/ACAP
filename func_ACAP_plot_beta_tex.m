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

function func_ACAP_plot_beta_tex(DATA_FIT)
boxplot(DATA_FIT.BETA_FIX_FIT,DATA_FIT.TEXTURE2,...
    'Colors',[0,0,0],...
    'MedianStyle','line',...%target or line
    'Widths',0.3,...
    'Notch','off',...
    'Symbol','k',...
    'Jitter',0);%symbol and color of outliers
fig1= gcf;

fig1.PaperUnits='centimeters';
fig1.PaperPosition=[0,0,8.5,7.5];
fig1.PaperSize=[8.5,7.5];

axe = get(fig1,'CurrentAxes');
% get(gcf,'CurrentAxes');
axe.FontName = 'Arial';
axe.FontWeight = 'bold';
axe.FontSizeMode = 'manual';
axe.FontSize = 7;

axe.XGrid = 'off';
axe.YGrid = 'on';
axe.YMinorGrid = 'on';
axe.GridAlpha = 1;
axe.MinorGridAlpha = 1;
axe.GridColor = [0.5,0.5,0.5];
axe.MinorGridColor = [0.65,0.65,0.65];
axe.YTick = [1,1.5,2,2.5];

xtickangle(axe,45);

fig1.Color=[1,1,1];
fig1.Units='centimeters';
fig1.Clipping='off';
fig1.Resize='off';
fig1.Position=[5,10,8,7];


ylabel('\beta_{fix.tex}');
end
