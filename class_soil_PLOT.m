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

classdef class_soil_PLOT
    %soil_PLOT Class to plot soil objects (in any of the models)
    %   Example:
    %     soilAP81 = soil_AP81(...
    %     'GSD',F0.GSD,...
    %     'WRC',F0.WRC,...
    %     'npor',F0.npor,...
    %     'alfa',1.38);
    %
    % soilplot = soil_PLOT(soilAP81)
    
    properties
        soil
        ppmin_from_GSD
        ppmax_from_GSD
        ppmax_from_WRC
        ppmin_from_WRC
        
        Dpmin_from_GSD
        Dpmax_from_GSD
        Dpmin_from_WRC
        Dpmax_from_WRC
        Dp_range_WRC
        Dp_range_GSD
        
        pp_range_WRC
        pp_range_GSD
        
        Sei_from_WRC
        Sei_from_GSD
        
        thi_from_WRC
        
        Dpi_from_WRC
        ppi_from_WRC
        
        Plot
        Axes
    end
    
    methods
        function obj = class_soil_PLOT(soil)
            %class_soil_PLOT Constructor of an instance of the class_soil_PLOT class.
            
            %soil.WRC = soil.WRC(soil.WRC(:,1)>0,:);
            soil.WRC(soil.WRC(:,1)<1E-6,1) = 1e-6;
            %soil.GSD=soil.GSD(soil.GSD(:,1)>0,:);
            soil.GSD(soil.GSD(:,1)<1E-6,1) = 1e-6;
            
            obj.soil = soil;
            obj.Sei_from_WRC = soil.Se_pp_from_WRC(soil.WRC(:,1));
            obj.thi_from_WRC = soil.WRC(:,2);
            obj.Sei_from_GSD = soil.GSD(:,2);
            
            obj.ppmin_from_WRC = max(1E-6,min(soil.WRC(:,1)));
            obj.ppmax_from_WRC = max(soil.WRC(:,1));
            
            obj.Dpi_from_WRC = soil.Dp_pp_from_WRC(soil.WRC(:,1));
            obj.ppi_from_WRC = soil.WRC(:,1);
            
            obj.Dpmin_from_WRC = min(obj.Dpi_from_WRC);
            obj.Dpmax_from_WRC = max(obj.Dpi_from_WRC);
            
            DPI_from_GSD = soil.Dp_Dg_from_GSD(soil.GSD(:,1));
            DPI_from_GSD = DPI_from_GSD(DPI_from_GSD>0);
            pp_from_GSD  = soil.pp_Dp_from_WRC(DPI_from_GSD);
            
            obj.Dpmin_from_GSD = min(DPI_from_GSD);
            obj.Dpmax_from_GSD = max(DPI_from_GSD);
            obj.ppmin_from_GSD = min(pp_from_GSD);
            obj.ppmax_from_GSD = max(pp_from_GSD);
            
            obj.Dp_range_WRC = logspace(log10(obj.Dpmin_from_WRC),log10(obj.Dpmax_from_WRC),100);
            obj.Dp_range_GSD = logspace(log10(obj.Dpmin_from_GSD),log10(obj.Dpmax_from_GSD),100);
            
            obj.pp_range_WRC = logspace(log10(obj.ppmin_from_WRC),log10(obj.ppmax_from_WRC),100);
            obj.pp_range_GSD = logspace(log10(obj.ppmin_from_GSD),log10(obj.ppmax_from_GSD),100);
            
        end
        
        function obj = show_PSD(obj,NAME)
            %show_PSD With this method a chart of the pore size
            %distribution is configured.
            
            obj.Plot = figure(...
                'Name',NAME,...
                'Color',[1,1,1],...
                'Units','centimeters',... % centimeters, pixels, points, characters
                'Position',[0,0,14,8],...%[left bottom width height]
                'PaperSize',[14,8]); %[w h] For printing purposes
            
            obj.Plot.CurrentAxes = axes(...
                'Parent',obj.Plot ,...
                'Units','normalized',... %normalized (max 1,1), inches centimeters points pixels
                'OuterPosition',[0.01,0,1.03,1.0]);%Accounting for text, titles...[left bottom width height]
            
            box(obj.Plot.CurrentAxes,'on');
            %grid(obj.Plot.CurrentAxes,'on');
            
            set(obj.Plot.CurrentAxes,...
                'FontSize',7,...
                'FontName','Arial',...
                'GridAlpha',1,...
                'GridColor',[0 0 0],...
                'GridLineStyle',':',...
                'XColor',[0 0 0],...
                'xscale','log',...
                'YColor',[0 0 0],...
                'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
                'YTickLabel',{'0%','','','','','50%','','','','','100%'});
            ylim([0,1]);
            xlim([0.01*obj.Dpmin_from_WRC,1e-2]);
            hold(obj.Plot.CurrentAxes,'on');
            s1=25;
            c1 = [0,0,0];
            
            scatter(obj.Plot.CurrentAxes,obj.Dpi_from_WRC,obj.Sei_from_WRC,s1,c1,...
                'MarkerFaceAlpha',1.0,...
                'MarkerEdgeAlpha',1.0,...
                'MarkerFaceColor',[0 0 0],...
                'MarkerEdgeColor',[0 0 0]);
            
            plot(obj.Plot.CurrentAxes,obj.Dp_range_WRC,obj.soil.FVp_Dp_from_WRC(obj.Dp_range_WRC),'LineStyle','-','LineWidth',1,'Color',[0 0 0]);
            
            %title(NAME);
            xlabel('pore size (m)');
            ylabel('Acumulated pore volume');
            
            %             legend('PSD-measured (tests)','PSD-measured (curve)','PSD-predicted (curve)','Location','southeast');
        end
        
        function obj = add_PSD_from_GSD(obj,soil,linstyle,linwith,color)
            %This method adds aditional Pore Size distribution from other
            %models (aditional model are input as soil oject argument).
            plot(obj.Plot.CurrentAxes, obj.Dp_range_GSD,soil.FVp_Dp_from_GSD(obj.Dp_range_GSD),'LineStyle',linstyle,'LineWidth',linwith,'Color',color);
        end
        
        function obj = show_WRC(obj,NAME)
            %show_PSD With this method a chart of the pore size
            %distribution is configured.
            
            obj.Plot = figure(...
                'Name',NAME,...
                'Color',[1,1,1],...
                'Units','centimeters',... % centimeters, pixels, points, characters
                'Position',[0,0,14,8],...%[left bottom width height]
                'PaperSize',[14,8]); %[w h] For printing purposes
            
            obj.Plot.CurrentAxes = axes(...
                'Parent',obj.Plot ,...
                'Units','normalized',... %normalized (max 1,1), inches centimeters points pixels
                'OuterPosition',[0.01,0,1.03,1.0]);%Accounting for text, titles...[left bottom width height]
            
            box(obj.Plot.CurrentAxes,'on');
            grid(obj.Plot.CurrentAxes,'on');
            
            set(obj.Plot.CurrentAxes,...
                'FontSize',7,...
                'FontName','Arial',...
                'GridAlpha',1,...
                'GridColor',[0 0 0],...
                'GridLineStyle',':',...
                'XColor',[0 0 0],...
                'xscale','log',...
                'YColor',[0 0 0]);%,...
%                 'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
%                 'YTickLabel',{'0%','','','','','50%','','','','','100%'});
            ylim([0,1.5*max(obj.thi_from_WRC)]);
            xlim([1E-2,1e+5]);
            hold(obj.Plot.CurrentAxes,'on');
            s1=10;
            c1 = [0,0,0];
            
            scatter(obj.Plot.CurrentAxes,obj.ppi_from_WRC,obj.thi_from_WRC,s1,c1,...
                'MarkerFaceAlpha',1.0,...
                'MarkerEdgeAlpha',1.0,...
                'MarkerEdgeColor',[0 0 0]);
            
            plot(obj.Plot.CurrentAxes,obj.pp_range_WRC,obj.soil.th_pp_from_WRC(obj.pp_range_WRC),'LineStyle','-','LineWidth',1,'Color',[0 0 0]);
            
            %title(NAME);
            xlabel('pore size (m)');
            ylabel('Acumulated pore volume');
            
            %             legend('PSD-measured (tests)','PSD-measured (curve)','PSD-predicted (curve)','Location','southeast');
        end
        
        function obj = add_WRC_from_GSD(obj,soil,linstyle,linwith,color)
            %This method adds aditional Pore Size distribution from other
            %models (aditional model are input as soil oject argument).
            plot(obj.Plot.CurrentAxes, obj.pp_range_GSD,soil.th_pp_from_GSD(obj.pp_range_GSD),'LineStyle',linstyle,'LineWidth',linwith,'Color',color);
        end
        
    end
end

