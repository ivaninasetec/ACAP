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

classdef class_soil_CNEAP
    %class_soil_CNEAP Class to build the CNEAP model (Nimmo et al, 2007).
    
    properties
        param_Ts       = 0.07197;  % Surface tension of water in 71.97 mN/m at 25�C
        param_rhow     = 997.0479; % Density of water in kg/m3 at 25�C
        param_g        = 9.81;         % Acceleration of gravity 9.81m/s2
        param_costheta = 1.00;  % Contact angle water-air-soil
        
        param_c = 0.233;
        param_d = 1.109;
        
        GSD=[]     %Grain Size Distribution (Size[m],Cum_vol[per unit])
        WRC=[]     %Water Retention Curve (optional) (Suction [m],Volumetric Water Content[m3/m3])
        npor=nan   %porosity index
        evoid = nan %void ratio
        alfa_fix = 1.38; %Best estimate for alfa in the Aria and Paris, 1981 model.
        di = [2E-6,53E-6,106E-6,0.180E-3,0.250E-3,0.355E-3,0.5E-3,0.71E-3,1E-3,2E-3] %Fixed diameters for the discretization
        
        dimean = nan % Mean diameter in each fraction
        dpimean = nan % Mean diameter of the pores
        incni = nan  % Number of particles in each fraction
        Wg0 = 0.001 % Unit weight in kg
        Vgi = nan     % Corresponding volumes on the discretization
        incVgi = nan
        rhos = 2650 %Density of particles in kg/m3
        outtable    %Table with some output results
        dmax = nan  %max diameter
        texture = nan %texture of the soil
        alfa = 1.38 %Parameter of the Arya and Paris, 1981 model.
        thres = 0 %Residual volumetric water content
        thsat = nan %Saturated volumetric water content
        
        Dg_max = nan  %Particle size at FVg=100%
        Dg_min = nan  %Particle size at FVg=0%
        Dp_max = nan  %Max pore size (linked to Dg_max)
        Dp_min = nan  %Min pore size (linked to Dg_min)
    end
    
    methods
        %% CONSTRUCTOR
        function obj = class_soil_CNEAP(varargin)
            %class_soil_CNEAP Construct an instance of this class
            def_GSD     = obj.GSD;
            def_npor    = obj.npor;
            def_alfa    = obj.alfa_fix;
            def_WRC     = obj.WRC;
            
            p=inputParser;
            addOptional(p,'GSD',def_GSD);
            addOptional(p,'npor',def_npor);
            addOptional(p,'alfa',def_alfa);
            addOptional(p,'WRC',def_WRC);
            
            parse(p,varargin{:});
            
            obj.npor    = p.Results.npor;
            obj.evoid   = obj.npor/(1-obj.npor);
            obj.thsat = obj.npor;
            
            if (not(isempty(p.Results.GSD)))
                obj.GSD = p.Results.GSD(p.Results.GSD(:,1)>0,:);
                obj.GSD = sortrows([max(p.Results.GSD(:,1)),1;obj.GSD;0,0]);%Add DGmax with 100% in volume, and 0% in Dg=0
                obj.GSD = obj.GSD(~isnan(obj.GSD(:,1)),:);
                obj.GSD = obj.GSD(~isnan(obj.GSD(:,2)),:);
                [~,ia,~] = unique(obj.GSD(:,1),'rows');
                obj.GSD = obj.GSD(ia,:); %->GSD
                obj.dmax = max(obj.GSD(:,1));
            end
            
            obj.texture = obj.get_texture();
            
            
            %Include measured WRC if exist
            if(not(isempty(p.Results.WRC))) %->WRC
                obj.WRC = p.Results.WRC;
                [~,ia,~] = unique(obj.WRC(:,1),'rows');
                obj.WRC = obj.WRC(ia,:);
            end
            
            obj.Dg_max  = obj.Dg_FVg_from_GSD(1);
            obj.Dg_min  = obj.Dg_FVg_from_GSD(0);
            obj.Dp_max  = obj.Dp_Dg_from_GSD(obj.Dg_max);
            obj.Dp_min  = obj.Dp_Dg_from_GSD(obj.Dg_min);
            
            obj.outtable = obj.get_tableoutput();
        end
        
        %% INCLUDE RESULTS PREDICTIONS AND WRC POINTS INTO A TABLE
        function tableout = get_tableoutput(obj)
            %Output a table with the results of the model in the points of
            %the Water Retention Curve.
            if(isempty(obj.WRC)||isempty(obj.GSD)||length(obj.WRC)<=3||length(obj.GSD)<=3)
                dg_i=nan;
                Se_i=nan;
                pp_i=nan;
                Dp_i=nan;
                Se_AP = nan;
                
            else
                Se_i = max(0,min(1,(obj.WRC(:,2)'-obj.thres)./(obj.thsat-obj.thres)));
                th_i = obj.WRC(:,2)';
                mask1 = Se_i>=0;
                mask2 = Se_i<=1;
                mask3 = obj.WRC(:,1)'>0;
                mask = (mask1.*mask2.*mask3)==1;
                Se_i = Se_i(mask);
                th_i = th_i(mask);
                if (isempty(Se_i))
                    dg_i=nan;
                    Se_i=nan;
                    pp_i=nan;
                    Dp_i=nan;
                    Se_CNEAP = nan;
                    th_CNEAP = nan;
                else
                    pp_i = obj.WRC(mask,1)';
                    dg_i = obj.Dg_FVg_from_GSD(Se_i);
                    Dp_i = obj.Dp_pp_from_WRC(pp_i);
                    
                    Se_CNEAP  = obj.FVp_Dp_from_GSD(Dp_i);
                    dpi_CNEAP = obj.Dp_Dg_from_GSD(dg_i);
                    th_CNEAP  = obj.th_pp_from_GSD(pp_i);
                end
            end
            tableout = array2table([dg_i',pp_i',Dp_i',Se_i',th_i' ,Se_CNEAP' ,dpi_CNEAP' ,th_CNEAP'],...
                'VariableNames',   {'DGI','PPI','DPI','SEI','THI','SEICNEAP','DPI_CNEAP','THI_CNEAP'});
            tableout.TEXTURE(:) = {obj.texture};
            
            
        end
        
        
        %% SOIL METHODS
        
        function output = get_texture(obj)
            %Calculate soil texture and clay, silt, sand and gravel
            %content.
            GRAVEL = max(0.0,1.0-obj.FVg_Dg_from_GSD(0.002));
            CLAY = obj.FVg_Dg_from_GSD(0.000002)/(1.0-GRAVEL);
            SILT = (obj.FVg_Dg_from_GSD(0.00005)-CLAY)/(1.0-GRAVEL);
            SAND = (1.0-SILT-CLAY)/(1.0-GRAVEL);
            
            if ((SILT + 1.5*CLAY) < 0.15)
                output='sand';
            elseif ((SILT + 1.5*CLAY >= 0.15) && (SILT + 2*CLAY < 0.30))
                output='loamy sand';
            elseif ((CLAY >= 0.07 && CLAY < 0.20) && (SAND > 0.52) && ((SILT + 2*CLAY) >= 0.30) || (CLAY < 0.07 && SILT < 0.50 && (SILT+2*CLAY)>=0.30))
                output='sandy loam';
            elseif ((CLAY >= 0.07 && CLAY < 0.27) && (SILT >= 0.28 && SILT < 0.50) && (SAND <= 0.52))
                output='loam';
            elseif ((SILT >= 0.50 && (CLAY >= 0.12 && CLAY < 0.27)) || ((SILT >= 0.50 && SILT < 0.80) && CLAY < 0.12))
                output='silt loam';
            elseif (SILT >= 0.80 && CLAY < 0.12)
                output='silt';
            elseif ((CLAY >= 0.20 && CLAY < 0.35) && (SILT < 0.28) && (SAND > 0.45))
                output='sandy clay loam';
            elseif ((CLAY >= 0.27 && CLAY < 0.40) && (SAND > 0.20 && SAND <= 0.45))
                output='clay loam';
            elseif ((CLAY >= 0.27 && CLAY < 0.40) && (SAND  <= 0.20))
                output='silt clay loam';
            elseif (CLAY >= 0.35 && SAND > 0.45)
                output='sandy clay';
            elseif (CLAY >= 0.40 && SILT >= 0.40)
                output='silty clay';
            elseif (CLAY >= 0.40 && SAND <= 0.45 && SILT < 0.40)
                output='clay';
            else
                output='';
            end
        end
        
        
        %% PREDICTIONS FROM GSD:
        %Grain Size Distribution (GSD)-------------------------------------
        function output = FVg_Dg_from_GSD(obj,Dg)
            %Return the agregated grain volume (per unit) given the grain
            %size, from the GSD.
            fittolog = true;
            if(size(obj.GSD,1)>=2)
                isinv=false;
                lder=false;
                output = func_soil.Interpolate_GSD(obj.GSD,Dg,lder,fittolog,isinv);
            else
                output = nan;
            end
            output = max(0,min(1,output));
        end
        
        function output = fVg_Dg_from_GSD(obj,Dg)
            %density probability function in volume for Dg
            %(mean volume of particles per total volume) (Derivative of the
            %GSD).
            fittolog = 1;
            if(size(obj.GSD,1)>=2)
                isinv=false;
                lder=true;
                output = func_soil.Interpolate_GSD(obj.GSD,Dg,lder,fittolog,isinv);
            else
                output = NaN;
            end
        end
        
        function output = Dg_FVg_from_GSD(obj,FVg)
            %Return the grain size given the aggregated particle volume
            %from the GSD.
            fittolog = true;
            if(size(obj.GSD,1)>=2)
                isinv=true;
                lder=false;
                output = func_soil.Interpolate_GSD(obj.GSD,FVg,lder,fittolog,isinv);
            else
                output = nan;
            end
        end
        
        
        %Pore Size Distribution (PSD)--------------------------------------
        
        function out_Dp = Dp_FVp_from_GSD(obj,FVp)
            %Predicted Pore size given the Cummulated Pore Volume (predicted from the GSD and porosity).
            Dg = obj.Dg_FVg_from_GSD(min(1.0,max(0.0,(FVp))));
            out_Dp = obj.Dp_Dg_from_GSD(Dg);
        end
        
        function output = FVp_Dp_from_GSD(obj,Dp)
            %Predicted Cummulated Pore Volume for a given pore diameter.
            %(Aggregated Pore Size Distribution) (predicted from the GSD and porosity)
            Dg = obj.Dg_Dp_from_GSD(Dp);
            output = min(1,max(0,obj.FVg_Dg_from_GSD(Dg)));
        end
        
        function output = Dp_Dg_from_GSD(obj,Dg)
            %Return the pore size linked to the grain size (in the ACAP
            %model for the given beta value.
            Dg=max(0.0,Dg);
            output = max(obj.Dp_min,min(obj.Dp_max,Dg.*(max(1E-10,3./(2.*obj.evoid) .*(2.*obj.Wg0./(pi.*obj.rhos) .* obj.fVg_Dg_from_GSD(Dg)./(Dg.^2))).^(obj.alfa-1) ).^(-1./2)));
            %             output = interp1(obj.dgimean',obj.dpimean',Dg','linear','extrap');
        end
        
        function out_Dg = Dg_Dp_from_GSD(obj,Dp)
            %Predicted grain size given the pore size (predicted from the GSD data and from the porosity)
            out_Dg = zeros(size(Dp));
            
            options = optimset;
            options = optimset(options,'Display', 'off');
            
            for i=1:length(Dp)
                if(Dp(i)>=obj.Dp_max)
                    out_Dg(i) = obj.Dg_max;
                elseif (Dp(i)<=obj.Dp_min)
                    out_Dg(i) = obj.Dg_min;
                else
                    out_Dg(i) = max(obj.Dg_min,min(obj.Dg_max,fzero(@(x)(obj.Dp_Dg_from_GSD(x)-Dp(i)),exp(log(obj.Dg_max+obj.Dg_min)./2),options)));
                end
            end
            
        end
        
        
        %Water Retention Curve (WRC)---------------------------------------
        
        function output = Se_pp_from_GSD(obj,pp)
            %Predicted relative saturation for a given suction (predicted from the GSD and porosity).
            Dp_temp = obj.Dp_pp_from_WRC(pp);
            Dg = obj.Dg_Dp_from_GSD(Dp_temp);
            output = obj.FVg_Dg_from_GSD(Dg);
        end
        
        function out_th = th_pp_from_GSD(obj,pp)
            %Predicted volumetric water content given the suction (predicted from the GSD and porosity).
            out_th = obj.thres+(obj.thsat-obj.thres).*obj.Se_pp_from_GSD(pp);
        end
        
        function out_pp = pp_Se_from_GSD(obj,Se)
            %Predicted suction given the relative saturation (predicted from the GSD and porosity).
            Dg = obj.Dg_FVg_from_GSD(min(1.0,max(0.0,Se)));
            out_pp = obj.pp_Dg_from_GSD(Dg);
        end
        
        function out_pp = pp_th_from_GSD(obj,th)
            %Predicted suction given the volumetric water content (predicted from the GSD and porosity).
            Se = min(1.0,max(0.0,(th-obj.thres)./(obj.thsat-obj.thres)));
            Dg = obj.Dg_FVg_from_GSD(Se);
            out_pp = obj.pp_Dg_from_GSD(Dg);
        end
        
        function out_pp = pp_Dg_from_GSD(obj,Dg)
            %Predicted suction given the linked grain size (predicted from the GSD and porosity).
            Dp = obj.Dp_Dg_from_GSD(Dg);
            out_pp = obj.pp_Dp_from_WRC(Dp);
        end
        
        
        %% PREDICTIONS BASED ON THE MEASURED WATER RETENTION CURVE
        
        %Pore Size Distribution (PSD)--------------------------------------
        function out_Dp = Dp_FVp_from_WRC(obj,FVp)
            %Pore size given the cumulated pore volume, calculated from
            th = obj.thres+(obj.thsat-obj.thres).*min(0.0,max(1.0,FVp));
            pp = obj.pp_th_from_WRC(th);
            out_Dp = obj.Dp_pp_from_WRC(pp);
        end
        
        function output = FVp_Dp_from_WRC(obj,Dp)
            %Cumulated pore volume givven the pore size, calculated from
            %the WRC.
            pp = obj.pp_Dp_from_WRC(Dp);
            output = min(1,max(0,obj.Se_pp_from_WRC(pp)));
        end
        
        function out_pp = pp_Dp_from_WRC(obj,Dp)
            %Suction given the pore diameter (Youngs Equation).
            out_pp = (4.*obj.param_Ts.*obj.param_costheta)./(obj.param_rhow.*obj.param_g.*Dp);
        end
        
        function out_Dp = Dp_pp_from_WRC(obj,pp)
            %Pore size given the suction (Youngs equation).
            out_Dp = (4.*obj.param_Ts.*obj.param_costheta)./(obj.param_rhow.*obj.param_g.*pp);
        end
        
        %Water Retention Curve (WRC)---------------------------------------
        function out_pp = pp_th_from_WRC(obj,th)
            %Suction given the volumetric water content from the WRC.
            derivative = false;
            isinv=true;
            thtemp = max(obj.thres,min(obj.thsat,th));
            out_pp = func_soil.Interpolate_WRC(obj.WRC,thtemp,derivative,obj.thres,obj.thsat,obj.options.fit_to_log,isinv);
        end
        
        function output = th_pp_from_WRC(obj,pp)
            %suction given the volumetric water content from the WRC.
            thres_local = 0; %This model considers thres=0
            fittolog = true;
            derivative = false;
            isinv=false;
            output = func_soil.Interpolate_WRC(obj.WRC,pp,derivative,thres_local,obj.thsat,fittolog,isinv);
        end
        
        
        function out_pp = pp_Se_from_WRC(obj,Se)
            %Suction given the relative saturation from the WRC.
            th = obj.thres+(obj.thsat-obj.thres).*max(0.0,min(1.0,Se));
            out_pp = obj.pp_th_from_WRC(th);
        end
        
        function out_Se = Se_pp_from_WRC(obj,pp)
            %Relative saturation given the suction from the WRC.
            thres_local = 0; %This model considers thres=0
            th_temp = obj.th_pp_from_WRC(pp);
            out_Se = min(1,max(0,(th_temp-thres_local)./(obj.thsat-thres_local)));
        end
        
    end
    
    
end




