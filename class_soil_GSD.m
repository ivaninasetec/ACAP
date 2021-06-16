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

classdef class_soil_GSD
    %class_soil Parent class that includes all common properties and
    %methods for a soil given its GSD and npor.
    
    properties
        
        GSD=[]        %Grain Size Distribution (Size[m],Cum_vol[per unit])
        npor=nan      %Porosity index
        evoid = nan   %Void ratio

        dmax = nan    %max diameter
        texture = nan %texture of the soil

        Dg_max = nan  %Particle size at FVg=100%
        Dg_min = nan  %Particle size at FVg=0%
        
        D10 = nan
        D20 = nan
        D30 = nan
        D40 = nan
        D50 = nan
        D60 = nan
        D80 = nan
        cu = nan
        cc = nan
        
        clay = nan
        silt = nan
        sand = nan
        gravel = nan
        
        options = struct('interpolation',struct(...
                            'fit_to_log',true,...
                            'extrapolate_dgmax',false),...
                         'liquid',struct(...
                            'param_Ts',0.07197,...
                            'param_rhow',0.07197,...
                            'param_g',0.07197,...
                            'param_costheta',0.07197)...%Extrapolate max Dg to 100% with log linear interpolation
                         ); %Struct with options of the ACAP model       
    end
    
    methods
        %% CONSTRUCTOR
        function obj = class_soil_GSD(varargin)
            %class_soil_AP81 Construct an instance of this class
            def_GSD     = obj.GSD;
            def_npor    = obj.npor;
            def_Dg_max    = nan;
            
            %Parse constructor arguments
            p=inputParser;
            addOptional(p,'GSD',def_GSD);
            addOptional(p,'Dg_max',def_Dg_max);
            addOptional(p,'npor',def_npor);
            
            parse(p,varargin{:});
            
            % Update properties from arguments
            % Void ratio
            obj.npor    = p.Results.npor;
            obj.evoid   = obj.npor/(1-obj.npor);           

            %Update GSD and Dg_max
            obj = obj.update_GSD(p.Results.GSD,p.Results.Dg_max);        
            
        end
        
        function obj = update_GSD(obj,GSD_INPUT,Dgmax)
            %Check the Grain Size Distribution and set some properties
            %dependin on the GSD.
            
            if (not(isempty(GSD_INPUT)))
                obj.GSD = GSD_INPUT;
                %Delete all NaN
                obj.GSD = obj.GSD(~isnan(obj.GSD(:,1)),:);
                obj.GSD = obj.GSD(~isnan(obj.GSD(:,2)),:);
                
                %Include Dg_max point if defined (Dg_max,1) and refresh
                if (not(isempty(Dgmax)||isnan(Dgmax)))   %->Dg_max
                    obj.GSD = [obj.GSD;Dgmax,1.0];
                end
                
                %All values of FVg over 1 are set to 1:
                obj.GSD(obj.GSD(:,2)>=1.0,2)=1.0; %All Dg values lower or equal to 0 are set to 0.
                obj.GSD(obj.GSD(:,2)==1.0,1)= min(obj.GSD(obj.GSD(:,2)==1.0,1)); %All FVg values with FVg=1 set to min Dg(FVg=1).
                
                %All values of FVg below 0 are set to 0:
                obj.GSD(obj.GSD(:,2)<=0.0,2)=0.0; %All Dg values lower or equal to 0 are set to 0
                obj.GSD(obj.GSD(:,1)<=0.0,1)=0.0; %All Dg values lower or equal to 0.0 are set to 0.0
                obj.GSD(obj.GSD(:,2)==0.0,1)= max(obj.GSD(obj.GSD(:,2)==0.0,1)); %All FVg values with FVg=0 set to max Dg(FVg=0).
                
                %Update now Dg_min and Dg_max:
                obj.Dg_max = max(obj.GSD(:,1));
                obj.Dg_min = min(obj.GSD(:,1));
                
                %Extend GSD to (Dg_max,1)
                if (obj.options.interpolation.extrapolate_dgmax)
                logDg = log(obj.GSD(obj.GSD(:,1)>0.0,1));
                %logFVDg = log(obj.GSD(obj.GSD(:,1)>0.0,2));
                FVDg = obj.GSD(obj.GSD(:,1)>0.0,2);
                mdl = fitlm(logDg,FVDg); %Linear fit of log
                lognewDgmax = log(obj.Dg_max)+(1.0-obj.Se_Dgmax).*mdl.Coefficients.Estimate(2);
                newDgmax = exp(lognewDgmax);
                obj.GSD = [obj.GSD;newDgmax,1.0]; %Dgmax with 1 is 10 times Dgmax (assumption)
                end
                
                %Extend to 0,0:
                obj.GSD = [0,0;obj.GSD];
                %Correct max and min FVg:
                obj.GSD(obj.GSD(:,2)==1.0,1)= min(obj.GSD(obj.GSD(:,2)==1.0,1)); %All FVg values with FVg=1 set to min Dg(FVg=1).
                obj.GSD(obj.GSD(:,2)==0.0,1)= max(obj.GSD(obj.GSD(:,2)==0.0,1)); %All FVg values with FVg=0 set to max Dg(FVg=0).
                
                %Remove duplicate values and order by
                %increasing Dg
                [~,ia,~] = unique(obj.GSD(:,1),'rows');
                obj.GSD = obj.GSD(ia,:); %->GSD
                obj.GSD = sortrows(obj.GSD,[1,2]);
                
            else
                if (not(isempty(Dgmax)||isnan(Dgmax)))
                    obj.Dg_max = Dgmax;
                else
                    obj.Dg_max = NaN;
                end
                obj.Dg_min = NaN;
            end
            
            % Update all GSD dependant properties:
            obj = obj.update_Dxx();
            obj = obj.update_texture();
        end
                
        function obj = update_Dxx(obj)
            %Update characteristics diameters
            if (not(isempty(obj.GSD)))
                obj.D10 = obj.Dg_FVg_from_GSD(0.1);
                obj.D20 = obj.Dg_FVg_from_GSD(0.2);
                obj.D30 = obj.Dg_FVg_from_GSD(0.3);
                obj.D40 = obj.Dg_FVg_from_GSD(0.4);
                obj.D50 = obj.Dg_FVg_from_GSD(0.5);
                obj.D60 = obj.Dg_FVg_from_GSD(0.6);
                obj.D80 = obj.Dg_FVg_from_GSD(0.8);
                obj.cu = obj.D60/obj.D10;
                obj.cc = (obj.D30^2)/(obj.D10*obj.D60);
            else
                obj.D10 = NaN;
                obj.D20 = NaN;
                obj.D30 = NaN;
                obj.D40 = NaN;
                obj.D50 = NaN;
                obj.D60 = NaN;
                obj.D80 = NaN;
                obj.cu =  NaN;
                obj.cc =  NaN;
            end
        end
 
        function obj = update_texture(obj)
            %Calculate and set the texture property of the soil. Also clay,
            %silt, sand and gravel content.
            if (not(isempty(obj.GSD)))
                GRAVEL = max(0.0,1.0-obj.FVg_Dg_from_GSD(0.002));
                CLAY = obj.FVg_Dg_from_GSD(0.000002)/(1.0-GRAVEL);
                SILT = (obj.FVg_Dg_from_GSD(0.00005)-CLAY)/(1.0-GRAVEL);
                SAND = (1.0-SILT-CLAY)/(1.0-GRAVEL);
                
                obj.clay   = obj.FVg_Dg_from_GSD(0.000002);                 %->clay
                obj.silt   = obj.FVg_Dg_from_GSD(0.00005)-obj.clay ;        %->silt
                obj.sand   = obj.FVg_Dg_from_GSD(0.002)-obj.clay-obj.silt ; %->sand
                obj.gravel = 1.0 - obj.sand-obj.clay-obj.silt ;             %->gravel
                
                if ((SILT + 1.5*CLAY) < 0.15)
                    obj.texture ='sand';
                elseif ((SILT + 1.5*CLAY >= 0.15) && (SILT + 2*CLAY < 0.30))
                    obj.texture ='loamy sand';
                elseif ((CLAY >= 0.07 && CLAY < 0.20) && (SAND > 0.52) && ((SILT + 2*CLAY) >= 0.30) || (CLAY < 0.07 && SILT < 0.50 && (SILT+2*CLAY)>=0.30))
                    obj.texture ='sandy loam';
                elseif ((CLAY >= 0.07 && CLAY < 0.27) && (SILT >= 0.28 && SILT < 0.50) && (SAND <= 0.52))
                    obj.texture ='loam';
                elseif ((SILT >= 0.50 && (CLAY >= 0.12 && CLAY < 0.27)) || ((SILT >= 0.50 && SILT < 0.80) && CLAY < 0.12))
                    obj.texture ='silt loam';
                elseif (SILT >= 0.80 && CLAY < 0.12)
                    obj.texture ='silt';
                elseif ((CLAY >= 0.20 && CLAY < 0.35) && (SILT < 0.28) && (SAND > 0.45))
                    obj.texture ='sandy clay loam';
                elseif ((CLAY >= 0.27 && CLAY < 0.40) && (SAND > 0.20 && SAND <= 0.45))
                    obj.texture ='clay loam';
                elseif ((CLAY >= 0.27 && CLAY < 0.40) && (SAND  <= 0.20))
                    obj.texture ='silt clay loam';
                elseif (CLAY >= 0.35 && SAND > 0.45)
                    obj.texture ='sandy clay';
                elseif (CLAY >= 0.40 && SILT >= 0.40)
                    obj.texture ='silty clay';
                elseif (CLAY >= 0.40 && SAND <= 0.45 && SILT < 0.40)
                    obj.texture ='clay';
                else
                    obj.texture ='';
                end
            else
                obj.texture=NaN;
            end
        end

        
        %%SOIL METHODS
        
        %% PREDICTIONS FROM GSD:
        %Grain Size Distribution (GSD)
        function output = FVg_Dg_from_GSD(obj,Dg)
            %Return the agregated grain volume (per unit) given the grain
            %size, from the GSD.
            fittolog = obj.options.interpolation.fit_to_log;
            if(size(obj.GSD,1)>=2)
                isinv=false;
                lder=false;
                output = func_soil.Interpolate_GSD(obj.GSD,Dg,lder,fittolog,isinv);
            else
                output = nan;
            end
            output = max(0,min(1,output));
        end
        
        function output = Dg_FVg_from_GSD(obj,FVg)
            %Return the grain size given the aggregated particle volume
            %from the GSD.
            fittolog = obj.options.interpolation.fit_to_log;
            if(size(obj.GSD,1)>=2)
                isinv=true;
                lder=false;
                output = func_soil.Interpolate_GSD(obj.GSD,FVg,lder,fittolog,isinv);
            else
                output = nan;
            end
        end       
    
        function output = fVg_Dg_from_GSD(obj,Dg)
            %Density probability function given the grain size (Derivative
            %of the GSD).
            if(size(obj.GSD,1)>=2)
                isinv=false;
                lder=true;
                output = func_soil.Interpolate_GSD(obj.GSD,Dg,lder,obj.options.fit_to_log,isinv);
            else
                output = NaN;
            end
        end
        
        %Pore size distribution and WRC from Youngs-Laplace in capillary
        %approach
        
        %Pore Size Distribution (PSD) from WRC
        
        function out_Dp = Dp_pp(obj,pp)
            %Pore size given the suction (Youngs equation).
            param_Ts=obj.options.liquid.param_Ts;
            param_costheta=obj.options.liquid.param_costheta;
            param_rhow=obj.options.liquid.param_rhow;
            param_g=obj.options.liquid.param_g;
            out_Dp = (4.*param_Ts.*param_costheta)./(param_rhow.*param_g.*pp);
        end
        
        function out_pp = pp_Dp_from_WRC(obj,Dp)
            %Suction given the pore diameter (Youngs Equation).
            param_Ts=obj.options.liquid.param_Ts;
            param_costheta=obj.options.liquid.param_costheta;
            param_rhow=obj.options.liquid.param_rhow;
            param_g=obj.options.liquid.param_g;
            out_pp = (4.*param_Ts.*param_costheta)./(param_rhow.*param_g.*Dp);
        end
    end
    
end



