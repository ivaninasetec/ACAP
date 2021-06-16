%-------------------------------------------------------------------------------
%MIT License
% 
% Copyright (c) 2020 Ivan Campos-Guereta Diez
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
%   'CAMPOS-GUERETA, I. et al 2020. An alternative continuous form of Arya and 
%   Paris (ACAP) model to predict the soil retention curve of a soil. MATLAB 
%   model and empirical fits. Journal of Hydrology (submitted for publication)'
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%-------------------------------------------------------------------------------

function [output] = Interpolate_WRC(WRC,xval,lder,thres,thsat,logaritmic,isinv)
%func_SOIL_interpolate_WRC Interpolate values xval in the GSD.
%   WRC: Array with the water retention curve points (suction[m],
%   volumetric water content [m3,m3])
%   xval: Array with values where to interpolate.
%   lder: If 'l' Then return the density function(fVg), else return usual cummulated grain size distribution
%   thres: Minimum admisible value for the volumetric water content
%   (residual or irreducible volumetric water content)
%   thsat: Maximum admisible value for the volumetric water content
%   (volumetric water content at saturation).
%   logarithmic: If '1' then interpolation is done with axis value in logarithmic coordinates.
%   isinv: If '1' Then return the Dg given the FVg (particle size given
%   volume fraction (per unit).

% Default value for thres=0, for thsat=max on the WRC.
if (isempty(thres)||isnan(thres))
    thres=0.0;
end

if (isempty(thsat)||isnan(thsat))
    thsat = max(WRC(:,2));
end

    WRC_PROCESSED = WRC;
    WRC_PROCESSED(WRC_PROCESSED(:,1)<1E-8,1)=1E-8; %All suction values below 1E-8 equal to 1E-8
    
    % Add the point for thsat.
    MINPP = min(WRC_PROCESSED(:,1));   %Minimum suction in the WRC
    WRC_PROCESSED = [max(1E-10,min(1E-5,MINPP./1000)),thsat;WRC_PROCESSED]; %Add a point for thsat at suction 1E_10

    % Add the point for thres
    MAXPP = max(WRC_PROCESSED(:,1));   %Max  suction in the WRC
    WRC_PROCESSED = [WRC_PROCESSED;max(1E5,MAXPP.*5),thres]; % Add a point for thres.
    
    %Select only points with volumetric water content below thsat, and with
    %suction over 1E-9 m (to avoid numerical issues)
    WRC_PROCESSED(WRC_PROCESSED(:,1)<1E-10,1)= 1E-10; %To avoid numerical issues all values of suction below 1E-10 will be included as 1E-10.
    WRC_PROCESSED(WRC_PROCESSED(:,2)>thsat,2)= thsat; %All values over thsat will be writen as thsat.
    WRC_PROCESSED(WRC_PROCESSED(:,2)<thres,2)= thres; %All values below thres will be writen as thres.

%Variable change to interpolate (inverse of suction, the same as changing
%sign of log(suction).
INVWRC = WRC_PROCESSED;
INVWRC(:,1) = 1./WRC_PROCESSED(:,1); 
INVWRC = flip(INVWRC); %Flip to put suctions ascending.

dini = 0; %Derivatives at the thres considered equal to 0
dend = 0; %Derivative at saturation considered equal to 0
if (isinv==1) 
    xvaltemp = xval; %No variable change on th.
else
    if (xval == 0)
        xvaltemp = 1./1E-10; %for the particular case of xval=0, then value is 1./1E-10 to avoid numerical problems.
    else
        xvaltemp = 1./xval; %In other case variable change for xval.
    end
end
output = real(func_math.Interpolate_Monotonic_Rational_Cuadratic(INVWRC,xvaltemp,dini,dend,logaritmic,true,0.5,lder,isinv));

if (isinv==1) %For inverse calculation, the output need to do the variable change.
    output = 1./output;
end
end



% function [output] = func_SOIL_interpolate_WRC(WRC,xval,lder,thres,thsat,logaritmic,isinv)
% %func_SOIL_interpolate_WRC Interpolate values xval in the GSD.
% %   WRC: Array with the water retention curve points (suction[m],
% %   volumetric water content [m3,m3])
% %   xval: Array with values where to interpolate.
% %   lder: If 'l' Then return the density function(fVg), else return usual cummulated grain size distribution
% %   thres: Minimum admisible value for the volumetric water content
% %   (residual or irreducible volumetric water content)
% %   thsat: Maximum admisible value for the volumetric water content
% %   (volumetric water content at saturation).
% %   logarithmic: If '1' then interpolation is done with axis value in logarithmic coordinates.
% %   isinv: If '1' Then return the Dg given the FVg (particle size given
% %   volume fraction (per unit).
% 
% % Default value for thres=0, for thsat=max on the WRC.
% 
% if isnan(thres)
%     thres=0.0;
% end
% if(isempty(thsat)||isnan(thsat))
%     mask = WRC(:,1)>0;
%     WRC_PROCESSED = WRC(mask==1,:);
%     MAXPP = max(WRC_PROCESSED(:,1)); %Max  suction in the WRC
%     WRC_PROCESSED = [WRC_PROCESSED; max(1E5,MAXPP.*5),thres]; % Add a point for thres.
%     MINPP = min(WRC_PROCESSED(WRC_PROCESSED(:,1)>0,1)); %Min  suction in the WRC
%     WRC_PROCESSED(WRC_PROCESSED(:,1)==0,1)=min(1E-3,MINPP./10); %To consider saturated state.
% else
%     mask1 = WRC(:,2)<thsat;
%     mask2 = WRC(:,1)>1E-9;
%     mask = mask1.*mask2;
%     WRC_PROCESSED = WRC(mask==1,:);
%     MINPP = min(WRC_PROCESSED(:,1)); %Min  suction in the WRC
%     MAXPP = max(WRC_PROCESSED(:,1)); %Max  suction in the WRC
%     WRC_PROCESSED = [min(1E-3,MINPP./5),min(thsat,max(WRC(:,2)));WRC_PROCESSED; max(1E5,MAXPP.*5),thres];
% end
% 
% %Variable change to interpolate (inverse of suction, the same as changing
% %sign of log(suction).
% INVWRC(:,1) = 1./WRC_PROCESSED(:,1); %Considering inverse for interpolation.
% INVWRC = flip(INVWRC);
% dini = 0;
% dend = 0;
% if (isinv==1)
%     xvaltemp = xval;
% else
%     if (xval == 0)
%         xvaltemp = 1./1E-10;
%     else
%         xvaltemp = 1./xval;
%     end
% end
% output = real(func_MATH_Interpolate_Monotonic_Rational_Cuadratic(INVWRC,xvaltemp,dini,dend,logaritmic,true,0.5,lder,isinv));
% if (isinv==1)
%     output = 1./output;
% end
% end
