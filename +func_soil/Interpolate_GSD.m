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

function [output] = Interpolate_GSD(GSD,xval,lder,logaritmic,isinv)
%func_SOIL_interpolate_GSD Interpolate values xval in the GSD.
%   If lder=1: Then return the density function(fVg), else return usual cummulated grain size distribution
%   If isinv=1: Then return the Dg given the FVg (particle size given
%   volume fraction (per unit).
%   If logarithmic=1: Then interpolation is done with axis value in logarithmic coordinates.

Data = [0,0;GSD(GSD(:,1)>0,:)]; % GSD extended to grain size of 0 at volume fraction 0%.
Data(Data(:,2)>1,2)=1; %Values of volume fractions FVg over 1 are converted to 1.
Data(Data(:,2)<0,2)=0; %The same for values below 0
dini = 0; %Initial derivative in the GSD assume to be horizontal in the minimum fraction
dend = 0; %Final derivative in the GSD assume to be horizontal in the maximum fraction

output = func_math.Interpolate_Monotonic_Rational_Cuadratic(Data,xval,dini,dend,logaritmic,true,0.5,lder,isinv);
end

