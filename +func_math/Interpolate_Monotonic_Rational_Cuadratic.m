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

function [output] = Interpolate_Monotonic_Rational_Cuadratic(Data,xv,dini,dend,plog,pinc,psup,lder,isinv,isxascending)
%Interpolate_Monotonic_Cuadratic Interpolate between Points 'Data', in xval
%with a piecewise monotonic rational cuadratic increasing (or decreasing
%pinc=false) with logaritmic axis (if plog=true) defining slope at ini
%(dini, can be empty) or slope at end (dend, can be empty).

%PREVIOUS PREPARATION OF THE DATA
xvnew = xv;
%Variable change in x if plog=1 (interpolate in logarithmic scale)
if (plog)
dini = 0;
Data = Data(Data(:,1)>=0,:);
Data(Data(:,1)==0,1) = 1E-15; %To be able to log interpolate in 0
if(~isinv) xvnew = log(xv);end
Data(:,1) = log(Data(:,1));
end

%Order data in ascending x if isxascending=false(is faster with x ordered):
if(exist('isxascending','var')&&~isxascending)
Data = sortrows(Data,1,'ascend');
end

%CALCULATION OF THE MAXIMUM BOUNDARY OF THE POINTS
valsup=Data;
nrowsb = size(Data,1);

for j=1:nrowsb
if (pinc)
valsup(j,2) = max(Data(1:j,2));
else
valsup(j,2) = max(Data(j:nrowsb,2));
end
end

%CALCULATION OF THE MINUMUM BOUNDARY OF THE POINTS
valinf=Data;
for j=1:nrowsb
if (pinc)
valinf(j,2) = min(Data(j:nrowsb,2));
else
valinf(j,2) = min(Data(1:j,2));
end
end

%CORRECTION OF NON MONOTONIC POINTS BY A PERCENTAGE (psup) BETWEEN MAXIMUM
%AND MINIMUM BOUNDARY.
Data(:,2) = psup.*valsup(:,2)+(1-psup).*valinf(:,2);

%Data is already x ascending, order to be y increasing:
if(~pinc)
    Data=flip(Data);
end


%INTERPOLATE NOW WITH PREPARED DATA:
%Interpolate with good data:
output = Interpolate_Monotonic_Cuadratic_aux(Data,xvnew,dini,dend,plog,pinc,lder,isinv);
%output = reshape(output,size(xvnew)); %Correction for when xvnew is in column.

% Undo variable change:
if (lder)
    if (plog) %multiply by log derivative
        masksup=(xvnew<=max(valsup(:,1)));
        maskinf=(xvnew>=min(valsup(:,1)));
        maskxvnotcero = (xv>1E-15);
        output(masksup.*maskinf.*maskxvnotcero==1) = output(masksup.*maskinf.*maskxvnotcero==1)./xv(masksup.*maskinf.*maskxvnotcero==1);
        output(masksup~=1) = 0; %Outside data, derivative is cero, no extrapolation.
        output(maskinf~=1) = 0; %Outside data, derivative is cero, no extrapolation.
    end   
end

if (isinv)
    if (plog)
        output = exp(output);
    end
end
%     plot(xv,output,'r-')
%     hold off

function [output] = Interpolate_Monotonic_Cuadratic_aux(Data,xv,dini,dend,plog,pinc,lder,isinver)
%Interpolate_Monotonic_Cuadratic_aux Function that interpolate prepared
%data with a monotonic rational cuadratic function.

% Initializations
nrows = size(Data,1);
nxval = length(xv);

isxinc = (Data(nrows)>Data(1));

%CALCULATE DERIVATIVES FOR POINTS: 
if (plog)
    if (pinc)
       dini = dini/exp(Data(1,1));
       dend=dend/exp(Data(nrows,1));
    else
       dini = dini/exp(Data(nrows,1));
       dend=dend/exp(Data(1,1));
    end
end

x1 = Data(2:nrows,1);
x0 = Data(1:nrows-1,1);  
y1=Data(2:nrows,2);
y0=Data(1:nrows-1,2);
y0b=Data(1:nrows-2,2);
y1b=Data(2:nrows-1,2);
y2b=Data(3:nrows,2);

hi = x1-x0; %Distance between points (xi+1-xi)
der = (y1-y0)./hi; %Slope of line between two points (fi-1-fi)/hi
di(1:nrows) = 0.; %Derivative in each point

%For the points inside: Arithmetic mean
%Armonic mean of the slopes in both segments:ponderated by fb-fa
fdelta1 = y2b-y1b; %f(2)-f(1)
fdelta0 = y1b-y0b; %f(1)-f(0)
sumfdelta = fdelta1+fdelta0;

di(2:nrows-1) = (der(1:nrows-2).*fdelta1+der(2:nrows-1).*fdelta0)./sumfdelta; %All derivatives as armonic mean
mask = [false;sumfdelta==0;false]; %When sumdelta is cero, derivative is assigned as cero.
di(mask)=0;


%For the first point, derivative will be dini or the slope in first segment
if (isempty(dini))
    if(isxinc)
        di(1) = der(1);
    else
        di(nrows) = der(nrows-1);
    end
else
    if(isxinc)
        di(1) = dini;
    else
        di(nrows) = dini;
    end
end


%For the end point:
if (isempty(dend))
 if(isxinc)
    di(nrows) = der(nrows-1);
else
    di(1) = der(1);
end

else
if(isxinc)
    di(nrows) = dend;
else
    di(1) = dend;
end
end

%Al points horizontal at left or horizontal at right has derivative cero
mask = [der==0;false];
di(mask)=0;
mask = [false;der==0];
di(mask)=0;

%INTERPOLATE POINTS
output(1:nxval)=min(Data(:,2)); %all to min value
iseg(1:nrows)=0;

%CALCULATE NUMBER OF POINTS BELOW YDATA
%As Data(:,2) is ascending:
if(isinver)
    mask = xv<Data(1,2); %below YData
else
    xyinf = min(Data(Data(:,2)==Data(1,2),1)); %Minimum x for YData
    if(pinc)
     mask = xv<xyinf; %Below xyinf    
    else
     mask = xv>xyinf; %Over xyinf  
    end
end

iseg(1) = sum(mask); %idx for xyinf

%Iteratively select segments in Y
nsegments = length(Data(:,2))-1;
for i=1:nsegments
    yinf = Data(i,2);
    ysup = Data(i+1,2);
    xyinf = Data(i,1);
    xysup = Data(i+1,1);
    iyinf = i;
    iysup = i+1;
    diyinf = di(iyinf); %derivative for yinf
    diysup = di(iysup); %derivative for ysup

    xmin = min(xyinf,xysup); %min value of x
    xmax = max(xyinf,xysup); %max value of x
    if (xmin == xyinf)
        ymin = yinf; %Value of y for xmin
        ymax = ysup; %Value of y for xmax
        dimin = diyinf; %Value of derivative for xmin
        dimax = diysup; %Value of derivative for xmin
    else
        ymin = ysup;
        ymax = yinf;
        dimin = diysup; %Value of derivative for xmin
        dimax = diyinf; %Value of derivative for xmin
    end
    

if(isinver)
    if(pinc)
    chk_overeqfirst = (xv>=ymin);
    chk_oversecond = (xv<ymax);
    chk_equallast = (xv==ymax).*(i==nsegments);
    else
    chk_overeqfirst = (xv<=ymin);
    chk_oversecond = (xv>ymax);
    chk_equallast = (xv==ymax).*(i==nsegments);        
    end
else
    chk_overeqfirst = (xv>=xmin);
    chk_oversecond = (xv<xmax);
    chk_equallast = (xv==xmax).*(i==nsegments);
end

    mask = chk_overeqfirst.*(chk_oversecond+chk_equallast); %Check if xval is in [xmin,xmax) or [xmin(nsegments),xmax(nesegments)]
    xsegment=xv(mask==1);
    iseg(i+1)=length(xsegment);
    if (lder)
        output(mask==1) = piecewise_monotonic_cuadratic_der(xsegment, xmin, xmax, ymin, ymax, dimin, dimax,pinc);
    else
        if(isinver)
            output(mask==1) = piecewise_monotonic_cuadratic_inv(xsegment, xmin, xmax, ymin, ymax, dimin, dimax,pinc);
        else
            output(mask==1) = piecewise_monotonic_cuadratic(xsegment,  xmin, xmax, ymin, ymax, dimin, dimax,pinc);
        end
    end
end 

%For the points over the data range: derivative=0, values= boundary
if(isinver)
    ymin=min(Data(:,2));
    ymax=max(Data(:,2));
    xymin=min(Data(Data(:,2)==ymin,1));
    xymax=max(Data(Data(:,2)==ymax,1));
    
    if(lder)minval=0; else minval=xymin;end
    if(lder)maxval=0; else maxval=xymax;end
    
    output(xv<=ymin) = minval;
    output(xv>=ymax) = maxval;     

else
    xmin=min(Data(:,1));
    xmax=max(Data(:,1));
    yxmin=min(Data(Data(:,1)==xmin,2));
    yxmax=max(Data(Data(:,1)==xmax,2));
    
    if(lder)minval=0; else minval=yxmin;end
    if(lder)maxval=0; else maxval=yxmax;end
 
    output(xv<=xmin) = minval;
    output(xv>=xmax) = maxval; 

end

if (size(output)~=size(xv))
output = reshape(output,size(xv));
end

end

end

%% PIECEWISE FUNCTIONS
function [output] = piecewise_monotonic_cuadratic(x,xa,xb,ffa,ffb,dda,ddb,isinc)
%piecewise_monotonic_cuadratic Piecewise monotonic rational cuadratic
%function between two points.

if(xa==xb)
    output(1:length(x)) = (dda+ddb)/2;
else

%   Detailed explanation goes here
if (isinc)
fa = ffa;
fb = ffb;
h = xb-xa; %Domain of x
HF = (fb-fa);
da = h*dda;
db = max(1E-15,h*ddb);
th = (x-xa)/h; %Relative x
else
fa = ffb;
fb = ffa;
h = xb-xa; %Domain of x
HF = (fb-fa);
th = (xb-x)/h; %Relative x
da = -h*ddb;
db = max(1E-15,-h*dda);
end


%A good valuer for m is when the inflection point is in the middle:

dalimit = 4*HF^2/(3*HF+db);
insideroot = (16*HF^4-4*HF^2*(2*HF + 5*da)*db+(-7*HF^2+14*HF*da+da^2)*db^2+4*da*db^3);

if (insideroot>0)
    if (da>dalimit)
        m=-(-4*HF^2+3*HF*db+da*db+sqrt(insideroot))/(8*HF^2-6*HF*da-2*da*db);
    else
        m=-(-4*HF^2+3*HF*db+da*db-sqrt(insideroot))/(8*HF^2-6*HF*da-2*da*db);
    end
else
    m=1;
end

if (m<0)
    m=1;
end

F = 1;
E = (m*da+db)/(fb-fa);
D = m;
B = (m*da*fb+db*fa)/(fb-fa);
C = fb;
A = m.*fa;

if ( h == 0 || fa == fb)
    output = mean([fa,fb]);
else
    p = A.*(1-th).^2+B.*(1-th).*th+C.*th.^2;
    q = D.*(1-th).^2+E.*(1-th).*th+F.*th.^2;
    output = p./q;
end
output(th<0)=fa;
output(th>1)=fb;
end
end


function [output] = piecewise_monotonic_cuadratic_der(x,xa,xb,ffa,ffb,dda,ddb,isinc)
%piecewise_monotonic_cuadratic_der Derivative of piecewise monotonic rational cuadratic
%function between two points.
if(xa==xb)
    output(1:length(x)) = (dda+ddb)/2;
else
if (isinc)
    fa = ffa;
    fb = ffb;
    h = xb-xa; %Domain of x
    HF = (fb-fa);
    da = h*dda;
    db = max(1E-15,h*ddb);
    th = (x-xa)/h; %Relative x
else
    fa = ffb;
    fb = ffa;
    h = xb-xa; %Domain of x
    HF = (fb-fa);
    th = (xb-x)/h; %Relative x
    da = -h*ddb;
    db = max(1E-15,-h*dda);
end

%A good valuer for m is when the inflection point is in the middlethat
%happens when insideroot>0 and m>0 (in other casse assign m:
dalimit = 4*HF^2/(3*HF+db);
insideroot = (16*HF^4-4*HF^2*(2*HF + 5*da)*db+(-7*HF^2+14*HF*da+da^2)*db^2+4*da*db^3);

if (insideroot>0)
    if (da>dalimit)
        m=-(-4*HF^2+3*HF*db+da*db+sqrt(insideroot))/(8*HF^2-6*HF*da-2*da*db);
    else
        m=-(-4*HF^2+3*HF*db+da*db-sqrt(insideroot))/(8*HF^2-6*HF*da-2*da*db);
    end
else
    m=1;
end

if (m<0)
    m=1;
end

F = 1;
E = (m*da+db)/(fb-fa);
D = m;
B = (m*da*fb+db*fa)/(fb-fa);
C = fb;
A = m.*fa;

if ( h == 0 || fa == fb)
    output = 0.0;
else
    dp = (B*D-A*E).*(1-th).^3+...
        (B*D+2*C*D-A*E-2*A*F).*(1-th).^2.*th+...
        (2*C*D+C*E-2*A*F-B*F).*(1-th).*th.^2+...
        (C*E-B*F).*th.^3;
    q = D.*(1-th).^2+E.*(1-th).*th+F.*th.^2;
    if(isinc)
        output = (1/h).*dp./(q.^2);
    else
        output = -(1/h).*dp./(q.^2);
    end
end
output(th<0)=0;
output(th>1)=0;
end
end

function [output] = piecewise_monotonic_cuadratic_inv(s,xa,xb,ffa,ffb,dda,ddb,isinc)
%piecewise_monotonic_cuadratic_inv Inverse of piecewise monotonic rational cuadratic
%function between two points.
if (ffa==ffb)
    output(1:length(s)) = (xa+xb)/2;
else

if (isinc)
fa = ffa;
fb = ffb;
h = xb-xa; %Domain of x
HF = (fb-fa);
da = h*dda;
db = max(1E-15,h*ddb);
else
fa = ffb;
fb = ffa;
h = xb-xa; %Domain of x
HF = (fb-fa);
da = -h*ddb;
db = max(1E-15,-h*dda);
end

if (fb<=fa)
    output = fa+0.*s;
else

%A good valuer for m is when the inflection point is in the middle:
dalimit = 4*HF^2/(3*HF+db);
insideroot = (16*HF^4-4*HF^2*(2*HF+5*da)*db+(-7*HF^2+14*HF*da+da^2)*db^2+4*da*db^3);

if (insideroot>0)
    if (da>dalimit)
        m=-(-4*HF^2+3*HF*db+da*db+sqrt(insideroot))/(8*HF^2-6*HF*da-2*da*db);
    else
        m=-(-4*HF^2+3*HF*db+da*db-sqrt(insideroot))/(8*HF^2-6*HF*da-2*da*db);
    end
else
    m=1;
end

if (m<0)
    m=1;
end

F = 1;
E = (m*da+db)/(fb-fa);
D = m;
B = (m*da*fb+db*fa)/(fb-fa);
C = fb;
A = m.*fa;

numeratorroot = B^2-4*A*C+(4*C*D).*s-(2*B*E).*s+(4*A*F).*s+(E^2).*s.^2-(4*D*F).*s.^2;
numerator1 = 2*A-B-(2*D).*s+E.*s;
denominator = 2.*(A-B+C-D.*s+E.*s-F.*s);

output = (numerator1+sqrt(abs(numeratorroot)))./denominator;
if(isinc)
output = xa + h.*output;
else
output = xb - h.*output;
end

    
end
%For the cases outside the domain, return NaN, no extrapolation admited
output(s>fb)=NaN;
output(s<fa)=NaN;
end

end


