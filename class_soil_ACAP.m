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

classdef class_soil_ACAP
    %class_soil_ACAP Class for the ACAP model to predict the pore size
    %distribution and the Water Retention Curve from Grain Size
    %Distribution data.
    %
    %EXAMPLE:
    % GSD = [ 0                       ,0;...
    %         2.00000000000000e-06    ,0.0200000000000000;...
    %         5.00000000000000e-05    ,0.130000000000000;...
    %         0.000106000000000000    ,0.256000000000000;...
    %         0.000250000000000000    ,0.743000000000000;...
    %         0.000500000000000000    ,0.916000000000000;...
    %         0.00100000000000000     ,0.985000000000000;...
    %         0.00200000000000000     ,0.997000000000000];
    %
    % npor = 0.45;
    %
    % WRC = [ 0,                      0.361500000000000;...
    %         0.200000000000000      ,0.324000000000000;...
    %         0.300000000000000      ,0.265500000000000;...
    %         0.400000000000000      ,0.194500000000000;...
    %         0.500000000000000      ,0.160000000000000;...
    %         0.700000000000000      ,0.111000000000000;...
    %         1                      ,0.0795000000000000;...
    %         2                      ,0.0570000000000000;...
    %         4.30000000000000       ,0.0445000000000000;...
    %         6.40000000000000       ,0.0365000000000000;...
    %         10.5000000000000       ,0.0265000000000000];
    %
    %   optionsACAP = struct(...
    % 	'fit_to_log',true,... % If true: Interpolations are done in the logDg instead on Dg
    % 	'fit_max_Se',0.95,... % Maximum value of relative saturation were to take points in the WRC for fitting purposes (to avoid horizontal line in high saturation)
    % 	'fit_min_Se',0.1,...  % Minimum value of relative saturation were to take points in the WRC for fitting purposes (to avoid horizontal line in high saturation)
    % 	'saturation_vs_porosity',1.0,... % The saturated water content will be considered as porosity multiplied by this coefficient (1.0 for thsat=npor)
    % 	'opt_beta',5,... %1:fix, 2:tex, 3:fix_fit,4:fix_emp, 5:emp_dp, 6:fit_dp
    % 	'opt_Dg0',1);... %Option Dg0=1, 2, 3, 4, 5
    %
    % soil = class_soil_ACAP(...
    %     'GSD',GSD,...
    %     'WRC',WRC,...
    %     'npor',npor,...
    %     'options',optionsACAP);
    %
    % disp(soil.FVp_Dp_from_GSD(1E-4));
    %
    % disp(soil.pp_th_from_GSD(0.35));
    %
    % INPUT PARAMETERS:
    % Mandatory:
    % GSD[]: Grain size distribution array-(Size[m],cum volume[per unit])
    % npor: Porosity index [m3/m3]
    % Optiona:
    % beta: Beta value if manually defined.
    % WRC[](optional): Water retention curve on drying array-(Suction[m], volumetric water
    % content [m3/m3])
    % Dg_max(optional): Max grain size in the sample [m] if known.
    
    %
    % OPTIONS FOR BETA:
    % opt_beta:	With beta is to be taken to calculate results
    %  0: Beta included as input (Included in property beta_set)
    %  1: Fixed (beta_fix)
    %  2: Fixed from texture (beta_fix_tex)
    %  3: Fixed fitted from WRC (beta_fix_fit)
    %  4: Fixed empirical (beta_fix_emp)
    %  5: Variable empirical dependant on Dp (beta_emp_Dp(Dp))
    %  6: Variable empirical dependant on Dp and on beta_fit_fix (beta_fit_Dp(Dp))
    %  7: Beta_fit_fix if possible, else beta_fix
    %
    % OPTIONS FOR Dg0:
    % optDg0: Option for considering Dg0
    % 1: Dg0 = D50/0.001m (Scaled diameter at D50 equivalent to a particle of 1mm
    % 2: Dg0 = D50/0.000074m (Scaled diameter at D50 equivalent to a particle of 0.074mm
    % 3: Dg0 = D50/0.000002m (Scaled diameter at D50 equivalent to a particle of 0.002mm
    % 4: Dg0 = 2650^(1/3)(Unit scaled volume equivalent to 1gr: Vg0/Dg0^3 = 1gr = 1/2650 m3, and considering unit Vg0: Dg0 = (2650.^(1/3)
    % 5: Dg0 = 1 (No Dg scaling, and values of Dg in m)
    % 6: Dg0 = Dg0 is empirically fitted at the same time as beta. For beta_fix=1.723184 -> lnDg0= (-2.734844+0.752894·lnD50-1.604248·lnD80+3.228156·npor)/(1-beta_fix);
    % Other empirical values for Dg0 selected if option_beta=2 (beta_fix_tex selected).
    %
    % OPTIONS IN optionsACAP:
    % 	'fit_to_log':(true of false): If true: Interpolations are done in the logDg instead on Dg
    % 	'fit_max_Se': (0-1): Maximum value of relative saturation were to take points in the WRC for fitting purposes (to avoid horizontal line in high saturation)
    % 	'fit_min_Se': (0-1): Minimum value of relative saturation were to take points in the WRC for fitting purposes (to avoid horizontal line in low saturation)
    % 	'saturation_vs_porosity': (double): The saturated water content will be considered as porosity multiplied by this coefficient (1.0 for thsat=npor)
    % 	'opt_beta': (int): See OPTIONS FOR BETA
    % 	'opt_Dg0': (int): See OPTIONS FOR DG0
    
    %% PARAMETERS AND OPTIONS:
    properties
        % Parameters of the wetting fluid.
        param_Ts       = 0.07197;  % Surface tension of water in 71.97 mN/m at 25ºC
        param_rhow     = 997.0479; % Density of water in kg/m3 at 25ºC
        param_g        = 9.81;     % Acceleration of gravity 9.81m/s2
        param_costheta = 1.00;     % Contact angle water-air-soil
    end
    
    %% DATA INPUTS:
    properties
        % Inputs (GSD,ht_dry,ht_wet,hk,Dg_max,npor,LL,ksat,optthetares,optthetasat,optbeta0)
        GSD=[]      %Grain Size Distribution (Size[m],Cum_vol[per unit])
        WRC=[]      %Drying Water Retention Curve measured in lab (drying) (suction[m],theta[per unit]) array
        npor=nan    %porosity index
        evoid = nan %void ratio
        beta_set = nan %beta included as input. Results calculated with this beta if opt_beta=0.
        
        options = struct('fit_to_log',true,'fit_max_Se',0.95,'fit_min_Se',0.1,'saturation_vs_porosity',1.0,'opt_beta',5,'opt_Dg0',1); %Struct with options of the ACAP model
        % 	'fit_to_log':(true of false): If true: Interpolations are done in the logDg instead on Dg
        % 	'fit_max_Se': (0-1): Maximum value of relative saturation were to take points in the WRC for fitting purposes (to avoid horizontal line in high saturation)
        % 	'fit_min_Se': (0-1): Minimum value of relative saturation were to take points in the WRC for fitting purposes (to avoid horizontal line in low saturation)
        % 	'saturation_vs_porosity': (double): The saturated water content will be considered as porosity multiplied by this coefficient (1.0 for thsat=npor)
        % 	'opt_beta': (int): See OPTIONS FOR BETA
        % 	'opt_Dg0': (int): See OPTIONS FOR DG0
    end
    
    %% PROPERTIES CALCULATED DURING CONSTRUCTION (READONLY)
    properties (SetAccess = private)
        %GSD calculated properties (read only):
        clay=nan    %Per unit clay content (<0.000002m)
        silt=nan    %Per unit silt content (<0.00005m)
        sand=nan    %Per unit sand content (
        gravel=nan  %Per gravel content content
        texture=nan %Texture name for the soil
        cu=nan      %coefficient of uniformity
        cc=nan      %coefficient of curvature
        Dg_max=NaN  %Max Diameter from GSD
        Dg_min=NaN  %Min Diameter from GSD
        
        Se_Dgmax  = nan %Relative saturation at Dgmax
        Se_Dgmin  = nan %Relative saturation at Dgmin
        Se_ppmin  = nan %Relative saturation at ppmin
        Se_ppmax  = nan %Relative saturation at ppmax
        Se_minfit = nan %Relative saturation at Se as defined in min fit
        Se_maxfit = nan %Relative saturation at Se as defined in max fit
        
        D10=nan     %Characteristic grain diameter for FVg=10%.
        D20=nan     %Characteristic grain diameter for FVg=20%.
        D30=nan     %Characteristic grain diameter for FVg=30%.
        D40=nan     %Characteristic grain diameter for FVg=40%.
        D50=nan     %Characteristic grain diameter for FVg=50%.
        D60=nan     %Characteristic grain diameter for FVg=60%.
        D80=nan     %Characteristic grain diameter for FVg=80%.
        
        DP10=nan     %Characteristic pore diameter for FVp=10%.
        DP20=nan     %Characteristic pore diameter for FVp=20%.
        DP30=nan     %Characteristic pore diameter for FVp=30%.
        DP40=nan     %Characteristic pore diameter for FVp=40%.
        DP50=nan     %Characteristic pore diameter for FVp=50%.
        DP60=nan     %Characteristic pore diameter for FVp=60%.
        DP80=nan     %Characteristic pore diameter for FVp=80%.
        
        %WRC calculated values:
        ppmin_fromWRC = nan %Minimum suction on the WRC
        ppmax_fromWRC = nan %Maximum suction on the WRC
        ppmin_from_maxSe = nan %Minimum suction from Se_maxfit
        ppmax_from_minSe = nan %Maximum suction from Se_minfit
        
        %Fitted beta from WRC data if available:
        beta_fix_fit        %Fixed value of beta fitted from ht_dry
        beta_fix_fit_RMSE   %Error of the fit
        beta_fix_fit_R2     %R2 of the fit
        
        % Table with results at WRC points (if available)
        results_at_WRC_points ; %Array with values of betai for every point of the WRC.[pp,betai]
        
    end
    
    %% PROPERTIES EVALUATED ON THE FLY BY THE GET METHOD:
    properties (Dependent)
        Dg0                 % Parameter defined from the optDg. To make adimensional the equations.
        beta_fix            % Fixed value of beta for all soils
        beta_fix_tex        % Fixed value of beta from each soil texture
        beta_fix_emp        % Fixed value of beta from empirical expr. 1
        thsat               % Volumetric water content at saturation
        thres               % Irreducible Volumetric water content (thres=0.0 in this model, but the shape of WRC could give other pseudo thres)
    end
    
    %% CALCULATIONS OF BETA
    properties (Access = private)
        %beta_tex_DG0: beta_fix_tex values for each soil texture and each Dg0 option
        %With median values:
        %   1                2                       3                       4                       5                       6                  7                    8                   9
        %C+SIC               SICL+CL                 SIL+SI                  SC+SCL                  L                       LS                 S                    SL                  OTHER
        beta_tex_DG0 = [1.39983702366954,1.06615021357796,1.17308234628200,1.45098384185803,1.28473918430092,1.17837441078356,1.17829228662542,1.25745825018460   ,1.20419957400611;...%Option_Dg0=1
            1.26559553276568,1.04842271132817,1.12346731627204,1.31905649880886,1.19546015589062,1.13243692318288,1.13211461997146,1.18665381529653   ,1.14824911332581;...%Option_Dg0=2
            1.18140930819496,1.03520386021952,1.08806801538045,1.22580487213817,1.13616953326054,1.09710434444152,1.09736967461926,1.13483073596038   ,1.10758470240943;...%Option_Dg0=3
            1.14783352907219,1.03169662100750,1.08551092394768,1.23865964115414,1.13961418342453,1.11201217619984,1.12104963199143,1.15175790420880   ,1.12017079390749;...%Option_Dg0=4
            1.18204391938494,1.03854221179243,1.10744056348961,1.30971826756314,1.18027032220174,1.14432996002729,1.15904969027968,1.19648274189929   ,1.15477471729463;...%Option_Dg0=5
            1.165660        ,1.984559        ,2.153743        ,2.090194        ,1.865421        ,1.656472        ,1.908166        ,1.685937           ,1.723184]%Option_Dg0=6
        
        %With mean values:
        %                     %   1                2                       3                       4                       5                       6                  7                    8                   9
        %                     %C+SIC               SICL+CL                 SIL+SI                  SC+SCL                  L                       LS                 S                    SL                  OTHER
        %     beta_tex_DG0 = [1.43119134201493,1.09882299905252,1.16627797702636,1.50313911580693,1.28747630135579,1.16449678975110,1.22946545439922,1.28737251679139 ,1.24152437887249;...
        %                     1.27946289753972,1.06741959679819,1.12112305695728,1.35078486114769,1.20661496349190,1.12324342666633,1.16884231784166,1.20826603422045 ,1.17456841991024;...
        %                     1.18800138345415,1.04672200562951,1.08821632527422,1.24760217782631,1.14903489777766,1.09109234875706,1.12376206905618,1.15100425415386 ,1.12634392869167;...
        %                     1.15239028643728,1.04262889890627,1.08725871475019,1.29169745400462,1.15040347604031,1.10576121118971,1.15850579201811,1.16820661150259 ,1.14462325573752;...
        %                     1.18892077280415,1.05364164190862,1.10838112700808,1.39299825777708,1.18915511842778,1.13582977098161,1.21202249063378,1.21652169078297 ,1.18872115269961;...
        %                     1.165660            ,1.984559               ,2.153743               ,2.090194               ,1.865421               ,1.656472           ,1.908166           ,1.685937           ,1.723184]       
        
    end
    
    methods
        function out = get.thres(obj)
            %Return value of irreducible volumetric water content
            % (residual saturation) in this model is 0 always.
            out = 0.0;
        end
        
        function out = get.thsat(obj)
            %Return value of the volumetric water content at saturation
            % in this model is equal to the porosity.
            out = obj.options.saturation_vs_porosity.*obj.npor;
        end
        
        function out = get.beta_fix(obj)
            %Return the value of beta_fix
            switch obj.options.opt_Dg0 %-> Dg0
                case 1 % Scaled diameter at D50 equivalent to particle of 1mm: (D50/Dg0 = 0.001m)
                    out = obj.beta_tex_DG0(1,9);
                case 2
                    out = obj.beta_tex_DG0(2,9);
                case 3
                    out = obj.beta_tex_DG0(3,9);
                case 4
                    out = obj.beta_tex_DG0(4,9);
                case 5
                    out = obj.beta_tex_DG0(5,9);
                case 6
                    out = obj.beta_tex_DG0(6,9);
                otherwise
                    out = NaN;
            end
        end
        
        function out = get.beta_fix_tex(obj)
            %Return the value of beta_fix_fit (depending on soil texture)
            switch obj.texture
                case 'sand'
                    ntex=7;
                case 'loamy sand'
                    ntex=6;
                case 'sandy loam'
                    ntex=8;
                case 'loam'
                    ntex=5;
                case 'silt loam'
                    ntex=3;
                case 'silt'
                    ntex=3;
                case 'sandy clay loam'
                    ntex=4;
                case 'clay loam'
                    ntex=2;
                case 'silt clay loam'
                    ntex=2;
                case 'sandy clay'
                    ntex=4;
                case 'silty clay'
                    ntex=1;
                case 'clay'
                    ntex=1;
                otherwise
                    ntex=9;
            end
            out = obj.beta_tex_DG0(obj.options.opt_Dg0,ntex);
        end
        
        function out = get.beta_fix_emp(obj)
            %Return empirical beta_fix_emp, depending on characteristics
            %diameters (D10, D50, D80), porosity and coefficient of
            %uniformity.
            switch obj.options.opt_Dg0 %-> Dg0
                case 1 % Scaled diameter at D50 equivalent to particle of 1mm: (D50/Dg0 = 0.001m)
                    out =+1.287987+0.513320.*log(obj.D10)-0.506320.*log(obj.D50)-0.033267.*log(obj.D80)+0.564173.*log(obj.cu)-0.325973.*obj.npor-0.004290.*log(obj.D10).*log(obj.D80); %With all points of the WRC
                    %out =+9.746296e-01-0.326924.*obj.npor+0.555587.*log(obj.D10)-0.563082.*log(obj.D50)-0.077482.*log(obj.D80)+0.624559.*log(obj.cu)-0.006612.*log(obj.D10).*log(obj.D80); %With all samples
                    
                case 2
                    out =+1.239192+0.340354.*log(obj.D10)-0.335174.*log(obj.D50)-0.019311.*log(obj.D80)+0.375547.*log(obj.cu)-0.220306.*obj.npor-0.002935.*log(obj.D10).*log(obj.D80);%With all points of the WRC
                    %out =+1.031201e+00-0.209883.*obj.npor+0.373206.*log(obj.D10)-0.374234.*log(obj.D50)-0.051523.*log(obj.D80)+0.420235.*log(obj.cu)-0.004462.*log(obj.D10).*log(obj.D80);%With all samples
                case 3
                    out =+1.189132+0.232807.*log(obj.D10)-0.228772.*log(obj.D50)-0.011646.*log(obj.D80)+0.257375.*log(obj.cu)-0.151684.*obj.npor-0.002031.*log(obj.D10).*log(obj.D80);%With all points of the WRC
                    %out =+1.048097e+00-0.138280.*obj.npor+0.257604.*log(obj.D10)-0.255703.*log(obj.D50)-0.035285.*log(obj.D80)+0.290136.*log(obj.cu)-0.003062.*log(obj.D10).*log(obj.D80);%With all samples
                case 4
                    out =+1.537278+0.271187.*log(obj.D10)-0.246088.*log(obj.D50)+0.017847.*log(obj.D80)+0.288340.*log(obj.cu)-0.136281.*obj.npor-0.000910.*log(obj.D10).*log(obj.D80);%With all points of the WRC
                    %out  =+1.347031e+00-0.129250.*obj.npor+0.285986.*log(obj.D10)-0.266319.*log(obj.D50)-0.011757.*log(obj.D80)+0.313904.*log(obj.cu)-0.002278.*log(obj.D10).*log(obj.D80);%With all samples
                case 5
                    out =+1.854601+0.378810.*log(obj.D10)-0.333862.*log(obj.D50)+0.036696.*log(obj.D80)+0.395870.*log(obj.cu)-0.161140.*obj.npor-0.000472.*log(obj.D10).*log(obj.D80);%With all points of the WRC
                    %out =+1.570625e+00-0.166580.*obj.npor+0.387283.*log(obj.D10)-0.355202.*log(obj.D50)-0.003142.*log(obj.D80)+0.420264.*log(obj.cu)-0.002492.*log(obj.D10).*log(obj.D80);%With all samples
                case 6
                    out = 1.723184;
            end
        end
               
        function [out_beta] = beta_emp_pp(obj,pp)
            %Return empirical beta_emp_Dp for a given suction.
            Dp = obj.Dp_pp_from_WRC(pp);
            out_beta = obj.beta_emp_Dp(Dp);
        end
        
        function [out_beta] = beta_emp_Dp(obj,Dp)
            %Return empirical beta_emp_Dp for a given pore size
            switch obj.options.opt_Dg0
                case 1
                    out_beta=-0.554079+0.807073.*obj.beta_fix_emp-0.067649.*log(Dp)-0.002852.*log(Dp).*obj.beta_fix_emp;
                case 2
                    out_beta=-0.395661+0.827978.*obj.beta_fix_emp-0.051617.*log(Dp)-0.001702.*log(Dp).*obj.beta_fix_emp;
                case 3
                    out_beta =-0.252411+0.829195.*obj.beta_fix_emp-0.037486.*log(Dp)-0.002146.*log(Dp).*obj.beta_fix_emp;
                case 4
                    out_beta =-0.227131+0.762048.*obj.beta_fix_emp-0.025756.*log(Dp)-0.016651.*log(Dp).*obj.beta_fix_emp;
                case 5
                    out_beta =-0.388402+0.788209.*obj.beta_fix_emp-0.038649.*log(Dp)-0.015805.*log(Dp).*obj.beta_fix_emp;
                case 6
                    out_beta = 1.723184;
            end
        end
        
        function [out_beta] = beta_emp_Dg(obj,Dg)
            %Return empirical beta_emp_Dp for a given grain size
            C1 = log(sqrt(3.0/2.0)*pi/4.0)/3.0-log(obj.evoid)./2.0-log(obj.Dg0);
            C2 = log(pi/6.0)/3.0-log(obj.Dg0);
            K2 = C2+log(Dg);
            
            switch obj.options.opt_Dg0
                case 1
                    A = -0.554079+0.807073.*obj.beta_fix_emp;
                    B = -0.067649-0.002852.*obj.beta_fix_emp;
                    
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 2
                    A = -0.395661+0.827978.*obj.beta_fix_emp;
                    B = -0.051617-0.001702.*obj.beta_fix_emp;
                    out_beta =  C1./K2+(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 3
                    A = -0.252411+0.829195.*obj.beta_fix_emp;
                    B = -0.037486-0.002146.*obj.beta_fix_emp;
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 4
                    A = -0.227131+0.762048.*obj.beta_fix_emp;
                    B = -0.025756-0.016651.*obj.beta_fix_emp;
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 5
                    A = -0.388402+0.788209.*obj.beta_fix_emp;
                    B = -0.038649-0.015805.*obj.beta_fix_emp;
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 6
                    out_beta = 1.723184;
            end            
            
        end
        
        function [out_beta] = beta_emp_Se(obj,Se)
            %Return empirical beta_emp_Dp for a given volume fraction
            %(relative saturation, or cummulated pore volume, or cummulated
            %grain volume)
            Dgtemp = obj.Dg_FVg_from_GSD(Se);
            out_beta = obj.beta_emp_Dg(Dgtemp);
        end
        
        function [out_beta] = beta_fit_Dp(obj,Dp,betafitfix)
            %Return beta_fit_Dp given pore size, from fitted betafitfix or from defined beta_set (if beta is previously calibrated)
            if(not(isnan(obj.beta_set)))
                betafitfix=obj.beta_set;
            end
            switch obj.options.opt_Dg0
                
                case 1, out_beta  =-0.092664+0.672826.*betafitfix-0.029338.*log(Dp)-0.011744.*betafitfix.*log(Dp);
                case 2, out_beta  =+0.063253+0.619926.*betafitfix-0.016464.*log(Dp)-0.014795.*betafitfix.*log(Dp);
                case 3, out_beta  =+0.167522+0.594918.*betafitfix-0.007437.*log(Dp)-0.016102.*betafitfix.*log(Dp);
                case 4, out_beta  =+0.232958+0.546931.*betafitfix+0.005403.*log(Dp)-0.026727.*betafitfix.*log(Dp);
                case 5, out_beta  =+0.159121+0.575382.*betafitfix+0.001285.*log(Dp)-0.026652.*betafitfix.*log(Dp);
                    
                case 6, out_beta = 1.723184; %r2 = 0.5288, rmse= 0.70996
            end
            
        end
        
        function [out_beta] = beta_fit_pp(obj,pp,betafitfix)
            %Return beta_fit_Dp given the suction, from fitted betafitfix or from defined beta_set (if beta is previously calibrated)
            Dp = obj.Dp_pp_from_WRC(pp);
            out_beta = obj.beta_fit_Dp(Dp,betafitfix);
        end
        
        function [out_beta] = beta_fit_Dg(obj,Dg,betafitfix)
            %Return beta_fit_Dg given the grain size, from fitted betafitfix or from defined beta_set (if beta is previously calibrated)
            if(not(isnan(obj.beta_set)))
                betafitfix=obj.beta_set;
            end
            C1 = log(sqrt(3.0/2.0)*pi/4.0)/3.0-log(obj.evoid)./2.0-log(obj.Dg0);
            C2 = log(pi/6.0)/3.0-log(obj.Dg0);
            K2 = C2+log(Dg);
            
            switch obj.options.opt_Dg0
                case 1
                    A = -0.092664+0.672826.*betafitfix;
                    B = -0.029338-0.011744.*betafitfix;
                    
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 2
                    A = +0.063253+0.619926.*betafitfix;
                    B = -0.016464-0.014795.*betafitfix;
                    out_beta =  C1./K2+(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 3
                    A = +0.167522+0.594918.*betafitfix;
                    B = -0.007437-0.016102.*betafitfix;
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 4
                    A = +0.232958+0.546931.*betafitfix;
                    B = +0.005403-0.026727.*betafitfix;
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 5
                    A = +0.159121+0.575382.*betafitfix;
                    B = +0.001285-0.026652.*betafitfix;
                    out_beta =  C1./K2 +(K2.*A-C1)./ (K2.* (1 - K2.* B));
                case 6
                    out_beta = 1.723184;
            end
            
        end
        
        function [out_beta] = beta_fit_Se(obj,Se,betafix)
            %Return beta_fit_Dg given the relative saturation (or the cummulated pore volume or cummulated grain volume), from fitted betafitfix or from defined beta_set (if beta is previously calibrated)
            Dgtemp = obj.Dg_FVg_from_GSD(Se);
            out_beta = obj.beta_fit_Dg(Dgtemp,betafix);
        end
        
        function [beta] = beta_Dg(obj,varargin)
            %Returns beta, given Dg. Select the output depending on
            %options.opt_beta.
            
            switch obj.options.opt_beta
                case 0
                    beta = obj.beta_set;
                case 1 %beta_fix
                    beta = obj.beta_fix;
                case 2 %beta_fix_tex
                    beta = obj.beta_fix_tex;
                case 3 %beta_fix_fit
                    if(not(isnan(obj.beta_set)))
                        beta=obj.beta_set;
                    else
                        if(isnan(obj.beta_fix_fit))
                            beta = obj.beta_fix;
                        else
                            beta = obj.beta_fix_fit;
                        end
                    end
                case 4 %beta_fix_emp
                    beta = obj.beta_fix_emp;
                case 5 %beta_emp_dp(dp(dg))
                    if (nargin==0)
                        beta = NaN;
                    else
                        beta = obj.beta_emp_Dg(varargin{1}(:,:));
                    end
                case 6 %beta_fit_dp(dp(dg))
                    if (nargin==0)
                        beta = NaN;
                    else
                        beta = obj.beta_fit_Dg(varargin{1}(:,:),obj.beta_fix_fit);
                    end
                case 7 %beta_fix_fit if possible, else beta_fix
                    if(isnan(obj.beta_fix_fit))
                        beta = obj.beta_fix;
                    else
                        beta = obj.beta_fix_fit;
                    end
                    
                otherwise
                    beta = NaN;
            end
        end
        
        function [beta] = beta_Dp(obj,varargin)
            %Returns beta, given Dp. Select the output depending on
            %options.opt_beta.
            switch obj.options.opt_beta
                case 0
                    beta = obj.beta_set;
                case 1 %beta_fix
                    beta = obj.beta_fix;
                case 2 %beta_fix_tex
                    beta = obj.beta_fix_tex;
                case 3 %beta_fix_fit
                    if(not(isnan(obj.beta_set)))
                        beta=obj.beta_set;
                    else
                        if(isnan(obj.beta_fix_fit))
                            beta = obj.beta_fix;
                        else
                            beta = obj.beta_fix_fit;
                        end
                    end
                case 4 %beta_fix_emp
                    beta = obj.beta_fix_emp;
                case 5 %beta_emp_dp(dp(dg))
                    if (nargin==0)
                        beta = NaN;
                    else
                        beta = obj.beta_emp_Dp(varargin{1}(:,:));
                    end
                case 6 %beta_fit_dp(dp(dg))
                    if (nargin==0)
                        beta = NaN;
                    else
                        beta = obj.beta_fit_Dp(varargin{1}(:,:),obj.beta_fix_fit);
                    end
                case 7 %beta_fix_fit if possible, else beta_fix
                    if(isnan(obj.beta_fix_fit))
                        beta = obj.beta_fix;
                    else
                        beta = obj.beta_fix_fit;
                    end
                otherwise
                    beta = NaN;
            end
        end
        
        function [beta] = beta_pp(obj,varargin)
            %Returns beta, given the suction. Select the output depending on
            %options.opt_beta.
            if (narging==0)
                beta=obj.beta_Dp();
            else
                Dp = obj.Dp_pp_from_WRC(varargin(1));
                beta=obj.beta_Dp(Dp);
            end
        end
        
        function [beta] = beta_Se(obj,varargin)
            %Returns beta given the relative saturation (or the cummulated pore volume or cummulated grain volume). Select the output depending on
            %options.opt_beta.
            if (narging==0)
                beta=obj.beta_Dg();
            else
                Dg = obj.Dg_FVg_from_GSD(varargin(1));
                beta=obj.beta_Dg(Dg);
            end
        end
        
    end
    
    
    %% VALUES FOR DG0
    methods
        function out = get.Dg0(obj)
            %Return the value of Dg0 depending on options.opt_Dg0.
            switch obj.options.opt_Dg0 %-> Dg0
                case 1 % Scaled diameter at D50 equivalent to particle of 1mm: (D50/Dg0 = 0.001m)
                    out = obj.D50./0.001;
                case 2 % Scaled diameter at D50 equivalent to particle of 1mm: (D50/Dg0 = 0.000074m)
                    out = obj.D50./0.000074;
                case 3 % Scaled diameter at D50 equivalent to particle of 1mm: (Dg/D50 = 0.000002m)
                    out = obj.D50./0.000002;
                case 4 % Unit scaled volume equivalent to 1gr: Vg0/Dg0^3 = 1gr = 1/2650 m3, and considering unit GSD: Dg0 = (2650.^(1/3):
                    out = 2650.^(1/3);
                case 5 % No Dg scaling, values of Dg in m:
                    out = 1.0;
                case 6
                    if (obj.options.opt_beta == 2)
                        switch obj.texture % Fixed value of beta from texture
                            case 'sand'
                                out = exp((-2.746228e+00-1.193779.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'loamy sand'
                                out = exp((+7.955937e-01-0.047215.*log(obj.D50)-0.415815.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'sandy loam'
                                out = exp((+5.199648e+00+1.827450.*log(obj.D50)-1.918128.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'loam'
                                out = exp((-7.115535e+00+0.654145.*log(obj.D20)-0.178228.*log(obj.D50)-2.333656.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'silt loam'
                                out = exp((-9.620350e+00-3.005273.*log(obj.D50)+1.145536.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'silt'
                                out = exp((-9.620350e+00-3.005273.*log(obj.D50)+1.145536.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'sandy clay loam'
                                out = exp((-8.988035e+00+1.844248.*log(obj.D50)-4.028346.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'clay loam'
                                out = exp((-2.447604e+01-3.393106.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'silt clay loam'
                                out = exp((-2.447604e+01-3.393106.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'sandy clay'
                                out = exp((-8.988035e+00+1.844248.*log(obj.D50)-4.028346.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'silty clay'
                                out = exp((-3.827717e-01+2.379576.*log(obj.D50)-2.918302.*log(obj.D80))/(1-obj.beta_fix_tex));
                            case 'clay'
                                out = exp((-3.827717e-01+2.379576.*log(obj.D50)-2.918302.*log(obj.D80))/(1-obj.beta_fix_tex));
                            otherwise
                                out = exp((-2.734844+0.752894.*log(obj.D50)-1.604248.*log(obj.D80)+3.228156.*obj.npor)./(1-obj.beta_fix));
                        end
                    else
                        out = exp((-2.734844+0.752894.*log(obj.D50)-1.604248.*log(obj.D80)+3.228156.*obj.npor)./(1-obj.beta_fix));
                    end
                    
                otherwise
                    out=NaN;
            end
        end
    end
    
    
    %% METHODS TO UPDATE PROPERTIES
    methods
        function obj = GSD_update(obj,GSD_INPUT,Dgmax)
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
                obj.Se_Dgmin = max(obj.GSD(obj.GSD(:,1)==obj.Dg_min,2));
                obj.Se_Dgmax = min(obj.GSD(obj.GSD(:,1)==obj.Dg_max,2));
                
                %Extend GSD to (Dg_max,1)
                logDg = log(obj.GSD(obj.GSD(:,1)>0.0,1));
                lobFVDg = log(obj.GSD(obj.GSD(:,1)>0.0,2));
                mdl = fitlm(logDg,lobFVDg); %Linear fit of log
                lognewDgmax = log(obj.Dg_max)+(1.0-obj.Se_Dgmax).*mdl.Coefficients.Estimate(2);
                newDgmax = exp(lognewDgmax);
                obj.GSD = [obj.GSD;newDgmax,1.0]; %Dgmax with 1 is 10 times Dgmax (assumption)
                
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
                if (not(isempty(Dgmax)||Dgmax))
                    obj.Dg_max = Dgmax;
                else
                    obj.Dg_max = NaN;
                end
                obj.Dg_min = NaN;
            end
            
            % Update all GSD dependant properties:
            obj = obj.Dxx_update();
            obj = obj.texture_update();
        end
        
        function obj = WRC_update(obj,WRC_INPUT,thsat,thres)
            %Check the Water Retention Curve and set some properties
            %dependin on the WRC.
            %% Fill WRC data, make some transformations and calc ppmin_fromWRC and ppmax_fromWRC
            if (not(isempty(WRC_INPUT)))
                WRC_local = WRC_INPUT;
                %Delete all NaN
                WRC_local= WRC_local(~isnan(WRC_local(:,1)),:);
                WRC_local= WRC_local(~isnan(WRC_local(:,2)),:);
                %All suction values below 0 equal to 0:
                WRC_local(WRC_local(:,1)<=0.0,1)= 0.0;
                WRC_local(WRC_local(:,2)==0.0,2)= max(WRC_local(WRC_local(:,2)==0.0,2)); %th at suction = 1E-5 equal to max th.
                
                %update thsat:
                if (isempty(thsat)||isnan(thsat))
                    thsat_local = max(WRC_local(WRC_local(:,1)==obj.ppmax_fromWRC,2));
                else
                    thsat_local = thsat;
                end
                
                %update thres: (default=0.0)
                if (isempty(thres)||isnan(thres))
                    thres_local = 0.0;
                else
                    thres_local = thres;
                end
                
                %Remove duplicate values and order by
                %increasing Dg
                [~,ia,~] = unique(WRC_local(:,1),'rows');
                WRC_local = WRC_local(ia,:); %->GSD
                WRC_local= sortrows(WRC_local,[1,2]);
                
                obj.WRC = WRC_local;
                
                obj.ppmin_fromWRC = min(WRC_local(:,1)); %Minumum suction in WRC
                obj.ppmax_fromWRC = max(WRC_local(:,1)); %Max suction
                %Max and min Se:
                thmin = max(WRC_local(WRC_local(:,1)==obj.ppmin_fromWRC,2));
                thmax = min(WRC_local(WRC_local(:,1)==obj.ppmax_fromWRC,2));
                obj.Se_ppmin = max(0.0,min(1.0,(thmin-obj.thres)./(thsat_local-thres_local)));
                obj.Se_ppmax = max(0.0,min(0.0,(thmax-obj.thres)./(thsat_local-thres_local)));
                
                %Update characteristic pore sizes:
                obj = obj.Dpxx_update();
            end
            
        end
        
        
        function obj = texture_update(obj)
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
        
        function obj = Dxx_update(obj)
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
        
        
        function obj = boundaries_fittings_update(obj)
            %Calculate some boundaries to define the intervals where to fit
            %measured WRC to predicted WRC. Defined by percentages   obj.options.fit_min_Se and options.fit_max_Se.
            obj.Se_minfit = obj.Se_ppmax+obj.options.fit_min_Se*(obj.Se_ppmin-obj.Se_ppmax);
            obj.Se_maxfit = obj.Se_ppmax+obj.options.fit_max_Se*(obj.Se_ppmin-obj.Se_ppmax);
            obj.ppmin_from_maxSe = obj.pp_Se_from_WRC(obj.Se_minfit);
            obj.ppmax_from_minSe = obj.pp_Se_from_WRC(obj.Se_maxfit);
            
            
        end
        
        function obj = Dpxx_update(obj)
            %Set characteristics diameters of the soil.
            if(not(isempty(obj.WRC))) %->WRC
                obj.DP10 = obj.Dp_FVp_from_WRC(0.1);
                obj.DP20 = obj.Dp_FVp_from_WRC(0.2);
                obj.DP30 = obj.Dp_FVp_from_WRC(0.3);
                obj.DP40 = obj.Dp_FVp_from_WRC(0.4);
                obj.DP50 = obj.Dp_FVp_from_WRC(0.5);
                obj.DP60 = obj.Dp_FVp_from_WRC(0.6);
                obj.DP80 = obj.Dp_FVp_from_WRC(0.8);
            else
                obj.DP10 = NaN;
                obj.DP20 = NaN;
                obj.DP30 = NaN;
                obj.DP40 = NaN;
                obj.DP50 = NaN;
                obj.DP60 = NaN;
                obj.DP80 = NaN;
            end
        end
        
        function obj = beta_fix_fit_update(obj)
            %Calculate, fit and set the value of the property beta_fix_fit from the data
            %included in the WRC (if WRC exist).
            [obj.beta_fix_fit,obj.beta_fix_fit_RMSE,obj.beta_fix_fit_R2] = obj.get_beta_fit_Se_from_WRC_loglinear();
        end
        
    end
    
    methods
        %% CONSTRUCTOR:
        function obj = class_soil_ACAP(varargin)
            %Class constructor fo ACAP model.
            %Example:
            %soil=class_soil_ACAP('GSD',GSD,'npor',npor);
            %
            %Required arguments are:
            %GSD: the grain size distribution
            %npor: the porosity.
            %
            %Optional parameters are:
            %WRC: Water Retention Curve
            %Dg_max: Maximum grain size.
            %options: Options as defined in the begining of this code.
            
            %Read the inputs:
            def_GSD     = obj.GSD;
            def_npor    = obj.npor;
            
            def_beta_set= obj.beta_set;
            def_WRC     = obj.WRC;
            def_Dg_max  = obj.Dg_max;
            def_options = obj.options;
            
            
            p=inputParser;
            addOptional(p,'GSD',def_GSD);
            addOptional(p,'WRC',def_WRC);
            addOptional(p,'Dg_max',def_Dg_max);
            addOptional(p,'npor',def_npor);
            addOptional(p,'options',def_options);
            addOptional(p,'beta',def_beta_set);
            
            parse(p,varargin{:});
            
            obj.options = p.Results.options;
            
            %Update void ratios
            obj.npor    = p.Results.npor;
            obj.evoid   = p.Results.npor/(1-p.Results.npor);
            obj.beta_set = p.Results.beta;
            
            %Update GSD and Dg_max
            obj = obj.GSD_update(p.Results.GSD,p.Results.Dg_max);
            
            %Updata WRC_data
            obj = obj.WRC_update(p.Results.WRC,obj.thsat,obj.thres);
            
            %Calculate the boundaries where to fit the data:
            obj = obj.boundaries_fittings_update();

            %Calculate beta_fix_fit by fitting to WRC
            obj = obj.beta_fix_fit_update();
            
            % Get a table with fit results in same points as the WRC.
            obj.results_at_WRC_points =  obj.get_results_at_WRC_points();
        end
        
        
        %% PREDICTIONS BASED ON THE GRAIN SIZE DISTRIBUTION
        
        % Grain Size Distribution GSD)-------------------------------------
        
        function output = FVg_Dg_from_GSD(obj,Dg)
            %Return the agregated grain volume (per unit) given the grain
            %size, from the GSD.
            if(size(obj.GSD,1)>=2)
                isinv=false;
                lder=false;
                output = func_soil.Interpolate_GSD(obj.GSD,Dg,lder,obj.options.fit_to_log,isinv);
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
        
        function output = Dg_FVg_from_GSD(obj,FVg)
            %Return the grain size given the aggregated particle volume
            %from the GSD.
            if(size(obj.GSD,1)>=2)
                isinv=true;
                lder=false;
                output = func_soil.Interpolate_GSD(obj.GSD,FVg,lder,obj.options.fit_to_log,isinv);
            else
                output = nan;
            end
        end
        
        % Pore Size Distribution (PSD)-------------------------------------
        
        function out_Dp = Dp_Dg_from_GSD(obj,Dg)
            %Return the pore size linked to the grain size (in the ACAP
            %model.
            beta = obj.beta_Dg(Dg);
            out_Dp = obj.Dp_Dg_beta_from_GSD(Dg,beta);
        end
        
        function out_Dp = Dp_Dg_beta_from_GSD(obj,Dg,beta)
            %Return the pore size linked to the grain size (in the ACAP
            %model for the given beta value.
            out_Dp = 2.^(5/6-beta./3) .* 3.^(-1/6-beta./3) .* sqrt(obj.evoid) .* pi.^((beta-1)./3) .* (Dg./obj.Dg0).^beta .* obj.Dg0;
        end      
        
        function output = FVp_Dp_from_GSD(obj,Dp)
            %Predicted Cummulated Pore Volume for a given pore diameter.
            %(Aggregated Pore Size Distribution) (predicted from the GSD and porosity)
            output = obj.FVg_Dg_from_GSD(obj.Dg_Dp_from_WRC(Dp));
        end
        
        function out_FVp = FVp_Dp_beta_from_GSD(obj,Dp,beta)
            %Predicted Cummulated Pore Volume for a given pore diameter. (predicted from the GSD and porosity)
            %(Aggregated Pore Size Distribution) given the beta value.
            Dg = obj.Dg_Dp_beta_from_WRC(Dp,beta);
            out_FVp = obj.FVg_Dg_from_GSD(Dg);
        end
        
        function out_Dp = Dp_FVp_from_GSD(obj,FVp)
            %Predicted Pore size given the Cummulated Pore Volume (predicted from the GSD and porosity).
            Dg = obj.Dg_FVg_from_GSD(FVp);
            out_Dp = obj.Dp_Dg_from_GSD(Dg);
        end
        
        function out_Dp = Dp_FVp_beta_from_GSD(obj,FVp,beta)
            %Predicted Pore size given the Cummulated Pore Volume
            %(predicted from the GSD and porosity), for a given beta value.
            Dg = obj.Dg_FVg_from_GSD(FVp);
            out_Dp = obj.Dp_Dg_beta_from_GSD(Dg,beta);
        end
        
        % Water Retention Curve (WRC)-------------------------------------
        
        function output = Se_pp_from_GSD(obj,pp)
            %Predicted relative saturation for a given suction (predicted from the GSD and porosity).
            Dg = obj.Dg_pp_from_WRC(pp);
            output = obj.FVg_Dg_from_GSD(Dg);
        end
        
        function out_FSe = Se_pp_beta_from_GSD(obj,pp,beta)
            %Predicted relative saturation for a given suction, and given the beta value (predicted from the GSD and porosity).
            Dg = obj.Dg_pp_beta_from_WRC(pp,beta);
            out_FSe = obj.FVg_Dg_from_GSD(Dg);
        end
        
        function out_pp = pp_Se_from_GSD(obj,Se)
            %Predicted suction given the relative saturation (predicted from the GSD and porosity).
            Dg = obj.Dg_FVg_from_GSD(Se);
            out_pp = obj.pp_Dg_from_GSD(Dg);
        end
        
        function out_pp = pp_Se_beta_from_GSD(obj,Se,beta)
            %Predicted suction given the relative saturation, and given the beta value (predicted from the GSD and porosity).
            Dg = obj.Dg_FVg_from_GSD(Se);
            out_pp = obj.pp_Dg_beta_from_GSD(Dg,beta);
        end
        
        function out_pp = pp_th_from_GSD(obj,th)
            %Predicted suction given the volumetric water content (predicted from the GSD and porosity).
            Se = max(0.0,min(1.0,(th-obj.thres)./(obj.thsat-obj.thres)));
            Dg = obj.Dg_FVg_from_GSD(Se);
            out_pp = obj.pp_Dg_from_GSD(Dg);
        end
        
        function out_th = th_pp_from_GSD(obj,pp)
            %Predicted volumetric water content given the suction (predicted from the GSD and porosity).
            out_th = obj.thres+(obj.thsat-obj.thres).*obj.Se_pp_from_GSD(pp);
        end
        
        function out_th = th_pp_beta_from_GSD(obj,pp,beta)
            %Predicted volumetric water content given the suction and given the beta value (predicted from the GSD and porosity).
            out_th = obj.thres+(obj.thsat-obj.thres).*obj.Se_pp_beta_from_GSD(pp,beta);
        end
        
        function out_pp = pp_th_beta_from_GSD(obj,th,beta)
            %Predicted suction given the volumetric water content and given the beta value(predicted from the GSD and porosity).
            Se = max(0.0,min(1.0,(th-obj.thres)./(obj.thsat-obj.thres)));
            Dg = obj.Dg_FVg_from_GSD(Se);
            out_pp = obj.pp_Dg_beta_from_GSD(Dg,beta);
        end
        
        function out_pp = pp_Dg_from_GSD(obj,Dg)
            %Predicted suction given the linked grain size (predicted from the GSD and porosity).
            Dp = obj.Dp_Dg_from_GSD(Dg);
            out_pp = obj.pp_Dp_from_WRC(Dp);
        end
        
        function out_pp = pp_Dg_beta_from_GSD(obj,Dg,beta)
            %Predicted suction given the linked grain size and the beta value (predicted from the GSD and porosity).
            Dp = obj.Dp_Dg_beta_from_GSD(Dg,beta);
            out_pp = obj.pp_Dp_from_WRC(Dp);
        end
        
        %% PREDICTIONS BASED ON THE MEASURED WATER RETENTION CURVE
        
        %Grain Size Distribution (GSD)
        function out_FVg = FVg_Dg_from_WRC(obj,Dg)
            %Predicted cumulated grain volume given the grain size
            %(predicted from the WRC data)
            pp = obj.pp_Dg_from_GSD(Dg);
            out_FVg = obj.Se_pp_from_WRC(pp);
        end
        
        function out_Dg = Dg_FVg_from_WRC(obj,FVg)
            %Predicted grain size given the cumulated grain volume (predicted from the WRC)
            Dp = obj.Dp_FVp_from_WRC(FVg);
            out_Dg = obj.Dg_Dp_from_WRC(Dp);
        end
        
        function out_FVg = FVg_Dg_beta_from_WRC(obj,Dg,beta)
            %Predicted cumulated grain volume given the grain size
            % and given the beta value(predicted from the WRC data)
            pp = obj.pp_Dg_beta_from_GSD(Dg,beta);
            out_FVg = obj.Se_pp_from_WRC(pp);
        end
        
        function out_Dg = Dg_FVg_beta_from_WRC(obj,FVg,beta)
            %Predicted grain size given the cumulated grain volume and given the beta value (predicted from the WRC)
            Dp = obj.Dp_FVp_from_WRC(FVg);
            out_Dg = obj.Dg_Dp_beta_from_WRC(Dp,beta);
        end
        
        function out_Dg = Dg_Dp_from_WRC(obj,Dp)
            %Predicted grain size given the pore size (predicted from the GSD data and from the porosity)
            beta = obj.beta_Dp(Dp);
            out_Dg = obj.Dg_Dp_beta_from_WRC(Dp,beta);
        end
        
        function out_Dg = Dg_Dp_beta_from_WRC(obj,Dp,beta)
            %Predicted grain size given the pore size (predicted from the GSD data and from the porosity)
            out_Dg = 2.^(1/3-5./(6.*beta)) .* 3.^(1/6.*(2+1./beta)) .* obj.evoid.^(-1./(2.*beta)) .* pi.^((-1+1./beta)./3) .* (Dp./obj.Dg0).^(1./beta) .* obj.Dg0;
        end
        
        function out_Dg = Dg_pp_from_WRC(obj,pp)
            %Predicted grain size given suction (predicted from the GSD data and from the porosity)
            Dp = obj.Dp_pp_from_WRC(pp);
            out_Dg = obj.Dg_Dp_from_WRC(Dp);
        end
        
        function out_Dg = Dg_pp_beta_from_WRC(obj,pp,beta)
            %Predicted grain size given suction and the beta value (predicted from the GSD data and from the porosity)
            Dp = obj.Dp_pp_from_WRC(pp);
            out_Dg = obj.Dg_Dp_beta_from_WRC(Dp,beta);
        end
        
        %Pore Size Distribution (PSD)
        function out_FVp = FVp_Dp_from_WRC(obj,Dp)
            %Cumulated pore volume givven the pore size, calculated from
            %the WRC.
            pp = obj.pp_Dp_from_WRC(Dp);
            out_FVp = obj.Se_pp_from_WRC(pp);
        end
        
        function out_Dp = Dp_FVp_from_WRC(obj,FVp)
            %Pore size given the cumulated pore volume, calculated from
            %the WRC.
            th = obj.thres+(obj.thsat-obj.thres).*FVp;
            pp = obj.pp_th_from_WRC(th);
            out_Dp = obj.Dp_pp_from_WRC(pp);
        end
        
        function out_Dp = Dp_pp_from_WRC(obj,pp)
            %Pore size given the suction (Youngs equation).
            out_Dp = (4.*obj.param_Ts.*obj.param_costheta)./(obj.param_rhow.*obj.param_g.*pp);
        end
        
        %Water Retention Curve (WRC)
        function output = th_pp_from_WRC(obj,pp)
            %suction given the volumetric water content from the WRC.
            derivative = false;
            isinv=false;
            output = func_soil.Interpolate_WRC(obj.WRC,pp,derivative,obj.thres,obj.thsat,obj.options.fit_to_log,isinv);
        end
        
        function out_Se = Se_pp_from_WRC(obj,pp)
            %Relative saturation given the suction from the WRC.
            th = obj.th_pp_from_WRC(pp);
            out_Se = max(0.0,min(1.0,(th-obj.thres)/(obj.thsat-obj.thres)));
        end
        
        function out_pp = pp_th_from_WRC(obj,th)
            %Suction given the volumetric water content from the WRC.
            derivative = false;
            isinv=true;
            out_pp = func_soil.Interpolate_WRC(obj.WRC,th,derivative,obj.thres,obj.thsat,obj.options.fit_to_log,isinv);
        end
        
        function out_pp = pp_Se_from_WRC(obj,Se)
            %Suction given the relative saturation from the WRC.
            th = obj.thres+(obj.thsat-obj.thres).*Se;
            out_pp = obj.pp_th_from_WRC(th);
        end
        
        function out_pp = pp_Dp_from_WRC(obj,Dp)
            %Suction given the pore diameter (Youngs Equation).
            out_pp = (4.*obj.param_Ts.*obj.param_costheta)./(obj.param_rhow.*obj.param_g.*Dp);
        end
        
        %% METHODS TO FIT BETA
        function [outbeta,RMSE,R2] = get_beta_fit_Se_from_WRC_loglinear(obj)
            %Calcultaion of beta_fix_fit by a linear regresion on each
            %point of the WRC (in the defined boundaries with option.obtSemax y optSemin).
            %(linear fit with numerator and denominator, see code)
            
            Se_i = max(0.0,min(1.0,(obj.WRC(:,2)-obj.thres)./(obj.thsat-obj.thres)));
            pp_i = obj.WRC(:,1);
            
            %Define boundaries for the fit
            mask1 = Se_i>=obj.Se_minfit;
            mask2 = Se_i<=obj.Se_maxfit;
            mask3 = pp_i>0;
            mask = mask1.*mask2.*mask3;
            %Filter data with the mask, inside the boundaries
            Se_i = Se_i(mask==1);
            pp_i = pp_i(mask==1);
            Dp_i = obj.Dp_pp_from_WRC(pp_i);
            Dg_i = obj.Dg_FVg_from_GSD(Se_i);
            
            numerator   = 1./3.*log((sqrt(3).*pi)./(sqrt(2).*4.*(obj.evoid.^(3/2))))  +log(Dp_i)  - log(obj.Dg0);
            denominator = 1./3.*log(pi./6)                                            +log(Dg_i)  - log(obj.Dg0);
            
            switch length(numerator)
                case 1 %Only 1 point
                    outbeta = numerator./denominator;
                    RMSE = 0;
                    R2 = 0;
                case 0 %No points
                    outbeta = nan;
                    RMSE = nan;
                    R2 =nan;
                otherwise %Several points
                    lmfit = fitlm(denominator,numerator,'Intercept',false);
                    outbeta = lmfit.Coefficients{1,1};
                    betai = numerator./denominator;
                    betai_fitted = lmfit.Fitted./denominator;
                    [R2,~] = func_math.Rsquared_rmse(numerator,lmfit.Fitted);
                    [~,RMSE] = func_math.Rsquared_rmse(betai,betai_fitted);
            end
            
        end
        
        %% SOME RESULTS INTO A TABLE
        function [results_table] = get_results_at_WRC_points(obj)
            %Output a table with the results of the model in the points of
            %the Water Retention Curve.
            variablenames = {'DGI' ,'BETAI','PPI','DPI','SEI','THI','SEI_BETAFIX' ,'SEI_BETAFIXTEX','SEI_BETAFIXEMP'  ,'SEI_BETAFIXFIT' ,'SEI_BETA_DP_FIT','SEI_BETA_DP_EMP','DPI_BETAFIX' ,'DPI_BETAFIXTEX' ,'DPI_BETAFIXEMP'  ,'DPI_BETAFIXFIT'  ,'DPI_BETA_DP_FIT','DPI_BETA_DP_EMP','THI_BETAFIX' ,'THI_BETAFIXTEX','THI_BETAFIXEMP'  ,'THI_BETAFIXFIT' ,'THI_BETA_DP_FIT','THI_BETA_DP_EMP','INTO_BOUND'};
            
            
            if(isempty(obj.WRC)||isempty(obj.GSD)||length(obj.WRC)<=3||length(obj.GSD)<=3)
                results_table=table;
                results_table.VariableNames = variablenames;
            else
                
                
                Se_i = max(0.0,min(1.0,(obj.WRC(:,2)'-obj.thres)./(obj.thsat-obj.thres)));
                mask_A = obj.WRC(:,1)'>0;
                
                mask1 = Se_i<=obj.Se_maxfit; %for fitting purposes
                mask2 = Se_i>=obj.Se_minfit;
                mask_B = (mask1.*mask2)==1;
                
                Se_i = Se_i(mask_A);
                mask_B = mask_B(mask_A);
                
                if (isempty(Se_i))
                    results_table=table;                   
                else
                    pp_i = obj.WRC(mask_A,1)';
                    th_i = obj.WRC(mask_A,2)';
                    dg_i = obj.Dg_FVg_from_GSD(Se_i);
                    Dp_i = obj.Dp_pp_from_WRC(pp_i);
                    
                    beta_i = obj.get_betai_Dg_Dp_from_GSDWRC(dg_i,Dp_i);
                    
                    beta_i_beta_fit_dp = obj.beta_fit_Dp(Dp_i,obj.beta_fix_fit);
                    beta_i_beta_emp_dp = obj.beta_emp_Dp(Dp_i);
                    
                    Se_beta_fix         = obj.FVp_Dp_beta_from_GSD(Dp_i,obj.beta_fix);
                    Se_beta_fix_tex     = obj.FVp_Dp_beta_from_GSD(Dp_i,obj.beta_fix_tex);
                    Se_beta_fix_emp     = obj.FVp_Dp_beta_from_GSD(Dp_i,obj.beta_fix_emp);
                    Se_beta_fix_fit     = obj.FVp_Dp_beta_from_GSD(Dp_i,obj.beta_fix_fit);
                    Se_beta_emp_dp      = obj.FVp_Dp_beta_from_GSD(Dp_i,beta_i_beta_emp_dp);
                    Se_beta_fit_dp      = obj.FVp_Dp_beta_from_GSD(Dp_i,beta_i_beta_fit_dp);
                    
                    dpi_beta_fix        = obj.Dp_Dg_beta_from_GSD(dg_i,obj.beta_fix);
                    dpi_beta_fix_tex    = obj.Dp_Dg_beta_from_GSD(dg_i,obj.beta_fix_tex);
                    dpi_beta_fix_emp    = obj.Dp_Dg_beta_from_GSD(dg_i,obj.beta_fix_emp);
                    dpi_beta_fix_fit    = obj.Dp_Dg_beta_from_GSD(dg_i,obj.beta_fix_fit);
                    
                    dpi_beta_emp_dp     = obj.Dp_Dg_beta_from_GSD(dg_i,beta_i_beta_emp_dp);
                    dpi_beta_fit_dp     = obj.Dp_Dg_beta_from_GSD(dg_i,beta_i_beta_fit_dp);
                    
                    
                    th_beta_fix         = obj.th_pp_beta_from_GSD(pp_i,obj.beta_fix);
                    th_beta_fix_tex     = obj.th_pp_beta_from_GSD(pp_i,obj.beta_fix_tex);
                    th_beta_fix_emp     = obj.th_pp_beta_from_GSD(pp_i,obj.beta_fix_emp);
                    th_beta_fix_fit     = obj.th_pp_beta_from_GSD(pp_i,obj.beta_fix_fit);
                    th_beta_emp_dp      = obj.th_pp_beta_from_GSD(pp_i,beta_i_beta_emp_dp);
                    th_beta_fit_dp      = obj.th_pp_beta_from_GSD(pp_i,beta_i_beta_fit_dp);
                    
                    
                    results_table = array2table([dg_i',beta_i',pp_i',Dp_i',Se_i',th_i',Se_beta_fix'  ,Se_beta_fix_tex' ,Se_beta_fix_emp' ,Se_beta_fix_fit' ,Se_beta_fit_dp'  ,Se_beta_emp_dp'  ,dpi_beta_fix' ,dpi_beta_fix_tex',dpi_beta_fix_emp' ,dpi_beta_fix_fit' ,dpi_beta_fit_dp' ,dpi_beta_emp_dp',th_beta_fix'  ,th_beta_fix_tex' ,th_beta_fix_emp' ,th_beta_fix_fit' ,th_beta_fit_dp'  ,th_beta_emp_dp' , mask_B'],...
                        'VariableNames',     variablenames);
                    
                end
            end
        end
        
        function out_betai = get_betai_Dg_Dp_from_GSDWRC(obj,Dg_i,Dp_i)
            %Value of beta given grain size from the GSD and pore size from
            %the PSD at the same Cummulated Volume (linked grain size and pore size).
            numerator   = 1./3.*log((sqrt(3).*pi)./(sqrt(2).*4.*(obj.evoid.^(3/2))))  +log(Dp_i)  - log(obj.Dg0);
            denominator = 1./3.*log(pi./6)                                            +log(Dg_i)  - log(obj.Dg0);
            out_betai = numerator./denominator;
        end
        
    end
end




