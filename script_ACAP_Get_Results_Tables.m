%%script_ACAP_Get_Results_Tables Generate output tables of ACAP model for each option of Dg0
% Data is included in the table 'DATA_INPUT'.
% Output tables are:
% DATA_FIT_DG1_ACAP
% DATA_FIT_I_DG1_ACAP
% DATA_FIT_DG2_ACAP
% DATA_FIT_I_DG2_ACAP
% DATA_FIT_DG3_ACAP
% DATA_FIT_I_DG3_ACAP
% DATA_FIT_DG4_ACAP
% DATA_FIT_I_DG4_ACAP
% DATA_FIT_DG5_ACAP
% DATA_FIT_I_DG5_ACAP
%
% DATA_FIT_DG*_ACAP Include results for each soil data sample (so the
% number of recors match with the number of records of DATA_INPUT.
% DATA_FIT_I_DG*_ACAP Include results for each point in the water retention
% curve of each soil data sample.
% These result tables let us to do fittings and empirical aproximations of
% the ACAP model parameter beta. (live_script's included).

%% INPUTS:
% Table with input data (see example included in DATA_ACAP). The fields
% are:
% N: Number of the record
% COD: Code of the soil sample
% NPOR: Porosity
% GSD: list with an array with the Grain Size Distribution
% HT_DRY: Water Retention Curve (drying curve), needed for fitting beta to
% real values, get empirical expressions...

% Load DATA_INPUT:
load('DATA_ACAP.mat', 'DATA_INPUT');

DATATABLE = DATA_INPUT;

% Select different opt_fitminSe and opt_fitmaxSe to filter the points of
% the WRC to fit beta values (it is better to avoid saturated part of the
% WRC and the residual part of the WRC) for that reason these boundaries
% are selected.
opt_fitmaxSe = 0.90; %Select only data below this percentage on the WRC to avoid the horizontal saturated part of the WRC.
opt_fitminSe = 0.20; %Select only data above this percentage on the WRC to avoid the residual part of the curve.

% This parameter don't need to be changed as results for each beta arr
% given.
optbeta = 1; %1:Fix, 2:fix_tex, 3:fit, 4:fix_emp, 5:emp_Dg(Dp), 6:fit_Dg(Dp)

for i=1:5
%Names of the output tables:
DATATABLEOUT_str    = strcat('DATA_FIT_DG',int2str(i),'_ACAP');
DATATABLEOUT_I_str  = strcat('DATA_FIT_I_DG',int2str(i),'_ACAP');

optDg0 = i  ; 

%Calculate results for all records.
[DATATABLEOUT,DATATABLEOUT_I] = RESULTS_IN_ALL_RECORDS(DATATABLE,optbeta,optDg0,opt_fitmaxSe,opt_fitminSe);

%Change name of the tables:
eval([DATATABLEOUT_str,'=DATATABLEOUT;']);
eval([DATATABLEOUT_I_str,'=DATATABLEOUT_I;']);
end
clear('DATATABLE','DATATABLEOUT','DATATABLE_95','DATATABLEOUT_95','DATATABLEOUT_I','DATATABLEOUT_I_95','DATATABLEOUT_str','DATATABLEOUT_I_str','optbeta','optVg0','optDg0','i','opt_fitmaxSe','opt_fitminSe');



function [DATATABLEOUT,DATATABLEOUT_I] = RESULTS_IN_ALL_RECORDS(TABLEUNSAT,optbeta,optDg0,opt_fitmaxSe,opt_fitminSe)
    options = struct(...
	'fit_to_log',true,...
	'fit_max_Se',opt_fitmaxSe,...
	'fit_min_Se',opt_fitminSe,...
	'saturation_vs_porosity',1.0,...
	'opt_beta',optbeta,...
	'opt_Dg0',optDg0);
DATATABLEOUT = array2table(TABLEUNSAT.COD,'VariableName',{'COD'}) ;
DATATABLEOUT_I = table;

% Loop on each record (each soil sample)
for i = 1:size(TABLEUNSAT,1)
    disp(['Dg0 option: ',num2str(optDg0),' Soil record: ',num2str(i)]);
    
    %Construct an instance of the model for record i
    soil = class_soil_ACAP(...
	'GSD',TABLEUNSAT.GSD{i}(:,:),...
	'WRC',TABLEUNSAT.HT_DRY{i}(:,:),...
	'npor',TABLEUNSAT.NPOR(i),...
	'options',options);

%BEGIN TO FILL THE TABLE 
DATATABLEOUT.N(i)       = i;
DATATABLEOUT.COD(i)     = TABLEUNSAT.COD(i);
DATATABLEOUT.NPOR(i)    = TABLEUNSAT.NPOR(i);
DATATABLEOUT.EVOID(i)   = DATATABLEOUT.NPOR(i)./(1-DATATABLEOUT.NPOR(i));
DATATABLEOUT.GSD(i)     = {soil.GSD};
DATATABLEOUT.WRC_DRY(i) = {soil.WRC};
%List with a table with the results for this soil:
DATATABLEOUT.RESULTS_AT_WRC_PT(i) = {soil.results_at_WRC_points};

DATATABLEOUT.DG_MIN(i) = soil.Dg_min;
DATATABLEOUT.D10(i) = soil.D10;
DATATABLEOUT.D20(i) = soil.D20;
DATATABLEOUT.D30(i) = soil.D30;
DATATABLEOUT.D40(i) = soil.D40;
DATATABLEOUT.D50(i) = soil.D50;
DATATABLEOUT.D60(i) = soil.D60;
DATATABLEOUT.D80(i) = soil.D80;
DATATABLEOUT.D80(i) = soil.D80;
DATATABLEOUT.DMAX(i) = soil.Dg_max;
DATATABLEOUT.CC(i) = soil.cc;
DATATABLEOUT.CU(i) = soil.cu;
DATATABLEOUT.LOGDG_MIN(i) = log(soil.Dg_min);
DATATABLEOUT.LOGD10(i) = log(soil.D10);
DATATABLEOUT.LOGD20(i) = log(soil.D20);
DATATABLEOUT.LOGD30(i) = log(soil.D30);
DATATABLEOUT.LOGD40(i) = log(soil.D40);
DATATABLEOUT.LOGD50(i) = log(soil.D50);
DATATABLEOUT.LOGD60(i) = log(soil.D60);
DATATABLEOUT.LOGD80(i) = log(soil.D80);
DATATABLEOUT.LOGDMAX(i) = log(soil.Dg_max);
DATATABLEOUT.LOGCC(i) = log(soil.cc);
DATATABLEOUT.LOGCU(i) = log(soil.cu);

DATATABLEOUT.DP10(i) = soil.DP10;
DATATABLEOUT.DP20(i) = soil.DP20;
DATATABLEOUT.DP30(i) = soil.DP30;
DATATABLEOUT.DP40(i) = soil.DP40;
DATATABLEOUT.DP50(i) = soil.DP50;
DATATABLEOUT.DP60(i) = soil.DP60;
DATATABLEOUT.DP80(i) = soil.DP80;

DATATABLEOUT.LOGDP10(i) = log(soil.DP10);
DATATABLEOUT.LOGDP20(i) = log(soil.DP20);
DATATABLEOUT.LOGDP30(i) = log(soil.DP30);
DATATABLEOUT.LOGDP40(i) = log(soil.DP40);
DATATABLEOUT.LOGDP50(i) = log(soil.DP50);
DATATABLEOUT.LOGDP60(i) = log(soil.DP60);
DATATABLEOUT.LOGDP80(i) = log(soil.DP80);

%%Beta values
DATATABLEOUT.BETA_FIX(i)           = soil.beta_fix;
DATATABLEOUT.BETA_FIX_TEX(i)       = soil.beta_fix_tex;
DATATABLEOUT.BETA_FIX_FIT(i)       = soil.beta_fix_fit;
DATATABLEOUT.BETA_FIX_FIT_RMSE(i)  = soil.beta_fix_fit_RMSE;
DATATABLEOUT.BETA_FIX_FIT_R2(i)    = soil.beta_fix_fit_R2;
DATATABLEOUT.BETA_FIX_EMP(i)       = soil.beta_fix_emp;

%thres and thsat
DATATABLEOUT.THRES(i) = soil.thres;
DATATABLEOUT.THSAT(i) = soil.thsat;

%Soil classification.
DATATABLEOUT.CLAY(i) = soil.clay;
DATATABLEOUT.SILT(i)  = soil.silt;
DATATABLEOUT.SAND(i)  = soil.sand;
DATATABLEOUT.GRAVEL(i) = soil.gravel;
DATATABLEOUT.TEXTURE(i) = {soil.texture};

%Dg0 parameter:
DATATABLEOUT.DG0(i) = soil.Dg0;


%BEGIN TO FILL TABLE WITH RESULTS FOR EACH WATER RETENTION CURVE POINT.    
    thtemp = DATATABLEOUT.RESULTS_AT_WRC_PT(i); 
    thtemp = thtemp{1};
    if(~(isempty(thtemp)))
        T_TEMP = thtemp;
        
        T_TEMP.N(:)         = DATATABLEOUT.N(i);
        T_TEMP.COD(:)       = DATATABLEOUT.COD(i);
        
        T_TEMP.LOGDGI(:) = log(T_TEMP.DGI);
        T_TEMP.LOGPPI(:) = log(T_TEMP.PPI);
        T_TEMP.LOGDPI(:) = log(T_TEMP.DPI);

        T_TEMP.CC(:)        = DATATABLEOUT.CC(i);
        T_TEMP.CU(:)        = DATATABLEOUT.CU(i);
        T_TEMP.D10(:)       = DATATABLEOUT.D10(i);
        T_TEMP.D20(:)       = DATATABLEOUT.D20(i);
        T_TEMP.D30(:)       = DATATABLEOUT.D30(i);
        T_TEMP.D40(:)       = DATATABLEOUT.D40(i);
        T_TEMP.D50(:)       = DATATABLEOUT.D50(i);
        T_TEMP.D60(:)       = DATATABLEOUT.D60(i);
        T_TEMP.D80(:)       = DATATABLEOUT.D80(i);
        T_TEMP.DMAX(:)      = DATATABLEOUT.DMAX(i);
        T_TEMP.LOGD10(:)    = log(DATATABLEOUT.D10(i));
        T_TEMP.LOGD20(:)    = log(DATATABLEOUT.D20(i));
        T_TEMP.LOGD30(:)    = log(DATATABLEOUT.D30(i));
        T_TEMP.LOGD40(:)    = log(DATATABLEOUT.D40(i));
        T_TEMP.LOGD50(:)    = log(DATATABLEOUT.D50(i));
        T_TEMP.LOGD60(:)    = log(DATATABLEOUT.D60(i));
        T_TEMP.LOGD80(:)    = log(DATATABLEOUT.D80(i));
        T_TEMP.LOGDMAX(:)   = log(DATATABLEOUT.DMAX(i));
        T_TEMP.LOGCU(:)   = log(DATATABLEOUT.CU(i));
        T_TEMP.LOGCC(:)   = log(DATATABLEOUT.CC(i));
        
        T_TEMP.DP10(:)       = DATATABLEOUT.DP10(i);
        T_TEMP.DP20(:)       = DATATABLEOUT.DP20(i);
        T_TEMP.DP30(:)       = DATATABLEOUT.DP30(i);
        T_TEMP.DP40(:)       = DATATABLEOUT.DP40(i);
        T_TEMP.DP50(:)       = DATATABLEOUT.DP50(i);
        T_TEMP.DP60(:)       = DATATABLEOUT.DP60(i);
        T_TEMP.DP80(:)       = DATATABLEOUT.DP80(i);
        
        T_TEMP.LOGDP10(:)    = log(DATATABLEOUT.DP10(i));
        T_TEMP.LOGDP20(:)    = log(DATATABLEOUT.DP20(i));
        T_TEMP.LOGDP30(:)    = log(DATATABLEOUT.DP30(i));
        T_TEMP.LOGDP40(:)    = log(DATATABLEOUT.DP40(i));
        T_TEMP.LOGDP50(:)    = log(DATATABLEOUT.DP50(i));
        T_TEMP.LOGDP60(:)    = log(DATATABLEOUT.DP60(i));
        T_TEMP.LOGDP80(:)    = log(DATATABLEOUT.DP80(i));
        
        T_TEMP.THRES(:)     = DATATABLEOUT.THRES(i);
        T_TEMP.THSAT(:)     = DATATABLEOUT.THSAT(i);
               
        T_TEMP.BETA_FIX(:)          = DATATABLEOUT.BETA_FIX(i);
        T_TEMP.BETA_FIX_TEX(:)      = DATATABLEOUT.BETA_FIX_TEX(i);
        T_TEMP.BETA_FIX_FIT(:)      = DATATABLEOUT.BETA_FIX_FIT(i);
        T_TEMP.BETA_FIX_FIT_RMSE(:)= DATATABLEOUT.BETA_FIX_FIT_RMSE(i);
        T_TEMP.BETA_FIX_FIT_R2(:)= DATATABLEOUT.BETA_FIX_FIT_R2(i);
        T_TEMP.BETA_FIX_EMP(:)     = DATATABLEOUT.BETA_FIX_EMP(i);

        T_TEMP.CLAY(:)      = DATATABLEOUT.CLAY(i);
        T_TEMP.SILT(:)      = DATATABLEOUT.SILT(i);
        T_TEMP.SAND(:)      = DATATABLEOUT.SAND(i);
        T_TEMP.GRAVEL(:)    = DATATABLEOUT.GRAVEL(i);
        T_TEMP.NPOR(:)      = DATATABLEOUT.NPOR(i);
        T_TEMP.EVOID(:)      = DATATABLEOUT.NPOR(i)./(1-DATATABLEOUT.NPOR(i));
        
        T_TEMP.TEXTURE(:)   = DATATABLEOUT.TEXTURE(i);
        T_TEMP.DG0(:) = DATATABLEOUT.DG0(i);

        %Append data to final table_i
        if(i==1)
            DATATABLEOUT_I = T_TEMP;
        else
           DATATABLEOUT_I = [DATATABLEOUT_I;T_TEMP];   
        end
    end
end

DATATABLEOUT.TEXTURE(isnan(DATATABLEOUT.D50),:)={''};
DATATABLEOUT_I.TEXTURE(isnan(DATATABLEOUT_I.D50),:)={''};
DATATABLEOUT.TEXTURE = categorical(DATATABLEOUT.TEXTURE);
DATATABLEOUT_I.TEXTURE = categorical(DATATABLEOUT_I.TEXTURE);

DATATABLEOUT_I=reclasify_textures(DATATABLEOUT_I);
DATATABLEOUT=reclasify_textures(DATATABLEOUT);
DATATABLEOUT.TEXTURE2 = categorical(DATATABLEOUT.TEXTURE2);
DATATABLEOUT_I.TEXTURE2 = categorical(DATATABLEOUT_I.TEXTURE2);
end


function obj = reclasify_textures(F0)
%Create a new TEXTURE2 field to group some textures with few data.
F0.TEXTURE2=F0.TEXTURE;
F0.TEXTURE2(F0.TEXTURE=='clay')='C+SIC';
F0.TEXTURE2(F0.TEXTURE=='silty clay')='C+SIC';
F0.TEXTURE2(F0.TEXTURE=='silt clay loam')='SICL+CL';
F0.TEXTURE2(F0.TEXTURE=='clay loam')='SICL+CL';
F0.TEXTURE2(F0.TEXTURE=='silt loam')='SIL+SI';
F0.TEXTURE2(F0.TEXTURE=='silt')='SIL+SI';
F0.TEXTURE2(F0.TEXTURE=='sandy clay')='SC+SCL';
F0.TEXTURE2(F0.TEXTURE=='sandy clay loam')='SC+SCL';
F0.TEXTURE2(F0.TEXTURE=='loam')='L';
F0.TEXTURE2(F0.TEXTURE=='loamy sand')='LS';
F0.TEXTURE2(F0.TEXTURE=='sand')='S';
F0.TEXTURE2(F0.TEXTURE=='sandy loam')='SL';
F0.TEXTURE2=removecats(F0.TEXTURE2);
obj = F0;
end
