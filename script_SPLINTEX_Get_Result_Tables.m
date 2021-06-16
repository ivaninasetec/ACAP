%%script_SPLINTEX_Get_Result_Tables Generate output tables of SPLINTEX model.
% Data is included in the table 'DATA_INPUT' and 'SPLINTEX_OUTPUT_SIM_A'.
% Output tables are:
% DATA_FIT_SPLINTEX
% DATA_FIT_I_SPLINTEX
%
% DATA_FIT_SPLINTEX Include results for each soil data sample (so the
% number of recors match with the number of records of DATA_INPUT.
% DATA_FIT_I_SPLINTEX Include results for each point in the water retention
% curve of each soil data sample.

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
load('DATA_ACAP_SPLINTEX.mat', 'SPLINTEX_OUTPUT_SIM_A');

[DATA_FIT_SPLINTEX_A, DATA_FIT_I_SPLINTEX_A] = soil_SPLINTEX_FIT(DATA_INPUT,SPLINTEX_OUTPUT_SIM_A);
DATA_FIT_I_SPLINTEX_A.TEXTURE = categorical(DATA_FIT_I_SPLINTEX_A.TEXTURE);
%%
DATA_FIT_I_SPLINTEX_A = DATA_FIT_I_SPLINTEX_A(DATA_FIT_I_SPLINTEX_A.SEI<1,:);
DATA_FIT_I_SPLINTEX_A = DATA_FIT_I_SPLINTEX_A(DATA_FIT_I_SPLINTEX_A.SEI_VG89<1,:);
DATA_FIT_I_SPLINTEX_A = DATA_FIT_I_SPLINTEX_A(DATA_FIT_I_SPLINTEX_A.SEI>0,:);
DATA_FIT_I_SPLINTEX_A = DATA_FIT_I_SPLINTEX_A(DATA_FIT_I_SPLINTEX_A.SEI_VG89>0,:);


load('DATA_ACAP_SPLINTEX.mat', 'SPLINTEX_OUTPUT_SIM_B');

[DATA_FIT_SPLINTEX_B, DATA_FIT_I_SPLINTEX_B] = soil_SPLINTEX_FIT(DATA_INPUT,SPLINTEX_OUTPUT_SIM_B);
DATA_FIT_I_SPLINTEX_B.TEXTURE = categorical(DATA_FIT_I_SPLINTEX_B.TEXTURE);
%%
DATA_FIT_I_SPLINTEX_B = DATA_FIT_I_SPLINTEX_B(DATA_FIT_I_SPLINTEX_B.SEI<1,:);
DATA_FIT_I_SPLINTEX_B = DATA_FIT_I_SPLINTEX_B(DATA_FIT_I_SPLINTEX_B.SEI_VG89<1,:);
DATA_FIT_I_SPLINTEX_B = DATA_FIT_I_SPLINTEX_B(DATA_FIT_I_SPLINTEX_B.SEI>0,:);
DATA_FIT_I_SPLINTEX_B = DATA_FIT_I_SPLINTEX_B(DATA_FIT_I_SPLINTEX_B.SEI_VG89>0,:);


function [TABLE_OUT,TABLE_OUT_I] = soil_SPLINTEX_FIT(TABLEUNSAT,TABLE_SPLINTEX)

TABLE_OUT = array2table(TABLEUNSAT.COD,'VariableName',{'COD'}) ;
TABLE_OUT_I = table;

for i = 1:size(TABLEUNSAT,1)
    disp(['soil record: ',num2str(i)]);
    
    VG_thres = TABLE_SPLINTEX.th_res(i);
    VG_thsat = TABLE_SPLINTEX.th_sat(i);
    VG_alpha = TABLE_SPLINTEX.alpha(i);
    VG_n = TABLE_SPLINTEX.n(i);
    VG_m = TABLE_SPLINTEX.m(i);

    soil = class_soil_VanGenuchten(...
    'GSD',TABLEUNSAT.GSD{i}(:,:),...
    'WRC',TABLEUNSAT.HT_DRY{i}(:,:),...
    'npor',TABLEUNSAT.NPOR(i), ...
    'thres',VG_thres, ...
    'thsat',VG_thsat, ...
    'alpha',VG_alpha, ...
    'n',VG_n, ...
    'm',VG_m ...
    );

%FILL TABLE WITH RESULTS FOR EACH SOIL SAMPLE
TABLE_OUT.N(i) = i;
TABLE_OUT.COD(i)        = TABLEUNSAT.COD(i);
TABLE_OUT.NPOR(i)       = TABLEUNSAT.NPOR(i);
TABLE_OUT.EVOID(i)      = TABLE_OUT.NPOR(i)./(1-TABLE_OUT.NPOR(i));
TABLE_OUT.GSD(i)        = {soil.GSD};
TABLE_OUT.WRC_DRY(i)    = {soil.WRC};
 
%FILL TABLE WITH RESULTS FOR EACH POINT OF THE WATER RETENTION CURVE.   
    thtemp = soil.outtable;
    if(~isempty(thtemp))
        T_TEMP = thtemp;
        
        T_TEMP.LOGDGI(:) = log(T_TEMP.DGI);
        T_TEMP.LOGPPI(:) = log(T_TEMP.PPI);
        T_TEMP.LOGDPI(:) = log(T_TEMP.DPI);

        T_TEMP.N(:)       = TABLE_OUT.N(i);
        T_TEMP.COD(:)     = TABLE_OUT.COD(i);

        T_TEMP.NPOR(:)    = TABLE_OUT.NPOR(i);
        T_TEMP.EVOID(:)   = TABLE_OUT.NPOR(i)./(1-TABLE_OUT.NPOR(i));

        if(i==1)
            TABLE_OUT_I = T_TEMP;
        else
           TABLE_OUT_I = [TABLE_OUT_I;T_TEMP];   
        end
    end
end



end



