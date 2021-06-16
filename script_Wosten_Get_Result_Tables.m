%%script_Wosten_Get_Result_Tables Generate output tables of 
% Wösten, J.H.M. et al 2001 pedotransfer functions.
% Data is included in the table 'DATA_INPUT'.
% Output tables are:
% DATA_FIT_Wosten
% DATA_FIT_I_Wosten
%
% DATA_FIT_Wosten Include results for each soil data sample (so the
% number of recors match with the number of records of DATA_INPUT.
% DATA_FIT_I_Wosten Include results for each point in the water retention
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

[DATA_FIT_Wosten, DATA_FIT_I_Wosten] = soil_Wosten_FIT(DATA_INPUT,2);
DATA_FIT_I_Wosten.TEXTURE = categorical(DATA_FIT_I_Wosten.TEXTURE);
%%
DATA_FIT_I_Wosten = DATA_FIT_I_Wosten(DATA_FIT_I_Wosten.SEI<1,:);
DATA_FIT_I_Wosten = DATA_FIT_I_Wosten(DATA_FIT_I_Wosten.SEI_WO99<1,:);
DATA_FIT_I_Wosten = DATA_FIT_I_Wosten(DATA_FIT_I_Wosten.SEI>0,:);
DATA_FIT_I_Wosten = DATA_FIT_I_Wosten(DATA_FIT_I_Wosten.SEI_WO99>0,:);

function [TABLE_OUT,TABLE_OUT_I] = soil_Wosten_FIT(TABLEUNSAT,model)

TABLE_OUT = array2table(TABLEUNSAT.COD,'VariableName',{'COD'}) ;
TABLE_OUT_I = table;

for i = 1:size(TABLEUNSAT,1)
    disp(['soil record: ',num2str(i)]);

    soil = class_soil_Wosten(...
    'GSD',TABLEUNSAT.GSD{i}(:,:),...
    'WRC',TABLEUNSAT.HT_DRY{i}(:,:),...
    'npor',TABLEUNSAT.NPOR(i));

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



