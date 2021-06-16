# Matlab model and dataset for the Alternative Arya and Paris model (ACAP) to predict the water retention curve of a soil.
[![Generic badge](https://img.shields.io/badge/Build-passing-green.svg)](https://shields.io/)
[![Version](https://img.shields.io/badge/Version-v.1.0.0-blue.svg)](https://shields.io/)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/Naereen/StrapDown.js/blob/master/LICENSE)

ACAP is a model to predict the water retention curve of a soil from its grain size distribution and its porosity. The model follows the formulation included in the paper 'Campos-Guereta, I. et al, 2020' with the same title as the heading, which is an evolution of the 'Arya and Paris, 1981' model.
The model is very convenient as it uses only one parameter beta, that can be defined by a constant value (beta_fit=1.2041, for Dg0 option 1), or by more accurate values as constant values depending on soil classification (beta_fix_tex), constant beta with an empirical expression depending on characteristics of the GSD and porosity (beta_fix_emp) or more acurate empirical beta depending on the location of the water retention curve (i.e. beta_emp_dp(DP) or beta_emp_dg(Dg)).
The code is in MATLAB rev. 2007b.
This code is supplementary content for the paper 'Campos-Guereta et al., 2021', it is free to use under MIT license but the original paper need to be referred.
The paper is available at: [https://doi.org/10.1016/j.advwatres.2021.103968](https://doi.org/10.1016/j.advwatres.2021.103968)
In the paper ACAP model is compared to previous models based on the original Aria and Paris, 1981 model that constitutes an evolution of the original model. Included into this code are therefore the followings models.
Models to predict the WRC based on the Arya and Paris, 1981 model included in the repository (same MIT license):
- ACAP model (Campos-Guereta et al., 2021)
- Arya and Paris, 1981 model
- Arya and Paris, 1999 model
- Arya and Paris, 2008 model
- CNEAP model (Nimmo et al, 2007)
- Haverkamp, R. & Parlange, J.Y., 1986 model.
- Mohammadi, M.H. & Vanclooster, M. 2011 model.
- Vereecken, H. et al 1989 PTFs.
- Wösten, J.H.M. et al 2001 PTFs.
- Weynants, M. et al 2009 PTFs.
- Rosetta model (Schaap, M.G. et al 2001)


## Table of Contents
1. [A simple example](#Simple_example)
2. [Files included](#Files_included)
3. [Files included](#Files_included)
 1. [Classes](#Classes)
 2. [Input data](#Input_data)
 3. [Scripts](#Scripts)
 4. [Live scripts](#live_scripts)
 5. [To load table with data inputs](#load_data_inputs)
 6. [To load table with data inputs and results](#load_data_results)
4. [License](#license)

## A simple example <a name="Simple_example"></a>
All the code is written in **MATLAB rev. 2007b**, the core of the model is the file **'class_soil_ACAP.m'**, to build an instance of an object for a particular soil is very simple, having the grain size distribution and the porosity, an example can be checked by running the live script: **'live_script_SAMPLE_ACAP_model.mlx'**.

Mandatory inputs are the Grain Size Distribution (in particle volume) (GSD) and the porosity of the soil:
```matlab
% GRAIN SIZE DISTRIBUTION
%       Part. Diam  Aggregated Particle Volume passing
%       [m]         [m3·m-3]
GSD = [ 0           ,0;...
        2.00E-06    ,0.020;...
        5.00E-05    ,0.130;...
        0.000106    ,0.256;...
        0.000250    ,0.743;...
        0.000500    ,0.916;...
        0.001000    ,0.985;...
        0.002000    ,0.997];
% POROSITY
npor = 0.45;
```
The default options of the model can also be modified:
```matlab
    optionsACAP = struct(...
	'fit_to_log',true,... % If true: Interpolations are done in the logDg instead on Dg
	'saturation_vs_porosity',1.0,... % The saturated water content will be considered as porosity multiplied by this coefficient (1.0 for thsat=npor)
	'opt_beta',5,... %0: beta_set 1:fix, 2:tex, 3:fix_fit,4:fix_emp, 5:emp_dp, 6:fit_dp
	'opt_Dg0',1);
```
To build an instance of the model:
```matlab
soil = class_soil_ACAP(...
    'GSD',GSD,...
    'npor',npor,...
    'options',optionsACAP);
```
With this, an object 'soil' is created, with several properties (as the soil texture, characteristic diameters,...) and with several methods to calculate properties of the Grain Size Distribution Curve, Pore Size distribution and Water Retention Curve of the Soil.

**Examples of possibles outputs of the model:**
```matlab
%Grain Size in the GSD at 20% (D20) in [m]:
soil.Dg_FVg_from_GSD(0.2)

%Pore Size for the quantile 40% [m] (for the aggregated pore volume of 40% in the PSD) predicte from the GSD:
soil.Dp_FDp_from_GSD(0.4)

%Pore Size corresponding to the Particle diameter of 1E-4 [m] predicted with the model:
soil.Dp_Dg_from_GSD(1E-4)

%Pore Volume for the Pore size of 1E-4[m], on the PSD predicted with the model [per unit]:
soil.FVp_Dp_from_GSD(1E-4)

%Volumetric Water content for a suction of 1E-7[m] predicted by the model:
soil.th_pp_from_GSD(1E-7)

%Suction [m] for a volumetric water content of 0.35 m3/m3:
soil.pp_th_from_GSD(0.35)

%Suction [m] prediction for a relative saturation of 60%:
soil.pp_Se_from_GSD(0.60)
```
## Files included<a name="Files_included"></a>
The repository include classes to build the ACAP model and the Arya and Paris based models. Also include Matlab 'live scripts' with examples of the models and also with the results of the empirical fittings for the model parameter 'beta'. And also some scripts to get tables with the results from input data.

### Classes<a name="Classes"></a>
The following classes include integrates each of the models. Also the class 'class_soil_PLOT' can be used to plot the results (i.e. pore size distribution) for any of the models of for some of them.
| Class | Description |
| ------ | ------ |
| [class_soil_ACAP.m](class_soil_ACAP.m) | Class with the ACAP model |
| [class_soil_AP81.m](class_soil_AP81.m) | Class with the Arya and Paris, 1981 model (AP81) |
| [class_soil_AP99.m](class_soil_AP99.m) | Class with the Arya and Paris, 1999 model (AP99) |
| [class_soil_AP08.m](class_soil_AP08.m) | Class with the Arya and Paris, 2008 model (AP08) |
| [class_soil_CNEAP.m](class_soil_CNEAP.m) | Class with the CNEAP model (CNEAP)|
| [class_soil_Haverkamp.m](class_soil_Haverkamp.m) | Class with the Haverkamp, R. & Parlange, J.Y., 1986 model |
| [class_soil_Mohammadi.m](class_soil_Mohammadi.m) | Class with the Mohammadi, M.H. & Vanclooster, M. 2011 model |
| [class_soil_Vereckeen.m](class_soil_Vereckeen.m) | Class with the Vereecken, H. et al 1989 PTFs |
| [class_soil_Wosten.m](class_soil_Wosten.m) | Class with the Wösten, J.H.M. et al 2001 PTFs |
| [class_soil_Weynants.m](class_soil_Weynants.m) | Class with the Weynants, M. et al 2009 PTFs |
| [class_soil_VanGenuchten.m](class_soil_VanGenuchten.m) | Class to implement any soil known Van Genucthen parameters |
| [class_soil_PLOT.m](class_soil_PLOT.m) | Class to plot the curves of the Pore Size Distribution on any of previous models |

### Input data<a name="Input_data"></a>
The input data is included in the file **DATA_ACAP.mat** under the name of **DATA_INPUT** include records for soils included in the UNSODA database (Leij et al, 1996), from Zapata et al (2000), and from Chiapponi (2017).
The table **ROSETTA_SCCBD** includes the Van Genuchten parameters obtained by running Rosetta model on this dataset. 
From these input data, and running each of the models with the scripts described in the following section, then tables with the results of the models are obtained. This table results are also included in the file **DATA_ACAP.mat**.
### Scripts<a name="Scripts"></a>
The following scripts need to be run in order to get the table of results (needed also to run the live scripts with the fittings for the beta parameter).
| Script | Description |
| ------ | ------ |
| [script_ACAP_Get_Results_Tables.m](script_ACAP_Get_Results_Tables.m) | Get results with the ACAP model from the data included in DATA_INPUT |
| [script_AP81_Get_Result_Tables.m](script_AP81_Get_Result_Tables.m) | Get results with the AP81 model from the data included in DATA_INPUT |
| [script_AP99_Get_Result_Tables.m](script_AP99_Get_Result_Tables.m) | Get results with the AP99 model from the data included in DATA_INPUT |
| [script_AP08_Get_Result_Tables.m](script_AP08_Get_Result_Tables.m) | Get results with the AP08 model from the data included in DATA_INPUT|
| [script_CNEAP_Get_Result_Tables.m](script_CNEAP_Get_Result_Tables.m) | Get results with the CNEAP model from the data included in DATA_INPUT|
| [script_Haverkamp_Get_Result_Tables.m](script_Haverkamp_Get_Result_Tables.m) | Get results with the Haverkamp, R. & Parlange, J.Y., 1986 model from the data included in DATA_INPUT|
| [script_Mohammadi_Get_Result_Tables.m](script_Mohammadi_Get_Result_Tables.m) | Get results with the Mohammadi, M.H. & Vanclooster, M. 2011 model from the data included in DATA_INPUT|
| [script_Vereecken_Get_Result_Tables.m](script_Vereecken_Get_Result_Tables.m) | Get results with the Wösten, J.H.M. et al 2001 PTFs from the data included in DATA_INPUT|
| [script_Wosten_Get_Result_Tables.m](script_Wosten_Get_Result_Tables.m) | Get results with the Wösten, J.H.M. et al 2001 PTFs from the data included in DATA_INPUT|
| [script_Weynants_Get_Result_Tables.m](script_Weynants_Get_Result_Tables.m) | Get results with the Weynants, M. et al 2009 PTFs from the data included in DATA_INPUT|
| [script_Rosetta_Get_Result_Tables.m](script_Rosetta_Get_Result_Tables.m) | Get results with the Rosetta model with the Van-Genuchten parameters in table **ROSETTA_SCCBD**. This table has been previously obtained by running the Rosetta model on the dataset included in **DATA_ACAP.mat**.|

### Live Scripts<a name="live_scripts"></a>
Some of the included live scripts perform samples for each of the models:
| Live Script | Description |
| ------ | ------ |
| [live_script_SAMPLE_ACAP_model.m](live_script_SAMPLE_ACAP_model.m) | Sample of application of the ACAP model.|
| [live_script_SAMPLE_AP81_model.m](live_script_SAMPLE_AP81_model.m) | Sample of application of the AP81 model. |
| [live_script_SAMPLE_AP99_model.m](live_script_SAMPLE_AP99_model.m) | Sample of application of the AP99 model. |
| [live_script_SAMPLE_AP08_model.m](live_script_SAMPLE_AP08_model.m) | Sample of application of the AP08 model.|
| [live_script_SAMPLE_CNEAP_model.m](live_script_SAMPLE_CNEAP_model.m) | Sample of application of the CNEAP model.|
| [live_script_SAMPLE_Haverkamp_model.m](live_script_SAMPLE_Haverkamp_model.m) | Sample of application of the Haverkamp, R. & Parlange, J.Y., 1986 model.|
| [live_script_SAMPLE_Mohammadi_model.m](live_script_SAMPLE_Mohammadi_model.m) | Sample of application of the Mohammadi, M.H. & Vanclooster, M. 2011 model.|
| [live_script_SAMPLE_Vereecken_model.m](live_script_SAMPLE_Vereecken_model.m) | Sample of application of the Vereecken, H. et al 1989 PTFs.|
| [live_script_SAMPLE_Wosten_model.m](live_script_SAMPLE_Wosten_model.m) | Sample of application of the Wösten, J.H.M. et al 2001 PTFs.|
| [live_script_SAMPLE_Weynants_model.m](live_script_SAMPLE_Weynants_model.m) | Sample of application of the Weynants, M. et al 2009 PTFs.|
| [live_script_SAMPLE_VanGenuchten_model.m](live_script_SAMPLE_VanGenuchten_model.m) | Sample of application of a soil knowns Van Genuchten parameters.|

Other livescripts are written to generate all the fittings fot the **beta** parameters from the output tables generated with the **scripts**.
| Live Script | Description |
| ------ | ------ |
| [live_script_FIT_beta_fix_emp_all_points.m](live_script_FIT_beta_fix_emp_all_points.m) | Fit for the **beta_fix_emp** parameter from all points of measured WRCs for the soils (fixed beta defined by an empirical expression). |
| [live_script_FIT_beta_fix_emp_all_samples.m](live_script_FIT_beta_fix_emp_all_samples.m) | Fit for the **beta_fix_emp** parameter from all **beta_fit** for each soil sample (fixed beta defined by an empirical expression). |
| [live_script_FIT_beta_fix_emp_Dp.m](live_script_FIT_beta_fix_emp_Dp.m) |Fit for the **beta_emp_dp** parameter from all points of measured WRCs for the soils (beta depending on pore size (or grain size or quantile fraction) defined by an empirical expression).|
| [live_script_FIT_beta_fix_fit_Dp.m](live_script_FIT_beta_fix_fit_Dp.m) | Fit for the **beta_fit_dp** parameter from all **beta_fit** for each soil sample and the points of the WRC (beta depending on pore size (or grain size or quantile fraction) defined by an empirical expression).|

*Note*: **beta_fix_fit** is the fixed beta parameter fitted to the WRC with the process defined in Campos-Guereta, I. et al, 2020. This value is calculated when the model instance is constructed with the WRC as an input.

Other livescripts are graphics included in the paper Campos-Guereta, I. et al, 2020:
| Live Script | Description |
| ------ | ------ |
| [live_script_PLOT_beta_tex.m](live_script_PLOT_beta_tex.m) | Graphics with the values of **beta_fix_tex** for each soil texture and depending on the option for DG0 selected.|
| [live_script_PLOT_beta_fit.m](live_script_PLOT_beta_fit.m) | Graphics with the sample on how the model fit the **beta_fix_fit** value if the WRC of the soil is given as an input.|

Finally a live script plot all the Pore Size Distributions predicted versus the Pore Size Distribution measured on each soil. (this could be done similarly and easily with the Water Retention Curves). The graphics shows the goodness of the fit, depending on the option of **beta** selected, and the option of **Dg0**. Also include the predicted vs measured PSD for older Arya and Paris models for comparison purposes.
| Live Script | Description |
| ------ | ------ |
| [live_script_PLOT_beta_tex.m](live_script_PLOT_beta_tex.m) | Graphics with the values of **beta_fix_tex** for each soil texture and depending on the option for DG0 selected.|

### To load table with data inputs<a name="load_data_inputs"></a>
``` matlab
load('DATA_ACAP.mat', 'DATA_INPUT', 'ROSETTA_SCCBD');
```
### To load table with data inputs and results<a name="load_data_results"></a>
``` matlab
load('DATA_ACAP.mat');
```

License<a name="license"></a>
----
MIT
Conditions:
- The license copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
- Paper must be quoted whith the following data when refering to this code or when the code is applied: '*CAMPOS-GUERETA, I., DAWSON, A. & THOM, N. 2021. An alternative continuous form of Arya and Paris model to predict the soil water retention curve of a soil. Advances in Water Resources, 154, 103968*'

# References
- CAMPOS-GUERETA, I., DAWSON, A. & THOM, N. 2021. An alternative continuous form of Arya and Paris model to predict the soil water retention curve of a soil. Advances in Water Resources, 154, 103968
- ARYA, L. M. & PARIS, J. F. 1981. A Physicoempirical Model to Predict the Soil Moisture Characteristic from Particle-Size Distribution and Bulk Density Data. Soil Science Society of America Journal, 45, 1023-1030.
- ARYA, L. M., LEIJ, F. J., VAN GENUCHTEN, M. T. & SHOUSE, P. J. 1999. Scaling Parameter to Predict the Soil Water Characteristic from Particle-Size Distribution Data. Soil Science Society of America Journal, 63, 510-519.
- ARYA, L. M., BOWMAN, D. C., THAPA, B. B. & CASSEL, D. K. 2008. Scaling Soil Water Characteristics of Golf Course and Athletic Field Sands from Particle-Size Distribution. Soil Science Society of America Journal, 72, 25-32.
- NIMMO, J. R., HERKELRATH, W. N. & LAGUNA LUNA, A. M. 2007. Physically Based Estimation of Soil Water Retention from Textural Data: General Framework, New Models, and Streamlined Existing Models. Vadose Zone Journal, 6, 766-773.
- LEIJ, F. J., ALVES, W. J., VAN GENUCHTEN, M. T. & WILLIAMS, J. R. 1996. The UNSODA Unsaturated Soil Hydraulic Database: User's Manual Version 1.0. National Risk Management Research Laboratory.
- ZAPATA, C. E., HOUSTON, W. N., HOUSTON, S. L. & WALSH, K. D. 2000. Soil-Water Characteristic Curve Variability. Advances in Unsaturated Geotechnics, 84-124.
- CHIAPPONI, L. 2017. Water retention curves of multicomponent mixtures of spherical particles. Powder Technology, 320, 646-655.
- HAVERKAMP, R. & PARLANGE, J. Y. 1986. Predicting the water-retention curve from particle-size distribution: 1. Sandy soils without organic matter. Soil Science, 142, 325-339.
- MOHAMMADI, M. H. & VANCLOOSTER, M. 2011. Predicting the Soil Moisture Characteristic Curve from Particle Size Distribution with a Simple Conceptual Model. Vadose Zone Journal, 10, 594-602.
- VEREECKEN, H., MAES, J., FEYEN, J. & DARIUS, P. 1989. Estimating the soil moisture retention characteristic from texture, bulk density, and carbon content. Soil science, 148, 389-403.
- WÖSTEN, J. H. M., PACHEPSKY, Y. A. & RAWLS, W. J. 2001. Pedotransfer functions: bridging the gap between available basic soil data and missing soil hydraulic characteristics. Journal of Hydrology, 251, 123-150.
- WÖSTEN, J. H. M., PACHEPSKY, Y. A. & RAWLS, W. J. 2001. Pedotransfer functions: bridging the gap between available basic soil data and missing soil hydraulic characteristics. Journal of Hydrology, 251, 123-150.
- WEYNANTS, M., VEREECKEN, H. & JAVAUX, M. 2009. Revisiting Vereecken Pedotransfer Functions: Introducing a Closed-Form Hydraulic Model. Vadose Zone Journal, 8, 86-95.
- SCHAAP, M. G., LEIJ, F. J. & VAN GENUCHTEN, M. T. 2001. rosetta : a computer program for estimating soil hydraulic parameters with hierarchical pedotransfer functions. Journal of Hydrology, 251, 163-176.
