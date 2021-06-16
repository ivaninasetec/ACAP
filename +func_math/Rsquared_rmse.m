function [out_R2,out_RMSE, out_ME] = Rsquared_rmse(MEASURED,PREDICTED)
mask1=not(isnan(MEASURED));
mask2=not(isnan(PREDICTED));
% mask3=isreal(MEASURED);
% mask4=isreal(PREDICTED);
mask = mask1.*mask2;

MEASURED = real(MEASURED(mask==1));
PREDICTED = real(PREDICTED(mask==1));

RESIDUAL = (MEASURED-PREDICTED);
NLENGTH = length(RESIDUAL);

MEAN_OBSERVED = mean(MEASURED);
SSTOT = sum((MEASURED-MEAN_OBSERVED).^2);
%SSREG = sum((PREDICTED-mean(MEASURED)).^2);
SSRES = sum(RESIDUAL.^2);

out_R2 = 1- SSRES./SSTOT;
out_RMSE  = sqrt(SSRES/NLENGTH);
out_ME = sum(abs(RESIDUAL))/NLENGTH;
end

