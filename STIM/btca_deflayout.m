function layout=btca_deflayout
% layout=btca_deflayout sets up the default layouts for augmented btc
% protocols
layout.meas_mults=[0.2 0.4 0.6 0.8 1.0]; %multipliers for in-sample (measured) directions
layout.meas_repts=1; %number of repeats for measured directions
layout.pred_mults=[1]; %multipliers for out-of-sample (prediction) directions
layout.pred_repts=2; %number of repeats for prediction directions
return
