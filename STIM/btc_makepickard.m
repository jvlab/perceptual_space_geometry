function [vout,minprob,maxdev]=btc_makepickard(vin,mode,dict)
% [vout,minprob,maxdev]=btc_makepickard(vin,mode,dict) completes a set of texture parameters so
% that a Pickard condition is satisfied
%
% mode='gamma_beta_card': specify gamma and a cardinal beta; other cardinal beta is assumed identical
%   vin=[gamma,beta_card]
%   vout=[beta_diag,theta]
% mode='gamma_beta_card_zero': specify gamma and a cardinal beta; other cardinal beta assumed zero
%   vin=[gamma,beta_card]
%   vout=[beta_diag,theta]
%
% minprob is minimum probability (should be >=0)
% maxdev is maximum deviation of a Pickard condition from 0
%
% See also:  BTC_LETCODE2VEC, BTC_VEC2CORRS, GETP2X2_CORRS, GETCORRS_P2X2
% 
if (nargin<3)
    dict=[];
end
if isempty(dict)
    dict=btc_define;
end
spec=btc_vec2letcode(zeros(1,length(dict.codel)),dict);
switch mode
    case 'gamma_beta_card'
        gamma=vin(1);
        beta_card=vin(2);
        qplus=(gamma+beta_card)^2/(1+gamma);
        qminus=(gamma-beta_card)^2/(1-gamma);
        beta_diag=(qplus+qminus)/2;
        theta=(qplus-qminus)/2;
        vout=[beta_diag,theta];
        spec.g=gamma(1);
        spec.b=beta_card(1);
        spec.c=beta_card(1);
        spec.e=beta_diag(1);
        spec.t=theta(1);
        spec.v=theta(1);
    case 'gamma_beta_card_zero'
        gamma=vin(1);
        beta_card=vin(2);
        qplus=gamma*(gamma+beta_card)/(1+gamma);
        qminus=gamma*(gamma-beta_card)/(1-gamma);
        beta_diag=(qplus+qminus)/2;
        theta=(qplus-qminus)/2;
        vout=[beta_diag,theta];
        spec.g=gamma(1);
        spec.b=beta_card(1);
        spec.c=0;
        spec.e=beta_diag(1);
        spec.t=theta(1);
        spec.v=theta(1);
end
% convert from spec format to vec format, and then get CIG values
vec=btc_letcode2vec(spec,dict);
corrs=btc_vec2corrs(vec,dict);
p2x2=getp2x2_corrs(corrs);
corrs=getcorrs_p2x2(p2x2,1); %suppress warnings
maxdev=max(max(abs(corrs.cig_conds(:,[1 4]))));
minprob=min(p2x2(:));
return

