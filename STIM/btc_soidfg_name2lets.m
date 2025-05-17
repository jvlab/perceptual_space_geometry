function lets=btc_soidfg_name2lets(param_name)
% lets=btc_soidfg_name2lets(param_name) is a utility function that extracts
% the coordinate letters (one or two of gbcdetuvwa) from a parameter name,
% such as fig_g, or fmg_b, or fmg_b_c, or fig_gnd_b_c, etc., as generated
% by btc_soidfg_model.  
%
% logic is to exclude known strings (fig,gnd,fmg), and then parse by
% underscores to find singletons
%
%   See also:  BTC_SOIDFG_MODEL.
%
tags={'fig','gnd','fmg'};
for itag=1:length(tags)
    param_name=strrep(param_name,tags{itag},'');
end
remain=param_name;
lets=[];
while ~isempty(remain)
    [let,remain]=strtok(remain,'_');
    if length(let)==1
        lets=[lets,let];
    end
end
return

