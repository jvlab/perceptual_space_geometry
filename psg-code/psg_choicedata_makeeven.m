function [data_even,opts_used]=psg_choicedata_makeeven(data,opts)
% [data_even,opts_used]=psg_choicedata_makeeven(data,opts)
% randomly deletes trials from a choice data array (typically from a multidimensional-scaling experiment),
% so that all triads have only an even number of trials
%
% data: choice data array, typically from psg_read_choicedata
%    size is [ntriads 5], columns are [ref s1 s2 choices(d(ref,s1)<d(ref,s2)) total trials]
% opts.if_log: 1 to log, defaults to 0, can be omitted
%
% data_even:  has size [ntriads 5], but coluumn 5 is always even and
%   column 4 randomly adjusted accordingly
% opts_used: options used
% 
% See also: PSG_READ_CHOICEDATA, PSG_UMI_TRIPLIKE_DEMO, PSG_TENTLIKE_DEMO.
%
if (nargin<2)
    opts=struct;
end
opts=filldefault(opts,'if_log',1);
%
col_trials=5;
col_closer=4;
%
opts_used=opts;
%
data_even=data;
which_odd=find(mod(data(:,col_trials),2)==1);
n_odd=size(which_odd,1);
%randomly delete a "closer" response based on fraction of number of "closer" responses
frac_closer=data(which_odd,col_closer)./data(which_odd,col_trials);
which_delete=which_odd(find(rand(n_odd,1)<frac_closer));
data_even(which_delete,col_closer)=data_even(which_delete,col_closer)-1;
data_even(which_odd,col_trials)=data_even(which_odd,col_trials)-1;
if (opts.if_log)
    disp(sprintf('psg_choicedata_makeeven: %6.0f triads, %6.0f have an odd number of trials, %6.0f changed.',...
        size(data,1),sum(mod(data(:,col_trials),2)==1),sum(data(:,col_closer)~=data_even(:,col_closer))));
end
return
