% psg_tent_stats_demo: demonstrate psg_tent_stats
%
% 20Feb23: add psg_data_dir
%
%   See also:  PSG_UMI_STATS, PSG_TENT_STATS, PSG_SPOKES_SETUP.
%
if ~exist('sessfiles') 
    sessfiles={'bc6pt9.mat','bcpm3pt9.mat','bdce3pt9.mat','bgca3pt9.mat','tvpm3pt9.mat'};
end
if ~exist('psg_data_dir')
    psg_data_dir='.\psg_data\';
end
if ~exist('if_log') if_log=1; end
if ~exist('if_cumulative') if_cumulative=1; end
umi_stats=cell(1,length(sessfiles));
%
for k=1:length(sessfiles)
    opts=psg_defopts();
    filename=cat(2,psg_data_dir,sessfiles{k});
    disp(' ');
    disp(sprintf('analyzing %s',filename));
    s=getfield(load(filename),'s');
    opts.cond_nstims=s.nstims;
    opts.cond_ncompares=size(s.sessions,2)-1;
    opts.cond_nsess=size(s.sessions,3);
    umi_stats{k}.filename=filename;
    umi_stats{k}.results=psg_tent_stats(s.sessions,setfields(opts,{'if_log','if_cumulative'},{if_log,if_cumulative}));    
end
    


