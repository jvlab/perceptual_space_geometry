%psg_cond_sess_split: script to split existing session configuration (csv) files into sub-sessions
%
% See also:  IRGB_PSG_SETUP, IRGB_PSG_SESS_SETUP, PSG_SPOKES_SETUP, PSG_COND_WRITE.
%
if ~exist('s') 
   setup_file=getinp('setup file (and path), e.g., *.mat','s',[],[]);
   s=getfield(load(setup_file),'s');
end
disp(s);
opts_psg=psg_defopts(s.opts_psg); %take from opts_psg as used by irgb_psg_sess_setup
nexamps=1+max(s.examps_used(:));
nstims=size(s.typenames,1);
nsess=length(s.session_cells);
ntrials=size(s.sessions,1);
%
disp(sprintf('number of  stimuli: %5.0f',nstims));
disp(sprintf('number of sessions: %5.0f',nsess));
disp(sprintf('number of trials per session: %5.0f',ntrials));
%
split_list=[];
while sum(split_list)~=ntrials
    if isempty(split_list)
        split_list=[floor(ntrials/2) ntrials-floor(ntrials/2)];
    end
    split_list=getinp(sprintf('split session sizes, must sum to %1.0f',ntrials),'d',[1 ntrials],split_list);
end
if isfield(s,'paradigm_name')
    filename_base=getinp('file name base (and path), _sess[#][a-z].csv will be appended for cond files','s',[],s.paradigm_name);
else
    filename_base=getinp('file name base (and path), _sess[#][a=z].csv will be appended for cond files','s',[],'split');
end
%
for isess=1:opts_psg.cond_nsess
    for isplit=1:length(split_list)
        split_range=[1+sum(split_list(1:isplit-1)) sum(split_list(1:isplit))];
        append=char(double('a')+isplit-1);
        disp(sprintf(' session %2.0f: split %2.0f (%s): trials %4.0f to %4.0f',isess,isplit,append,split_range));
        filename=cat(2,filename_base,'_sess',zpad(isess,opts_psg.sess_zpad),append);
        psg_cond_write(filename,s.session_cells{isess}(split_range(1):split_range(2),:),setfield(opts_psg,'if_log',1));
    end
end
