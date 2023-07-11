function params_new=irgb_spec_modify(spec_params)
% params_new=irgb_spec_modify(spec_params) modifies a stimulus specification structure for independently distributed rgb stimuli (irgb)
%
% See also:  IRGB_PSG_SETUP, IRGB_SPEC_MAKE, IRGB_STIM_MAKE, IRGB_SPEC_DEFAULTS.
%
nrgb=3;
%
if (nargin<=0)
    spec_params=struct;
end
% fill defaults
spec_params=irgb_spec_defaults(spec_params);
spec_params_orig=spec_params;
ifok=0;
while (ifok==0)
    spec_params.paradigm_name=getinp('paradigm name','s',[],spec_params.paradigm_name);
    spec_params=irgb_spec_util(spec_params,'paradigm_type');
    switch spec_params.paradigm_type
        case 'spokes'
            spec_params.mean_mults=getinp('one or more mean multipliers (number of points on a ray)','f',[-1 1],spec_params.mean_mults);
            nmean_dirs_old=size(spec_params.mean_dirs,1);
            nmean_dirs=getinp('number of mean directions (rays)','d',[1 36],nmean_dirs_old);
            for k=1:nmean_dirs
                vname=cat(2,'mean_dir_',sprintf('%1.0f',k));
                if k<=nmean_dirs_old
                    p.(vname)=spec_params.mean_dirs(k,:);
                else
                    p.(vname)=zeros(1,nrgb);
                end
                spec_params.mean_dirs(k,:)=getfield(irgb_spec_getrgb(p,vname),vname);
            end
            spec_params.mean_include_zero=getinp('1 to include a mean of zero','d',[0 1],spec_params.mean_include_zero);
            spec_params=irgb_spec_getrgb(spec_params,'mean_offset');
            spec_params=irgb_spec_util(spec_params,'cov_mode');
            spec_params.cov_mults=getinp('one or more covariance multipliers','f',[-1 1],spec_params.cov_mults);
        case 'distributions'
            %allow for changing distribution_weights, distribution_count
        otherwise
            warning(sprintf('unknown paradigm type (%s)',spec_params.paradigm_type));
    end
    ifok=getinp('1 if ok, -1 to revert (discard changes)','d',[-1 1]);
    if (ifok==-1)
        spec_params=spec_params_orig;
    end
end
params_new=spec_params;
return

function pnew=irgb_spec_util(p,vname)
vlist=cat(2,vname,'_list');
for k=1:length(p.(vlist));
    disp(sprintf('%s %1.0f->%s',strrep(vname,'_',' '),k,p.(vlist){k}));
end
vptr=strmatch(p.(vname),p.(vlist),'exact');
vptr=getinp('choice','d',[1 length(p.(vlist))],vptr);
pnew=p;
pnew.(vname)=p.(vlist){vptr};
return

function pnew=irgb_spec_getrgb(p,vname)
nrgb=3;
ifok=0;
pnew=p;
while ifok==0
    pnew.(vname)=getinp(sprintf('%s, %1.0f vals',strrep(vname,'_',' '),nrgb),'f',[-1 1],p.(vname));
    if any(size(pnew.(vname))~=[1 nrgb])
        disp(sprintf(' must have %1.0f entries in [-1 1]',nrgb));
    else
        ifok=1;
    end
end
return
