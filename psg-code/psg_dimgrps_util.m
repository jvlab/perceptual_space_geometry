function [ndgs,dim_groups,dim_sources,tstring_dimspecs]=psg_dimgrps_util(dim_max)
% [ndgs,dim_groups,dim_sources,tstring_dimspecs]=psg_dimgrps_util(dim_max) is a utility
% to subdivde a total number of dimensions into a subset of increasing dimensions, used for nesting
%
% ndgs: number of dimension groups
% dim_groups: specifier: 1-d array of length<=dim_max, summing to dim_max
%  dim_groups=[dim_max] (the default) nests all dimensions together
%  dim_groups=ones(1,dim_max) makes a model of dimension k to be the first k dimensions of higher models
%  dim_groups=[2 3] makes models 1-2 the first 2 dimensions of the 2-d model, models 3-5 the first dimensions of the 5-d model
% dim_sources: 1-d array of length dim_max, dim_sources(id) is the source dimension to use for the model of dimension id
% tstring_dimspecs: descriptor string
%
% See also: PSG_COORD_PIPE_PROC.
%
ifok=0;
while (ifok==0)
    dim_groups=[];
    tstring_dimspecs='';
    ndgs=0;
    tstring_dimspec=cell(0);
    ifok=1;
    while sum(dim_groups)<dim_max
        ndgs=ndgs+1;
        dim_groups(ndgs)=getinp(sprintf('size of dimension group %1.0f, which starts at dimension %1.0f',ndgs,1+sum(dim_groups)),...
            'd',[1 dim_max-sum(dim_groups)],dim_max-sum(dim_groups));
        tstring_dimspec{ndgs}=cat(2,'d[',sprintf('%1.0f ',[(1+sum(dim_groups(1:end-1))):sum(dim_groups)]),']');
        tstring_dimspecs=cat(2,tstring_dimspecs,tstring_dimspec{ndgs},' ');
    end
    tstring_dimspecs=deblank(tstring_dimspecs);
    disp(sprintf('dimension groups: %s',tstring_dimspecs))
    ifok=getinp('1 if ok','d',[0 1],ifok);
end
dim_sources=[zeros(1,dim_max)];
for id=1:dim_max
    dgp=min(find(id<=cumsum(dim_groups)));
    dim_sources(id)=sum(dim_groups(1:dgp));
end
return
