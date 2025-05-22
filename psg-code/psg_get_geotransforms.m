function [transforms,dims_avail,desc,opts_used]=psg_get_geotransforms(opts)
%[transforms,dims_avail,desc,opts_used]=psg_get_geotransforms(opts) gets a set of geometric transforms, typically from a file
%
% opts: options for file name, path name, model type, etc.
%   opts.if_uigetfile: 1 (default) to get file name via graphical interface
%     fields that are missing or NaN are prompted
%
% transforms{id}
%   structure to specify the transformation from a model of dimension id, as in psg_geomodels_apply, psg_geo_procrustes and psg_pwaffine_apply.
%   For linear transformations (additional params for projective, piecewise affine)
%    transforms{id}.b: scalar -- for compatibility with procrustes, but will always be 1 (scale absorbed into T)
%    transforms{id}.T: matrix of size [dim_max,dim_max]
%    transforms{id}.c: offset, row of size [dim_max 1], will be zeros if ref and adj are centered
%    new=transforms{id}.b*old*transforms{id}.T+repmat(transforms{id}.c,npts,1)
%  Note this is non-empty for id in dims_avail
% desc: text descriptor
% dims_avail: dimensions for which a transformation is available
% opts_used: options used
%   contains all the fields of opts, with interactive entries, and also
%   opts_used.fullname: full file name
%
%   See also: PSG_COORD_PIPE_PROC, PSG_GEOMODELS_APPLY, FILLDEFAULT, PSG_GET_TRANSFORM, HLID_GEOM_TRANSFORM_STATS, PSG_COORD_PIPE_PROC,
%     PSG_GEOMODELS_DEFINE.
%
if (nargin<1)
    opts=struct;
end
opts=filldefault(opts,'if_uigetfile',1);
opts=filldefault(opts,'pathname',NaN);
opts=filldefault(opts,'pathname_default',cat(2,'.',filesep));
opts=filldefault(opts,'filename',NaN);
opts=filldefault(opts,'fieldname',NaN);
opts=filldefault(opts,'fieldname_default','results');
opts=filldefault(opts,'model_type',NaN);
opts=filldefault(opts,'geo_dims',{'sub','preproc','embed'}); %labels for dimensions of results.geo from hlid_geom_transform_stats
ngd=length(opts.geo_dims);
labels=cell(1,ngd);
for id=1:ngd
    opts=filldefault(opts,cat(2,'geo_choice_',opts.geo_dims{id}),NaN);
end
model_types_def=psg_geomodels_define();
%
transforms=cell(0);
dims_avail=[];
desc=[];
%
ifok=0;
while (ifok==0)
    if isnan(opts.filename)
        if opts.if_uigetfile
            [filename,pathname]=uigetfile('*.mat','file containing transformations');
            opts.pathname=pathname;
            opts.filename=filename;
        else
            opts.pathname=getinp('path to file containing transformations','s',[],strrep(opts.pathname_default,filesep,'/'));
            opts.filename=getinp('file containing transformations','s',[]);
        end
        opts.filename=cat(2,opts.filename,'.mat');
        opts.filename=strrep(opts.filename,'.mat.mat','.mat');
        opts.fullname=cat(2,opts.pathname,filesep,opts.filename);
    end
    if isnan(opts.fieldname)
        opts.fieldname=getinp('field name containing transformation','s',[],opts.fieldname_default);
    end   
    if exist(opts.fullname,'file')
        f=load(opts.fullname);
        if isfield(f,opts.fieldname)
            r=f.(opts.fieldname);
            ifok=1;
        else
            disp(sprintf('field %s does not exist',opts.fieldname));
            opts.fieldname_default=opts.fieldname;
            opts.fieldname=NaN;
            opts.filename=NaN;
        end
    else
        disp(sprintf('file %s not found',opts.filename));
        opts.filename=NaN;
    end
    if ifok==0 & opts.if_uigetfile
        getinp('anything to proceed','d',[0 1]);
    end
end
%find the structure or substructure that contains transformatons at {ref_dim,adj_dim}
rt=[];
desc=[];
%
if iscell(r) %typical output of psg_geomodel_run
    rt=r;
elseif isfield(r,'geo') %typical ouptut of hlid_geom_transform_stats
    gv_choice=zeros(1,ngd);
    for id=1:ngd
        geo_dim=opts.geo_dims{id};
        gv_max(id)=r.(cat(2,'n',geo_dim,'s'));
    end
    ifok=0;
    while (ifok==0)
        if_ask=0;
        gv_index=0;
        labels=cell(1,ngd);
        for id=1:ngd
            geo_dim=opts.geo_dims{id};
            gv=cat(2,'geo_choice_',geo_dim);
            labels{id}=r.(cat(2,geo_dim,'_labels'));
            for k=1:length(labels{id})
                if isempty(labels{id}{k})
                    labels{id}{k}='[empty]';
                end
            end
            if isnan(opts.(gv))
                if_ask=1;
                for k=1:gv_max(id)
                    disp(sprintf(' choice %1.0f for %s: %s',k,geo_dim,labels{id}{k}));
                end
                gv_choice(id)=getinp('choice','d',[1 gv_max(id)]);
             else
                gv_choice(id)=opts.(gv);
            end
            %construct a multi-index
            if (id==1)
                gv_index=gv_choice(id);
            else
                gv_index=gv_index+prod(gv_max(1:id-1))*(gv_choice(id)-1);
            end
        end
        if (if_ask)
            disp(' selection:')
            for id=1:ngd
                geo_dim=opts.geo_dims{id};
                disp(sprintf(' choice for %s is %1.0f->%s',geo_dim,gv_choice(id),labels{id}{gv_choice(id)}));
            end
%            disp(sprintf('multi-index: %3.0f',gv_index));
            ifok=getinp('1 if ok','d',[0 1]);
         else
            ifok=1;
        end
    end
    rt=r.geo{gv_index};
    %make a descriptor
    for id=1:ngd
        geo_dim=opts.geo_dims{id};
        gv=cat(2,'geo_choice_',geo_dim);
        opts.(gv)=gv_choice(id);
        desc=cat(2,desc,'geo_dim=',labels{id}{gv_choice(id)},' ');
    end
    desc=cat(2,' (',deblank(desc),')');
end
%
if ~iscell(rt)
    warning('no transforms found');
    transforms=[];
    opts_used=opts;
    return
else
    %scan results for diagonal (ref_dim=adj_dim) terms, show ref file and adj file, and model types
    rt_diag=cell(1);
    dims_avail=[];
    model_list=cell(0);
    if_mtypes_ok=1;
    for ir=1:size(rt,1)
        for ia=1:size(rt,2)
            t=rt{ir,ia};
            if isfield(t,'ref_dim') & isfield(t,'adj_dim')
                if t.ref_dim==t.adj_dim
                    dims_avail=[dims_avail,t.ref_dim];
                    mtypes=t.model_types_def.model_types;
                    rt_diag{t.ref_dim}=t.transforms;
                    if isempty(model_list)
                        model_list=mtypes;
                    end
                    for k=1:length(mtypes)
                        if ~strcmp(mtypes{k},model_list{k})
                            if_mtypes_ok=0;
                        end
                    end
                end %on diag
            end %isfield
        end %ia
    end %ir
    if ~if_mtypes_ok
        warning('model type mismatch');
        disp(mtypes);
        disp(model_list);
        transforms=[];
        opts_used=opts;
    end
    model_type=opts.model_type;
    model_type_ptr=strmatch(model_type,model_list,'exact');
    if length(model_type_ptr)~=1
       model_type=NaN;
    end
    if isnan(model_type)
        for k=1:length(model_list)
            disp(sprintf('%2.0f-> model type %s',k,model_list{k}));
        end
        model_type_ptr=getinp('choice','d',[1 length(model_list)]);
    end
    opts.model_type=model_list{model_type_ptr};
    opts.model_class=model_types_def.(opts.model_type).class;
    transforms=cell(1,max(max(dims_avail),max(size(rt))));
    for k=1:length(dims_avail)
        transforms{dims_avail(k)}=rt_diag{dims_avail(k)}{model_type_ptr};
    end
    desc=cat(2,model_list{model_type_ptr},desc);
    opts_used=opts;
end
return
