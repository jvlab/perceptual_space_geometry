%psg_rayang_summ: summarize ray angle statistics (see psg_visualize_demo for plotting)
%
% tabulates, for multiple datasets
%
%angles between positive and negative rays on each axis (using psg_rayfit with bid=1)
%angles between best-fitting line for each pair of axes (using psg_rayfit with bid=0)
%
%  See also: PSG_GET_COORDSETS, PSG_READ_COORDDATA, PSG_FINDRAYS, PSG_DEFOPTS, BTC_DEFINE,
%  PSG_RAYFIT, PSG_RAYANGLES, PSG_VISUALIZE_DEMO.
%
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_fit') opts_fit=struct(); end %for psg_rayfit
if ~exist('opts_ang') opts_ang=struct(); end %for psg_rayangles
if ~exist('opts_qpred') opts_qpred=struct(); end %for psg_qformpred
if ~exist('data_fullname') data_fullname=[]; end
if ~exist('setup_fullname') setup_fullname=[]; end
%
opts_read=filldefault(opts_read,'if_log',1);
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,opts_qpred);
nsets=length(ds);
d_rayfit=cell(nsets,1);
d_bidfit=cell(nsets,1);
ray_ends=cell(nsets,1);
bid_ends=cell(nsets,1);
angles_ray=cell(nsets,1);
angles_bid=cell(nsets,1);
opts_fit_used=cell(nsets,1);
opts_ang_used=cell(nsets,1);
ray_labels=cell(nsets,1);
ray_pair_labels=cell(nsets,1);
for iset=1:nsets
    nrays=rayss{iset}.nrays;
    disp(' ');
    disp(sprintf('set %1.0f: %s',iset,sets{iset}.label))
    disp(sprintf('stimulus coordinates group along %2.0f rays',nrays));
    %
    ray_counts=full(sparse(rayss{iset}.whichray(rayss{iset}.whichray>0),ones(sum(rayss{iset}.whichray>0),1),1,nrays,1));
    for iray=1:nrays
        disp(sprintf('ray %2.0f: %2.0f points; endpoint: %s',iray,ray_counts(iray),sprintf('%5.2f',rayss{iset}.endpt(iray,:))));
    end
    %find stimulus label at end of ray, in positive direction when posssible
    ray_labels{iset}=cell(1,nrays);
    for iray=1:nrays
        mults=rayss{iset}.mult(rayss{iset}.whichray==iray);
        maxend=intersect(find(abs(rayss{iset}.mult)==max(abs(mults))),find(rayss{iset}.whichray==iray));
        maxend=maxend(find(rayss{iset}.mult(maxend)==max(rayss{iset}.mult(maxend)))); %choose positive direction if possible
        ray_labels{iset}{iray}=strrep(strrep(sas{iset}.spec_labels{maxend},' ',''),'=',''); %strip = and space
        disp(sprintf('ray %2.0f label: %s',iray,ray_labels{iset}{iray}));       
    end
    ray_pair_labels{iset}=cell(1,nrays*(nrays-1)/2);
    ilab=0;
    for iray=1:nrays-1
        for jray=iray+1:nrays
            ilab=ilab+1;
            ray_pair_labels{iset}{ilab}=cat(2,ray_labels{iset}{iray},':',ray_labels{iset}{jray});
        end
    end
    disp(cat(2,' dim    ',sprintf('%20s ',ray_labels{iset}{:}),sprintf('%30s ',ray_pair_labels{iset}{:}))); %header
    %
    %for each dimension model, find best-fitting signed and unsigned rays, including the origin
    %
    dim_list=sets{iset}.dim_list;
    model_dim_max=max(dim_list);
    d_rayfit{iset}=cell(1,model_dim_max); %coordinate structure of best-fitting rays
    d_bidfit{iset}=cell(1,model_dim_max); %coordinate structure of best-fitting bidirectional rays
    ray_ends{iset}=cell(1,model_dim_max); %unidirectional ray directions at max
    bid_ends{iset}=cell(1,model_dim_max); %bidirectional  ray directions at max
    %
    angles_ray{iset}=cell(1,model_dim_max); %angle data, unidirectional
    angles_bid{iset}=cell(1,model_dim_max); %angle data, bidirectional
    %
    opts_fit_used{iset}=cell(model_dim_max,2); %d1: model dim, d2: 1+if_bid
    opts_ang_used{iset}=cell(model_dim_max,2);
    %
    for idimptr=1:length(dim_list) %compute ray fits and angles
        idim=dim_list(idimptr);
        [d_rayfit{iset}{idim},ray_ends{iset}{idim},opts_fit_used{iset}{idim,1}]=psg_rayfit(ds{iset}{idim},rayss{iset},filldefault(opts_fit,'if_bid',0));
        [d_bidfit{iset}{idim},bid_ends{iset}{idim},opts_fit_used{iset}{idim,2}]=psg_rayfit(ds{iset}{idim},rayss{iset},filldefault(opts_fit,'if_bid',1));
        [angles_ray{iset}{idim},opts_ang_used{iset}{idim,1}]=psg_rayangles(ray_ends{iset}{idim},sas{iset},rayss{iset},opts_ang);
        [angles_bid{iset}{idim},opts_ang_used{iset}{idim,2}]=psg_rayangles(bid_ends{iset}{idim},sas{iset},rayss{iset},opts_ang);
        %make a nice table
        t=sprintf('%3.0f   ',idim);
        for iray=1:nrays
            t=cat(2,t,sprintf(' %15.4f',angles_ray{iset}{idim}.cosangs(iray,iray,1,2)),'     '); %cosine of angle between pos and neg direction on each axis
        end
        for iray=1:nrays-1
            for jray=iray+1:nrays
                t=cat(2,t,sprintf(' %25.4f',angles_bid{iset}{idim}.cosangs(iray,jray)),'     '); %cosine of angle between biirectional fits of two axes
            end
        end
        disp(t);
    end
end
