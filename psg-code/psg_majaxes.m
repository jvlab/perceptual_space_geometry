function [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts)
% [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts) analyzes
%     an affine geometric model to determine major axes
%
%  d_ref, sa_ref: data and sa (labeling) structure for reference dataset, typically from psg_read_coorddata
%  d_adj, sa_adj: data and sa (labeling) structure for adjustable dataset, typically from psg_read_coorddata
%       sa_ref,sa_adj: only relevant field is typenames; typically should match but this is not checked
%  results_geo: results field from psg_geomodels_run
%  opts: options (can be empty)
%     opts.if_log: 1 to log
%     opts.model_class_list: list of model classes, defaults to {'affine','pwaffine'};
%        should also run for {'affine','procrustes'} but un-interesting
%        since eigenvalues should all be 1, and eigenvectors are degenerate
%     opts.tol_neg: tolerance for negative eigenvalues
%     opts.tol_match: tolerance for matching eigenvalues
%
%   consistency of results_geo (i.e., same models for each dimension) is not checked
%
%  results: analysis results
%  opts_used: options used
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE.
%
if nargin<=5 opts=struct; end
%
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'model_class_list',{'affine','pwaffine'});
opts=filldefault(opts,'tol_neg',10^-6);
opts=filldefault(opts,'tol_match',10^-4);
%
model_types_def=psg_geomodels_define();
opts.model_types_def=model_types_def;
opts.warnings=[];
results=cell(size(results_geo));
%
%scan the geometry results
nds_ref=size(results_geo,1);
nds_adj=size(results_geo,2);
d_ref_list=[];
d_adj_list=[];
for id_ref=1:nds_ref
    for id_adj=1:nds_adj
        res=results_geo{id_ref,id_adj};
        if ~isempty(res)
            d_ref_list=[d_ref_list,res.ref_dim];
            d_adj_list=[d_adj_list,res.adj_dim];
            model_types=res.model_types_def.model_types;
            ref_file=res.ref_file;
            adj_file=res.adj_file;
            results{id_ref,id_adj}.ref_dim=res.ref_dim;
            results{id_ref,id_adj}.adj_dim=res.adj_dim;
            results{id_ref,id_adj}.ref_file=ref_file;
            results{id_ref,id_adj}.adj_file=adj_file;
        end
    end %id_adj
end %id_ref
d_ref_list=unique(d_ref_list);
d_adj_list=unique(d_adj_list);
if opts.if_log
    disp(sprintf('reference dataset: %s',ref_file));
    disp('dimensions available:');
    disp(d_ref_list);
    disp('stimuli');
    disp(sa_ref.typenames);
    disp(sprintf('adjusted dataset: %s',adj_file));
    disp('dimensions available:');
    disp(d_adj_list);
    disp('stimuli');
    disp(sa_adj.typenames);
    disp('model types used for fitting geometric transformation from adjusted to reference');
    disp(model_types);
    disp('model classes to be anayzed');
    disp(opts.model_class_list);
end
%calculations
im_ptr=0;
for im=1:length(model_types)
    model_type=model_types{im};
    mclass=strmatch(model_types_def.(model_type).class,opts.model_class_list,'exact');
    if mclass>0 
        im_ptr=im_ptr+1;
        if opts.if_log
            disp(sprintf('analyzing model %25s (class: %s)',model_type,opts.model_class_list{mclass}))
        end   
        for id_ref=1:nds_ref
            for id_adj=1:nds_adj
                if ~isempty(res)
                    ref_dim=results{id_ref,id_adj}.ref_dim;
                    adj_dim=results{id_ref,id_adj}.adj_dim;
                    results{id_ref,id_adj}.model_types{im_ptr}=model_type;
                    transform_struct=results_geo{id_ref,id_adj}.transforms{im};
                    results{id_ref,id_adj}.ref.typenames=sa_ref.typenames;
                    results{id_ref,id_adj}.adj.typenames=sa_adj.typenames;
                    b=transform_struct.b;
                    Tlist=transform_struct.T;
                    npw=size(Tlist,3);
                    for ipw=1:npw
                        T=Tlist(:,:,ipw); %size of T is adj_dim x max(adj_dim,ref_dim),
                        %and if adj_dim> ref_dim, then columns ref_dim+1:end should be zero
                        num_eigs=adj_dim; % this is min(size(T)), number of eigenvalues computed
                        num_eigs_nz=min(ref_dim,adj_dim); %expected number of nonzero eigenvalues
                        for iar=1:2
                            switch iar
                                case 1
                                    A=T*transpose(T); %think of T as a stretching matrix M * a rotation matrix R
                                    lab='adj'; %compute eigenvals and eigenvecs of T*Ttranspose
                                    d_coords=d_adj{adj_dim}; %coordinates in dataset that is adjusted
                                    %A is square, size is adj_dim
                                case 2
                                    lab='ref'; %compute eigenvals and eigenvecs of Ttranspose*T
                                    A=transpose(T)*T; %think of T as a rotation matrix R * a stretching matrix M
                                    d_coords=d_ref{ref_dim}; %coordinates in reference dataset
                                    %A is square, size is max(adj_dim,ref_dim)
                            end
                            [eivecs,eivals,opts]=psg_majaxes_eigs(A,sprintf('%s [%s ipw %1.0f]',lab,model_type,ipw),ref_dim,adj_dim,opts);
                            results{id_ref,id_adj}.(lab).magnifs{im_ptr}(:,ipw)=b*sqrt(eivals); % since A=(MR)*transpose(MR)
                            results{id_ref,id_adj}.(lab).eivecs{im_ptr}(:,:,ipw)=eivecs;
                            results{id_ref,id_adj}.(lab).magnif_ratio{im_ptr}(:,ipw)=sqrt(eivals(1)/eivals(num_eigs_nz)); %ratio of highest to lowest magnification factor
                            prj_dim=size(d_coords,2); %number of dimensions to project onto, less than size(A) if ref_dim<adj_dim lab='ref'
                            %disp(sprintf('ref_dim %2.0f adj_dim %2.0f lab %s prj_dim %1.0f size(A) %2.0f %2.0f',ref_dim,adj_dim,lab,prj_dim,size(A)));
                            eivecs_project=eivecs(1:prj_dim,1:prj_dim); %if id_ref<id_adj, remaining eigenvecs are units
                            %rows of d_coords are the coordinates for each stimulus, either in adj set or ref set
                            results{id_ref,id_adj}.(lab).projections{im_ptr}(:,:,ipw)=d_coords*eivecs_project; %projections of original coordinates onto eigenvectors
                        end
                        %check that magnification factors agree (sqrt of eigenvals)
                        eig_diff=max(abs(results{id_ref,id_adj}.ref.magnifs{im_ptr}(1:num_eigs)-results{id_ref,id_adj}.adj.magnifs{im_ptr}(1:num_eigs)));
                        if eig_diff>opts.tol_match
                            warning_text=sprintf('magnifications disagree [%s ipw %1.0f] ref_dim %2.0f adj_dim %2.0f, disparity is %15.12f',model_type,ipw,ref_dim,adj_dim,eig_diff);
                            if opts.if_log
                                disp(warning_text);
                            end
                            opts.warnings=strvcat(opts.warnings,warning_text);
                        end
                    end %each transformationmatrix
                end %non-empty
            end %id_adj
        end %id_ref
    end %model type is affine
end
%
opts_used=opts;
return

function [eivecs,eivals,opts_new]=psg_majaxes_eigs(A,label,ref_dim,adj_dim,opts)
%compute, sort, and check eigenvalues and eigenvecs
[eivecs,eivals]=eig(A);
eivals=real(diag(eivals)); %A is self-adjoint
if any(eivals<-opts.tol_neg)
    warning_text=sprintf('negative eigenvalue in %s set to zero for ref_dim %2.0f adj_dim %2.0f: %15.12f',label,ref_dim,adj_dim,min(eivals));
    if opts.if_log
        disp(warning_text);
    end
    opts.warnings=strvcat(opts.warnings,warning_text);
end
opts_new=opts;
eivals=max(eivals,0);
[eivals,sort_inds]=sort(eivals,'descend'); %obtain eigenvalues in descending order
eivecs=real(eivecs(:,sort_inds));
return

