function [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts)
% [results,opts_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts) analyzes
%     an affine geometric model to determine major axes
%
%  d_ref, sa_ref: data and sa (labeling) structure for reference dataset, typically from psg_read_coorddata
%  d_adj, sa_adj: data and sa (labeling) structure for adjustable dataset, typically from psg_read_coorddata
%  results_geo: results field from psg_geomodels_run
%  opts: options (can be empty)
%     opts.if_log: 1 to log
%     opts.model_class_list: list of model classes, defaults to {'affine','pwaffine'};
%        should also run for {'affine','procrustes'} but un-interesting
%        since eigenvalues should all be 1, and eigenvectors are degenerate
%     opts.tol: tolerance for negative eigenvalues
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
opts=filldefault(opts,'tol',10^-6);
%
model_types_def=psg_geomodels_define();
opts.model_types_def=model_types_def;
opts.warnings=[];
results=cell(size(results_geo));
%
%scan the results
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
    disp(sprintf('adjusted dataaset: %s',adj_file));
    disp('dimensions available:');
    disp(d_adj_list);
    disp('model types')
    disp(model_types);
end
%calculations
im_ptr=0;
for im=1:length(model_types)
    model_type=model_types{im};
    if strmatch(model_types_def.(model_type).class,opts.model_class_list,'exact')>0
        im_ptr=im_ptr+1;
        if opts.if_log
            disp(sprintf('analyzing %s',model_type))
        end   
        for id_ref=1:nds_ref
            for id_adj=1:nds_adj
                if ~isempty(res)
                    results{id_ref,id_adj}.model_types{im_ptr}=model_type;
                    transform_struct=results_geo{id_ref,id_adj}.transforms{im};
                    b=transform_struct.b;
                    Tlist=transform_struct.T;
                    npw=size(Tlist,3);
                    for ipw=1:npw
                        T=Tlist(:,:,ipw); %T is square,number of rows = adj dim, number of nonzero cols=ref dim
                        %
                        A=T*transpose(T); %think of T as a stretching matrix M * a rotation matrix R
                        %compute eigenvals and eigenvecs of T*Ttranspose,
                        [eivecs,eivals,opts]=psg_majaxes_eigs(A,'adj',id_ref,id_adj,opts);
                        results{id_ref,id_adj}.magnifs_adj{im_ptr}(:,ipw)=b*sqrt(eivals); % since A=(MR)*transpose(MR)
                        results{id_ref,id_adj}.eivecs_adj{im_ptr}(:,:,ipw)=eivecs;
                        %
                        A=transpose(T)*T; %think of T as a rotation matrix R * a stretching matrix M
                        %compute eigenvals and eigenvecs of Ttranspose*T
                        [eivecs,eivals,opts]=psg_majaxes_eigs(A,'ref',id_ref,id_adj,opts);
                        results{id_ref,id_adj}.magnifs_ref{im_ptr}(:,ipw)=b*sqrt(eivals); % since A=transpose(MR)*MR
                        results{id_ref,id_adj}.eivecs_ref{im_ptr}(:,:,ipw)=eivecs;
                        %
                        %magnif_ref and magnifs_adj should be identical other than possibly trailing zeros

                    end
                end
            end %id_adj
        end %id_ref
    end %model type is affine
end
%
opts_used=opts;
return

function [eivecs,eivals,opts_new]=psg_majaxes_eigs(A,label,id_ref,id_adj,opts)
%compute, sort, and check eigenvalues and eigenvecs
[eivecs,eivals]=eig(A);
eivals=real(diag(eivals)); %A is self-adjoint
opts_new=opts;
if any(eivals<-opts.tol)
    warning_text=sprintf('negative eigenvalue in %s calc set to zero for id_ref %2.0f id_adj %2.0f: %15.12f',label,id_ref,id_adj,min(eivals));
    if opts.if_log
        disp(warning_text);
    end
    opts_new.warnings=strvcat(opts_new.warnings,warning_text);
end
eivals=max(eivals,0);
[eivals,sort_inds]=sort(eivals,'descend'); %obtain eigenvalues in descending order
eivecs=real(eivecs(:,sort_inds));
return

