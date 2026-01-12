%btc_atg_alt_demo: look at alternative constructions for alpha, theta, gamma textures
%
%b=c, and they determine d,u,w
%g,a,e,t,v,a specified to be Pickard
%
% first does an exhaustive analysis of a particular (a,g) pair, optionally specifying one theta
%
% then does a survey of the (a,g) range, with:
%   theta and one beta-diag chosen to be Pickard, other beta-diag forced to match ; all thetas equal, beta_h=beta_v chosen by maxent
%   theta and one beta-diag chosen to be Pickard, other beta-diag forced to match ; two thetas=0,  beta_h=beta_v chosen by maxent
%   theta and one beta-diag chosen to be Pickard, other beta-diag forced to match ; two thetas are neg,  beta_h=beta_v chosen by maxent
%   thetas set to gamma^3, betas set to gamma^2 (standard btc_augcoords)
%
%  Note that the maxent solution when two thetas are zero could allow for the beta-diags to be unequal, but, if all four thetas are the same, or
%   if the 90-deg rotated thetas are their negatives, then the beta-diags for maxent must be equal.
%
% Also see symmetric_pickard_notes.docx
%
%   See also:  BTC_DEFINE, BTC_AUGCOORDS, BTC_PAIRDOMAIN_SYMPICK.
%
dict=btc_define;
nsamps=getinp('number of samples in [0 1] for search for b','d',[10 10^6],10^3);
nsamps_surv=getinp('number of samples in [0 1] for a,g survey','d',[2 1000],10);
if ~exist('mapsize') mapsize=64;end
%
a=getinp('alpha','f',[-1 1],0);
g=getinp('gamma','f',[-1 1],0);
if_specify_t=getinp('1 to specify t, 0 to determine by symmetry','d',[0 1]);
if if_specify_t
    t=getinp('theta','f',[-1 1],g^3);
    v=getinp('theta-opposite','f',[-t t],0);
    disp('d_opposite will be equal to d');
end
n=2*nsamps+1;
bvals=[-nsamps:nsamps]/nsamps;
b_valid=zeros(1,n);
pmins=zeros(1,n);
entropy=zeros(1,n);
pickard_NWSE=zeros(1,n);
pickard_NESW=zeros(1,n);
dvals=NaN(1,n);
uvals=NaN(1,n);
for ib=1:n
    b=bvals(ib);
    d_plus_u=(g+b).^2/(1+g);
    d_minus_u=(g-b).^2/(1-g);
    d=(d_plus_u+d_minus_u)/2;
    u=(d_plus_u-d_minus_u)/2;
    %
    spec=struct;
    %symmetrize: c=b, e=d, t=u=v=w
    spec.g=g;
    spec.b=b;
    spec.c=b;
    spec.d=d;
    spec.e=d;
    if if_specify_t       
        spec.t=t;
        spec.v=v;
    else
        spec.t=u;
        spec.v=u;
    end
    spec.u=u;
    spec.w=u;
    spec.a=a;
    %
    vec=btc_letcode2vec(spec,dict);
    corrs=btc_vec2corrs(vec,dict);
    p2x2=getp2x2_corrs(corrs);
    pmins(ib)=min(p2x2(:));
    dvals(ib)=d;
    uvals(ib)=u;
    if all(p2x2(:)>=0)
        b_valid(ib)=1;
        corrs=getcorrs_p2x2(p2x2);
        entropy(ib)=corrs.entropy;
        pickard_NWSE(ib)=max(max(abs(corrs.cig_conds(:,[1 4]))));
        pickard_NESW(ib)=max(max(abs(corrs.cig_conds(:,[2 3]))));
    end
end
if any(b_valid>0)
    bmin_ptr=min(find(b_valid==1));
    bmax_ptr=max(find(b_valid==1));
    disp(sprintf('range of validity for b: %7.4f to %7.4f',bvals(1,[bmin_ptr bmax_ptr])));
    if any(b_valid(bmin_ptr:bmax_ptr)==0)
        disp('discontinuous range');
    end
    bmaxent_ptr=find(entropy==max(entropy));
    disp(sprintf('max entropy =%7.4f at b=%7.4f',entropy(bmaxent_ptr),bvals(bmaxent_ptr)));
else
    disp('no valid values for b')
end
figure;
set(gcf,'Position',[100 100 1200 800]);
plot(bvals,[b_valid;pmins;entropy;pickard_NWSE;pickard_NESW;dvals;uvals]'),
legend({'If valid','min prob','entropy','Pickard NWSE','Pickard NESW','b_d_i_a_g','u and w'},'Location','Best');
xlabel('b');
set(gca,'XLim',[-1 1]);
set(gca,'YLim',[-1 1]);
if if_specify_t
    title(sprintf(' a=%7.5f g=%7.5f t=%7.5f v=%7.5f',a,g,t,v))
else
    title(sprintf(' a=%7.5f g=%7.5f, t determined by symmetry',a,g))

end
%
%now do a survey, with t always determined by symmetry
%
nts=4;
ns=2*nsamps_surv+1;
surv_vals=[-nsamps_surv:nsamps_surv]/nsamps_surv;
blo=NaN(ns,ns,nts);
bhi=NaN(ns,ns,nts);
maxent=NaN(ns,ns,nts);
%also will want to compute brange, and entropy of the canonical (g,a)
for it=1:nts
    switch it
        case 1
            t_type='sym';
            tmult=1;
        case 2
            t_type='zero';
            tmult=0;
        case 3
            t_type='antisym';
            tmult=-1;
        case 4
            t_type='gamma^3';
    end
    disp(sprintf('doing survey for t type: %s',t_type));
    for ia=1:ns
        a=surv_vals(ia);
        for ig=1:ns
            g=surv_vals(ig);
            b_valid=zeros(1,n);
            entropy=zeros(1,n);
            if it~=4
                for ib=1:n
                    b=bvals(ib);
                    d_plus_u=(g+b).^2/(1+g);
                    d_minus_u=(g-b).^2/(1-g);
                    d=(d_plus_u+d_minus_u)/2;
                    u=(d_plus_u-d_minus_u)/2;
                    %
                    spec=struct;
                    %always symmetrize b and c, 
                    spec.g=g;
                    spec.b=b;
                    spec.c=b;
                    spec.d=d;
                    spec.e=d;
                    spec.t=u*tmult;
                    spec.v=u*tmult;
                    spec.u=u;
                    spec.w=u;
                    spec.a=a;
                    %
                    vec=btc_letcode2vec(spec,dict);
                    corrs=btc_vec2corrs(vec,dict);
                    p2x2=getp2x2_corrs(corrs);
                    pmins(ib)=min(p2x2(:));
                    dvals(ib)=d;
                    uvals(ib)=u;
                    if all(p2x2(:)>=0)
                        b_valid(ib)=1;
                        corrs=getcorrs_p2x2(p2x2);
                        entropy(ib)=corrs.entropy;
                    end               
                    if any(b_valid>0)
                        bmin_ptr=min(find(b_valid==1));
                        bmax_ptr=max(find(b_valid==1));
                        bmaxent_ptr=min(find(entropy==max(entropy)));
                        blo(ia,ig,it)=bvals(bmin_ptr);
                        bhi(ia,ig,it)=bvals(bmax_ptr);
                        maxent(ia,ig,it)=entropy(bmaxent_ptr);
                    end
                end
            else % t=gamma^3
                spec.g=g;
                spec.b=g^2;
                spec.c=g^2;
                spec.d=g^2;
                spec.e=g^2;
                spec.t=g^3;
                spec.v=g^3;
                spec.u=g^3;
                spec.w=g^3;
                spec.a=a;
                %
                vec=btc_letcode2vec(spec,dict);
                corrs=btc_vec2corrs(vec,dict);
                p2x2=getp2x2_corrs(corrs);
                if all(p2x2>=0)
                    blo(ia,ig,it)=spec.b;
                    bhi(ia,ig,it)=spec.b;
                    corrs=getcorrs_p2x2(p2x2);
                    maxent(ia,ig,it)=corrs.entropy;
                end
            end
        end %ig
    end %ia
    figure;
    set(gcf,'Position',[100 100 1000 800]);
    for isub=1:4
        subplot(2,2,isub)
        switch isub
            case 1
                vplot=blo;
                tstring='b_m_i_n';
                yr=[-1 1];
            case 2
                vplot=bhi;
                tstring='b_m_a_x';
                yr=[-1 1];
            case 3
                vplot=bhi-blo;
                tstring='b_m_a_x - b_m_i_n';
                yr=[0 2];
            case 4
                vplot=maxent;
                tstring='entropy';
                yr=[0 1];
        end
        imagesc(surv_vals,surv_vals,vplot(:,:,it),yr);
        set(gca,'YDir','normal');
        xlabel('g');
        ylabel('a');
        axis tight;
        axis equal;
        colorbar;
        title(cat(2,tstring,' (t: ', t_type,')'));
    end
end %it

