%btc_atg_alt_demo: look at alternative constructions for alpha, theta, gamma textures
%
%b=c, and they determine d,u,w
%g,a,e,t,v,a specified to be Pickard
%
%   See also:  BTC_DEFINE, BTC_AUGCOORDS.
%
dict=btc_define;
nsamps=getinp('number of samples in [0 1]','d',[10 10^6],10^3);
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
for ib=1:n
    b=bvals(ib);
    d_plus_u=(g+b).^2/(1+g);
    d_minus_u=(g-b).^2/(1-g);
    d=(d_plus_u+d_minus_u)/2;
    u=(d_plus_u-d_minus_u)/2;
    %
    spec=struct;
    %here, symmetrize: c=b, e=d, t=u=v=w
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
    if all(p2x2(:)>=0)
        b_valid(ib)=1;
        corrs=getcorrs_p2x2(p2x2);
        entropy(ib)=corrs.entropy;
        pickard_NWSE(ib)=max(max(abs(corrs.cig_conds(:,[1 4]))));
        pickard_NESW(ib)=max(max(abs(corrs.cig_conds(:,[2 3]))));
        dvals(ib)=d;
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
plot(bvals,[b_valid;pmins;entropy;pickard_NWSE;pickard_NESW;dvals]'),
legend({'b_v_a_l_i_d','p_m_i_n_s','entropy','Pickard NWSE','Pickard NESW','b_d_i_a_g'},'Location','Best');
xlabel('b');
set(gca,'XLim',[-1 1]);
set(gca,'YLim',[-1 1]);
if if_specify_t
    title(sprintf(' a=%7.5f g=%7.5f t=%7.5f v=%7.5f',a,g,t,v))
else
    title(sprintf(' a=%7.5f g=%7.5f, t determined by symmetry',a,g))

end
