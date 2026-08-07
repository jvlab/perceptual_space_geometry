%rigidity_complete_demo: examine characteristics of random rigidity matrices from
%complete graphs
%
% See also:  INCIDENCE2RIGIDITY
%
if ~exist('if_frozen') if_frozen=1; end %1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed
if ~exist('dim_list') dim_list=[2:10]; end
if ~exist('nstims_list') nstims_list=[2:10]; end
if ~exist('ndraws') ndraws=100; end
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
recip_conds=zeros(length(dim_list),length(nstims_list),ndraws); %reciprocal condition numbers, highest eig/lowest
recip_part_ratios=zeros(length(dim_list),length(nstims_list),ndraws); %sum(abs(eigs))^2/(sum of eigs^2)
ranks=zeros(length(dim_list),length(nstims_list),ndraws);
%
for istim_ptr=1:length(nstims_list)
    nstims=nstims_list(istim_ptr);
    combs=nchoosek([1:nstims],2);
    imtx=zeros(nstims,size(combs,1)); %generate a complete graph
    for ie=1:size(combs,1)
        imtx(combs(ie,:),ie)=1;
    end
    for idim_ptr=1:length(dim_list)
        d=dim_list(idim_ptr);
        for idraw=1:ndraws
            rmtx=incidence2rigidity(imtx,d);
            [v,eivals]=eig(rmtx*rmtx');
            eivals=diag(eivals);
            recip_conds(idim_ptr,istim_ptr,idraw)=min(eivals)/max(eivals);
            recip_part_ratios(idim_ptr,istim_ptr,idraw)=sum(eivals).^2/sum(eivals.^2);
            ranks(idim_ptr,istim_ptr,idraw)=rank(rmtx);
        end
        disp(sprintf('did analysis with %2.0f stimuli and dimension %2.0f, %6.0f draws',nstims,d,ndraws));
    end %idim_ptr
end %istim_ptr

