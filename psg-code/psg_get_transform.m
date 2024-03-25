function [transform,opts_used]=psg_get_transform(dim_max,opts)
% [transform,opts_used]=psg_get_transform(dim_max,opts) gets a specified linear transformation
%
% dim_max: maximum dimension of transform
%
% opts.if_prompt: 1 (default) for a prompt
% opts.type: 'simple' (default) scale, permute, and offset, only option at present
% opts.if_permute: 1 to allow for a permutation
% opts.if_offset: 1 (default) to allow for scaling
% (scaling and permute
%
% transform: structure to specify the transformation, as in psg_geo_procrustes
%   transform.b: scalar, will be 1 if opts.if_scale=0
%   transform.T: matrix of size [dim_max,dim_max]
%   transform.c: offset, row of size [dim_max 1], will be zeros if ref and adj are centered
%   new=transform.b*old*transform.T+repmat(transform.c,npts,1)
%
%   See also: PSG_COORD_PIPE_PROC, PSG_GEO_PROCRUSTES, FILLDEFAULT.
if (nargin<=1)
    opts=struct;
end
opts=filldefault(opts,'type','simple');
opts=filldefault(opts,'if_prompt',1);
opts=filldefault(opts,'if_permute',1);
opts=filldefault(opts,'if_offset',1);
opts_used=opts;
%
transform=struct;
transform.T=eye(dim_max);
transform.b=1;
transform.c=zeros(1,dim_max);
if_ok=0;
while (if_ok==0)
    switch opts.type
        case 'simple'
        if opts.if_prompt
            disp('You are providing a simple linear transformation. Each coordinate is rescaled');
            if opts.if_permute
                disp('then optionally permuted')
            end
            if opts.if_offset
                disp('then optionally offset')
            end
        end
        %
        scales=getinp(sprintf('up to %2.0f scale factors (will be cycled if necessary)',dim_max),'f',[-Inf Inf],1);
        scales=scales(1+mod([1:dim_max]-1,length(scales))); %cyclically reuse
        transform.T=diag(scales);
        %
        if opts.if_permute
            if_okperm=0;
            while (if_okperm==0)
                disp('In the permutation, the kth element indicates the new dimension for old dimension k');
                p=getinp(sprintf('a permutation of [1 ... %1.0f]',dim_max),'d',[1 dim_max],[1:dim_max]);
                if length(p)==dim_max
                    if all(sort(p(:))==[1:dim_max]')
                        if_okperm=1;
                    end
                end
            end
            pmatrix=zeros(dim_max);
            for id=1:dim_max
                pmatrix(id,p(id))=1;
            end
            transform.T=transform.T*pmatrix;
        end
        %
        if opts.if_offset
            offsets=getinp(sprintf('up to %2.0f offset (will be cycled if necessary)',dim_max),'f',[-Inf Inf],0);
            offsets=offsets(1+mod([1:dim_max]-1,length(offsets))); %cyclically reuse
            transform.c=offsets;
        end
    otherwise
        disp('unknown transformation type, transformation set to identity');
    end
    %
    disp('T (post-multiplying matrix');
    disp(transform.T);
    disp('b (scale factor');
    disp(transform.b);
    disp('offset');
    disp(transform.c);
    %
    if_ok=getinp('1 if ok','d',[0 1])
end
return
