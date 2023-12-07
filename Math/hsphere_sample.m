function [pts,opts_used]=hsphere_sample(d,opts)
% [pts,opts_used]=hsphere_sample(d,opts) samples a hypersphere (surface)
% d: number of dimensions (3=ordinary sphere)
% opts: options
%  opts.method: 'random','axes','orthants','axes_and_orthants' 
%     defaults to random
%  opts.if_hemisphere: 0 (default) to sample whole sphere, 1 to sample only
%     a hemisphere (first nonzero coord is > 0)
%  opts.nsamps: number of points to sample (if random), otherwise set as follows:
%     for axes, 2*d/(1+if_hemisphere)
%     for orthants, 2^d/(1+if_hemisphere)
%     for axes_and_orthants, 2*d+2^2/(1+if_hemisphere)
%     for random, defaults to 2^d
%
% pts: array of size [nsamps d]
% opts_used: options used
%
%   See also:  FILLDEFAULT, INT2NARY.
%
if (nargin<2)
    opts=struct;
end
opts=filldefault(opts,'method','random');
opts=filldefault(opts,'nsamps',2^d);
opts=filldefault(opts,'if_hemisphere',0);
switch opts.method
    case 'axes'
        nsamps=(2*d)/(1+opts.if_hemisphere);
    case 'orthants'
        nsamps=(2^d)/(1+opts.if_hemisphere);
    case 'axes_and_orthants'
        nsamps=(2*d+2^d)/(1+opts.if_hemisphere);
    case 'random'
        nsamps=opts.nsamps;
    otherwise
        nsamps=0;
        warning(sprintf('method %s not recognized',opts.method));
        pts=[];
end
opts_used=opts;
if strfind(opts.method,'orthants')
    pts_orthants=fliplr([(1-2*int2nary([0:2^(d-opts.if_hemisphere)-1]',2,d))/sqrt(d)]);
end
if strfind(opts.method,'axes')
    pts_axes=eye(d);
    if opts.if_hemisphere==0
        pts_axes=[pts_axes;-eye(d)];
    end
end
switch opts.method
    case 'axes'
        pts=pts_axes;
    case 'orthants'
        pts=pts_orthants;
    case 'axes_and_orthants'
        if d>1
            pts=[pts_axes;pts_orthants];
        else
            pts=pts_axes;
            nsamps=2/(1+opts.if_hemisphere);
        end
    case 'random'
        pts=randn(nsamps,d);
        pts=pts./repmat(sqrt(sum(pts.^2,2)),[1 d]);
        if opts.if_hemisphere
            pts(:,1)=abs(pts(:,1));
        end
end
opts_used.nsamps=nsamps;
return
end