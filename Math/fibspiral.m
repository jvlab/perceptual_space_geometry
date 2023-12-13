function [x,opts_used]=fibspiral(n,d,opts) 
% x=fibspiral(n,d) produces a sampling of points on the surface of a d-dimensional hypersphere
% n: number of points
% d: dimension.  Should work for any number of dimensions >=1
%    d=1: alternates between +1 and -1
%    d=2: equally spaced around a circle
%    d>=3: interesting
% opts: options
%   opts.nangs_min sets the minimal sampling to find the axial ("theta") angle
%   opts.nangs_mult sets the minimal factor of the number of points in
%       sampling the axial angle
%   opts.mults: a list of incommensurate numbers, at least of length d-1, otherwise it
%       will be filled in with the square roots of lowest primes
% x: array of size[n d] of the sampled points, magnitude 1
% opts_used: options used
%
% Method inspired by https://math.stackexchange.com/questions/3291489/can-the-fibonacci-lattice-be-extended-to-dimensions-higher-than-3
%
%  See also: FILLDEFAULT, HSPHERE_SAMPLE.
%
if (nargin<=2)
    opts=struct;
end
opts=filldefault(opts,'nangs_min',1024);
opts=filldefault(opts,'nangs_mult',16);
opts=filldefault(opts,'mults',sqrt([2 3]));
pmax=3;
while length(opts.mults)<(d-1)
    opts.mults=sqrt(primes(pmax));
    pmax=2*pmax;
end
opts.mults=opts.mults([1:d-1]);
opts_used=opts;
%    
x=zeros(n,d);
%
switch d
    case 1
      x=(-1).^[0:n-1]';
    case 2
        thetas=2*pi*(1/2+[0:n-1])/n;
        x=[cos(thetas(:)) sin(thetas(:))];
    otherwise
        x=ones(n,d); %we will cumulatively multiply factors for the opening-up angles
        mult_list=opts.mults;
        nangs=max(opts.nangs_min,n*max(opts.nangs_mult,4*d));
        theta_list=pi*(1/2+[0:nangs-1])/nangs;
        for a=1:d-2 %determine the axial angles
            %choose angles equally weighted according to sin^(d-1-a)
            wts=cumsum(sin(theta_list).^(d-1-a));
            wts=[0 wts(1:end-1)/wts(end)]; %wts ranges from 0 to 1
            %find critical values for weights
            wt_crits=mod(mult_list(a)*([0:n-1]-0.5),1);
            thetas=zeros(1,n); %points to the angles wereh wts first exceeds wt_crits
            for ipt=1:n
                thetas(ipt)=theta_list(sum(double(wts<wt_crits(ipt))));
            end
            x(:,a)=x(:,a).*cos(thetas(:));
            x(:,a+1:end)=x(:,a+1:end).*repmat(sin(thetas(:)),1,d-a);
        end
        %final coordinates
        phis=mult_list(d-1)*2*pi*(1/2+[0:n-1])/n;
        x(:,d-1)=x(:,d-1).*cos(phis(:));
        x(:,d)=x(:,d).*sin(phis(:));
end
return
