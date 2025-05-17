function [symclasses,if_warn,opts_used]=btc_symclasses(opts,btc_dict)
% [symclasses,if_warn,opts_used]=btc_symclasses(opts,btc_dict) sets up the symmetry classes of btc coordinates
%
% opts: structure defined by btc_soidfg_define, [] if not provided
%   uses the following fields:
%      coords:  the coordinate set to include in the model (e.g.,'gbcdetuvwa' or 'bcde')
%      sym_type: one of opts.sym_type_avail (see btc_soidfg_define)
%         Note that rot90 and full are equivalent here but not for btc_soidfg_model, or if there were more than 2 gray levels.
%  btc_dict: structure returned by btc_define; can be omitted
%   
% symclasses: a cell array of disjoint subsets of the elements of coords, along with the coordinates
%    they can map to by one of the specified symmetries.
%  if_warn: 1 if symmetry classes are not contained in coords
%  opts_used: opts, with defaults filled in by btc_soidfg_define
%
% 10Jun20:  Pass btc_dict to btc_hflip, btc_vflip, btc_rotcode, btc_hvi, btc_hvirev to speed up
%
%    See also:  BTC_SOIDFG_DEFINE, BTC_SOIDFG_MODEL, BTC_PAIRSNEEDED, BTC_HFLIP, BTC_VFLIP, BTC_HVI, BTC_HVIREV, BTC_ROTCODE,
%     BTC_SOIDFG_MODEL_TEST.
%
if (nargin==0)
    opts=[];
end
if (nargin<=1)
    btc_dict=btc_define([]);
end
opts=btc_soidfg_define(opts);
if_warn=0;
symclasses=cell(0);
%
if_used=zeros(1,length(opts.coords));
for ilet=1:length(opts.coords)
    if if_used(ilet)==0
        let=opts.coords(ilet);
        switch opts.sym_type
            case 'none'
                let_class=let;
            case 'hflip'
                let_class=unique([let,btc_hflip(let,btc_dict)]);
            case 'vflip'
                let_class=unique([let,btc_vflip(let,btc_dict)]);
            case 'hvflip'
                let_class=unique([let,btc_hflip(let,btc_dict),btc_vflip(let,btc_dict),btc_rotcode(let,2,btc_dict)]);
            case 'rot90'
                let_class=unique([let,btc_rotcode(let,1,btc_dict),btc_rotcode(let,2,btc_dict),btc_rotcode(let,3,btc_dict)]);
            case 'rot180'
                let_class=unique([let,btc_rotcode(let,2,btc_dict)]);
            case 'diag'
                let_dflip=unique([let,btc_hvi(let,btc_dict)]);
                let_class=unique([let_dflip,btc_hvirev(let_dflip,btc_dict)]);
            case 'full'
                let_rot=unique([let,btc_rotcode(let,1,btc_dict),btc_rotcode(let,2,btc_dict),btc_rotcode(let,3,btc_dict)]);
                let_class=unique([let_rot,btc_hflip(let_rot,btc_dict)]);
            otherwise
                error(sprintf('unrecognized symmetry type: %s',opts.sym_type));
        end
        for s=1:length(let_class)
            coord_ptr=find(opts.coords==let_class(s));
            if isempty(coord_ptr)
                warning(sprintf('btc %s has a symmetry class %s,  which has coordinates that are not members of the coordinate list %s',let,let_class,opts.coords));
                if_warn=1;
            end
            if_used(coord_ptr)=1;
        end
        symclasses{end+1}=let_class;
    end 
end %icoord
opts_used=opts;
return
