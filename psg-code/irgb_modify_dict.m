function manip_dict=irgb_modify_dict()
% manip_dict=irgb_modify_dict() creates a structure that gives the
% correspondence between image manipulations and the images returned by irgb_modify
%
% manip_dict: a structure whose fields describe characteristics of the image manipulations
%
%   See also: IRGB_MODIFY, IRGB_PSG_IMGS_SETUP, IRGB_PSG_SESS_SETUP.
%
manip_dict=struct;
manip_dict.orig_bw.modify_name='gray';
manip_dict.bw_whiten.modify_name='whitened';
manip_dict.bw_randph.modify_name='randphase';
manip_dict.bw_randph.multiples=1;
manip_dict.filt_bw.modify_name='filt';
manip_dict.filt_bw_whiten.modify_name='whitened_filt';
manip_dict.filt_bw_randph.modify_name='randphase_filt';
manip_dict.filt_bw_randph.multiples=1;
%
fns=fieldnames(manip_dict);
for ifn=1:length(fns)
    fn=fns{ifn};
    if ~isfield(manip_dict.(fn),'multiples')
        manip_dict.(fn).multiples=0;
    end
end
return
