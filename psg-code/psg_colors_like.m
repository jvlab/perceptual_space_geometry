function c=psg_colors_like
% c=psg_colors_like sets up default colors and symbols for likelihood plots
%
% c: structure of default colors for paradigms, and symbols for subjects
%
% 12Aug23: added colors for btc selected paradigms, from psg_typenames2colors
% 23Aug23: added colors for bright paradigms
%
% See also: PSG_LIKE_ANALTABLE, PSG_TYPENAMES2COLORS.
%
paradigm_colors=struct;
paradigm_colors.texture=     [1.0 0.0 0.0];
paradigm_colors.texture_like=[0.5 0.4 0.0];
paradigm_colors.intermediate_texture=paradigm_colors.texture_like;
paradigm_colors.image_like=  [0.0 0.5 0.0];
paradigm_colors.intermediate_object=paradigm_colors.image_like;
paradigm_colors.image=       [0.0 0.0 0.7];
paradigm_colors.word=        [0.4 0.0 0.4];
%
paradigm_colors.bgca3pt=     [1.0 0.0 0.0];
paradigm_colors.bc6pt=       [0.7 0.0 0.7];
paradigm_colors.bcpm3pt=     [0.0 0.0 0.7];
paradigm_colors.bdce3pt=     [0.0 1.0 0.0];
paradigm_colors.tvpm3pt=     [0.6 0.6 0.0];
%
%btc selected paradigms
paradigm_colors.bc6pt_b_rand= psg_typenames2colors({'bp1000'});
paradigm_colors.bc6pt_bm_rand=psg_typenames2colors({'bp1000'});
paradigm_colors.bc6pt_bp_rand=psg_typenames2colors({'bp1000'});
%
paradigm_colors.bc6pt_c_rand= psg_typenames2colors({'cp1000'});
paradigm_colors.bc6pt_cm_rand=psg_typenames2colors({'cp1000'});
paradigm_colors.bc6pt_cp_rand=psg_typenames2colors({'cp1000'});
%
paradigm_colors.bgca3pt_a_rand= psg_typenames2colors({'ap1000'});
paradigm_colors.bgca3pt_am_rand=psg_typenames2colors({'ap1000'});
paradigm_colors.bgca3pt_ap_rand=psg_typenames2colors({'ap1000'});
%
paradigm_colors.bgca3pt_b_rand= psg_typenames2colors({'bp1000'});
paradigm_colors.bgca3pt_bm_rand=psg_typenames2colors({'bp1000'});
paradigm_colors.bgca3pt_bp_rand=psg_typenames2colors({'bp1000'});
%
paradigm_colors.bgca3pt_c_rand= psg_typenames2colors({'cp1000'});
paradigm_colors.bgca3pt_cm_rand=psg_typenames2colors({'cp1000'});
paradigm_colors.bgca3pt_cp_rand=psg_typenames2colors({'cp1000'});
%
paradigm_colors.bgca3pt_g_rand= psg_typenames2colors({'gp1000'});
paradigm_colors.bgca3pt_gm_rand=psg_typenames2colors({'gp1000'});
paradigm_colors.bgca3pt_gp_rand=psg_typenames2colors({'gp1000'});
%
%bright paradigm
paradigm_colors.c02s08_c01=[0.0 0.0 1.0];
paradigm_colors.c02s08_c02=[0.0 0.5 0.8];
paradigm_colors.c02s08_cXs=[0.0 0.0 0.0];
paradigm_colors.c02s08_s01Xs03Xs06Xs08=[1.0 0.0 0.0];
%
%reserved subject symbols, anomalous subjects in the word domain are unfilled: BL (intelligence), EFV (phonetics), and SA (pet)
subj_symbs_res.MC='s';
subj_symbs_res.SAW='d';
subj_symbs_res.ZK='^';
subj_symbs_res.BL='^';
subj_symbs_res.EFV='<';
subj_symbs_res.SA='>';
subj_symbs_res.AJ='o';
subj_symbs_res.SJ='v';
subj_symbs_res.SN='p';
subj_symbs_res.YCL='d';
subj_symbs_res.CME='x';
subj_symbs_res.JF='+';
subj_symbs_res.NK='.';
%
subj_fills_res.MC=1;
subj_fills_res.SAW=1;
subj_fills_res.ZK=1;
subj_fills_res.BL=0;
subj_fills_res.EFV=0;
subj_fills_res.SA=0;
subj_fills_res.AJ=1;
subj_fills_res.SJ=1;
subj_fills_res.SN=1;
subj_fills_res.YCL=1;
subj_fills_res.CME=0;
subj_fills_res.JF=0;
subj_fills_res.NK=0;
%
subj_symbs_unres='ovsdph<>h'; %other available symbols
subj_fills_unres=[zeros(1,6) ones(1,3)]; %fill options
apriori_symb='*';
%
c.paradigm_colors=paradigm_colors;
c.subj_symbs_res=subj_symbs_res;
c.subj_fills_res=subj_fills_res;
c.subj_symbs_unres=subj_symbs_unres;
c.subj_fills_unres=subj_fills_unres;
c.apriori_symb=apriori_symb;
end

