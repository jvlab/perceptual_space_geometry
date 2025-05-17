%btc_soidfg_model_test: test btc_soidfg_model, btc_symclasses.
%
%    See also:  BTC_SOIDFG_MODEL, BTC_SOIDFG_DEFINE, BTC_SYMCLASSES, BTC_DEFINE.
%
if ~exist('opts_def')
    opts_def=[];
end
dict=btc_define;
opts_def=btc_soidfg_define(opts_def);
if ~exist('coords_lists') coords_lists={dict.codel,'gbcde','bcde','tuvw'}; end
model_types=opts_def.model_type_avail;
model=cell(length(opts_def.sym_type_avail),length(coords_lists),length(model_types));
symclasses=cell(length(opts_def.sym_type_avail),length(coords_lists),length(model_types));
for isym_type=1:length(opts_def.sym_type_avail)
    sym_type=opts_def.sym_type_avail{isym_type};
    disp(' ')
    disp(sprintf('symmetry type %s',sym_type))
    for icoord=1:length(coords_lists)
        coords=coords_lists{icoord};
        disp(' ')
        for imodel_type=1:length(model_types)
            model_type=model_types{imodel_type};
            opts=setfields(opts_def,{'coords','sym_type','model_type'},{coords,sym_type,model_type});
            [symclasses{isym_type,icoord,imodel_type},if_warn]=btc_symclasses(opts);
            model{isym_type,icoord,imodel_type}=btc_soidfg_model(opts);
            disp(sprintf(' coord list %10s, model type %15s: number of classes: %3.0f, number of params: %3.0f',...
                coords,model_type,length(symclasses{isym_type,icoord,imodel_type}),model{isym_type,icoord,imodel_type}.nparams));
        end %imodel_type
    end %icoord
end %isym_type
