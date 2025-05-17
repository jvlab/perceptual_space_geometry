function [mapout,isgood] = btc_makeonemap(values,mapsize,pairnumber)

%%inputs are:correlation parameter values (1x2 pair), size(1x2 pair), 
%%pairnumber(1-45,scalar)

btc_surv_filename='btc_surv_13Aug09.mat';
s=load(btc_surv_filename);
pmethsum=getfield(s,'pmethsum');
dict=getfield(s,'dict');
clear s;

isgood=1;
ifshow=1;
verbose=0;
nmaps=1;
%
opts_augcoords=[];
%
opts_makemaps=[];
opts_makemaps.show=0;
opts_makemaps.area= mapsize;
opts_makemaps.nmaps=nmaps;
opts_makemaps.verbose=verbose;
%
npar=2;
coordnums_list=nchoosek([1:length(dict.codel)],npar);
allowed_pairs=[];
disp(' id expt   method(s).....');
for icomb=1:size(coordnums_list,1)
    exptname=btc_exptname(dict.codel(coordnums_list(icomb,:)),dict);
    if (isfield(pmethsum,exptname))
        allowed_pairs=[allowed_pairs,icomb];
    end
end
xy=values;
icomb=pairnumber;
    exptname=btc_exptname(dict.codel(coordnums_list(icomb,:)),dict);
    methods=getfield(pmethsum,exptname);
    for iv=1:npar
        icn(iv)=dict.codel(coordnums_list(icomb,iv));
    end
    ic=char(icn);
    maps=zeros(mapsize(1),mapsize(2));
    im=1; %use first method only
        tstring=sprintf('%2.0f->%s %s %s',icomb,exptname,methods{im}.name,methods{im}.variant_lab);
        disp(cat(2,'processing ',tstring));
            spec=[];
            for iv=1:npar
                spec=setfield(spec,ic(iv),xy(iv));
            end
            rv=btc_augcoords(spec,dict,opts_augcoords);
            if (rv.method{im}.ok_norm==1) & (rv.method{im}.ok_probs==1)
                [maps(:,:),ou,errs]=btc_makemaps(rv.method{im},opts_makemaps,dict);
                design=maps(:,:);
                if ~isempty(errs)
                    disp('cannot make a map with these values')
                    disp(errs);
                    disp(ou);
                end
            else
                disp('cannot find extreme param values')
                isgood=0;
            end
            %%%UNCOMMENT TO SHOW MAPS AS THEY'RE MADE
%         if (ifshow)
%             % now show the first map generated in each set
%             figure;
%             set(gcf,'Position',[200 200 800 800]);
%             set(gcf,'NumberTitle','off');
%             set(gcf,'Name',tstring);
%             imagesc(maps(:,:),[0 1]);
%             axis equal;axis tight;colormap('gray');
%             xlabel(sprintf('%2s=[%5.2f %5.2f]',ic,xy));
%             set(gca,'XTick',[]);
%             set(gca,'YTick',[]);
%         end %ifshow
        a=[];
        a.maps=maps;
        mapout=maps;
        %
    clear a
clear ic icn icomb iextreme im iv tstring xy
clear corrs counts counts_map counts_ubi
clear maps vecs vecs_emp entropy_mrfs entropy_mrfs_emp
clear imap indices methods ou outstring p2x2
clear rv spec which_pairs
