function [planes_needed,complete]=btc_plan(planes_have,if_show)
% [planes_needed,complete]=btc_plan(planes_have,if_show) dtermines which planes are needed for
% fitting quadratic forms
%
% planes_have: a cell array of two-character designations (btc_exptname convention) of planes
%   that one already has data in, defaults to blank
% if_show: 1 to make a nice display on the console, defaults to 1
%
% planes_needed: planes_needed.letcodes is a cell array of the planes still needed to fit
%   a quadratic form based on letcode
% complete: each field of complete is a letcode of a quadratic form for which all planes
%   needed are in planes_have.  complete.letcodes is a cell array of the required planes.
%  
%   See also:  BTC_PAIRSNEEDED, BTC_DEFINE.
%
if (nargin==0) planes_have=[]; end
if (nargin<=1) if_show=1; end
planes_needed=[];
complete=[];
%
dict=btc_define;
val_a='a ';
val_g='g ';
list_bhv={'   ','b  ','bc '};
list_bdi={'   ','d  ','de '};
list_theta={'     ','t    ','tu   ','tv   ','tuvw '};
allplanes=sortrows(btc_pairsneeded(dict.codel,dict));
%
%remove planes_have from allplanes
%
imiss=ones(size(allplanes,1),1);
for iall=1:size(allplanes,1)
    if strmatch(allplanes(iall,:),char(planes_have))>0
        imiss(iall)=0;
    end
end
allplanes=allplanes(find(imiss==1),:);
%
headers=([allplanes,repmat(' ',size(allplanes,1),1)])';
headers=headers(:)';
if (if_show)
    disp('Table of planes needed for various quadratic forms');
end
for i_theta=1:length(list_theta)
    if (if_show)
        if length(deblank(list_theta{i_theta}))>0
            disp(sprintf(' thetas used: %s',list_theta{i_theta}));
        else
            disp(' thetas used: [none]');
        end
        disp(sprintf('                 %s',headers));
    end
    for i_bdi=1:length(list_bdi)
        for i_bhv=1:length(list_bhv)
            s=[];
            uletcodes=cat(2,val_a,list_theta{i_theta},list_bhv{i_bhv},list_bdi{i_bdi},val_g);
            letcodes=strrep(uletcodes,' ','');
            exptname=btc_exptname(letcodes,dict);
            needed=btc_pairsneeded(exptname,dict);
            matches=zeros(1,size(allplanes,1));
            for ineeded=1:size(needed,1)
                pneeded=strmatch(needed(ineeded,:),allplanes);
                if (pneeded>0)
                    matches(pneeded)=matches(pneeded)+1;
                    s=strvcat(s,needed(ineeded,:));
                end
            end
            if isempty(s)
                planes_needed.(exptname)=cell(0);
            else
                planes_needed.(exptname)=cellstr(s)';
            end
            mstring=repmat(' ',1,3*size(allplanes,1));
            mstring(3*[1:size(allplanes,1)])='0';
            mstring(3*find(matches>0))='1';
            if (if_show)
                disp(sprintf('%15s %s  total planes needed: %2.0f',uletcodes,mstring,sum(matches)));
            end
            if sum(matches)==0
                complete.(exptname)=cellstr(needed)';
            end
        end %bhv
    end %bdi
end %theta
return
