% psg_typenames2colors_test is a test for psg_typenames2colors
%
%    See also:  BTC_DEFINE.
%
dict=btc_define;
codel=dict.codel;
nbtc=length(codel);
mplet={'m','p'};
itn=0;
figure;
hl=cell(0);
ht=[];
if ~exist('opts') opts=struct;end
tn_all=cell(1,length(codel)*2);
for icodel=1:nbtc
    for imp=1:2
        itn=itn+1;
        tn_base=cat(2,codel(icodel),mplet{imp});
        tn=cat(2,tn_base,'0500');
        tn_all{itn}=tn;
        disp(sprintf('typename %2.0f: %s',itn,tn));
        [rgb,symb,vecs,opts_used]=psg_typenames2colors({tn},opts);
        hp=plot(2*imp-3+[-0.5 0 0.5],repmat(icodel,1,3),'LineWidth',1);
        hl=[hl;hp];
        ht=strvcat(ht,tn_base);
        hold on;
        set(hp,'color',rgb);
        set(hp,'Marker',symb);
        set(gca,'XLim',[-2 4]);
        set(gca,'YLim',[-0.5 0.5]+[1 nbtc]);
    end
    set(gca,'YDir','reverse');
end
legend(hl,ht);
[rgb_all,symb_all,vecs_all,opts_used_all]=psg_typenames2colors(tn_all,opts);
vecs_all


