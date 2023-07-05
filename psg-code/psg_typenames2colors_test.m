% psg_typenames2colors_test is a test for psg_typenames2colors
%
% 04Jul23: add a test for faces_mpi
%
%    See also:  BTC_DEFINE, PSG_TYPENAMES2COLORS, FACES_MPI_INVENTORY.
%
%
if ~exist('opts')
    opts=struct;
end
opts=psg_defopts(opts);
%
%btc test
disp('testing renderings for btc');
dict=btc_define;
codel=dict.codel;
nbtc=length(codel);
mplet={'m','p'};
itn=0;
figure;
set(gcf,'Name','btc');
set(gcf,'NumberTitle','off');
hl_btc=cell(0);
ht_btc=[];
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
        hl_btc=[hl_btc;hp];
        ht_btc=strvcat(ht_btc,tn_base);
        hold on;
        set(hp,'color',rgb);
        set(hp,'Marker',symb);
        set(gca,'XLim',[-2 4]);
        set(gca,'YLim',[-0.5 0.5]+[1 nbtc]);
    end
    set(gca,'YDir','reverse');
end
legend(hl_btc,ht_btc);
[rgb_all,symb_all,vecs_all,opts_used_all]=psg_typenames2colors(tn_all,opts);
vecs_all
%
disp('testing renderings for faces_mpi');
if ~exist(opts.faces_mpi_inventory_filename,'file')
    disp(sprintf('cannot find faces inventory in %s',opts.faces_mpi_inventory_filename));
else
    figure;
    set(gcf,'Name','faces_mpi');
    set(gcf,'NumberTitle','off');
    set(gcf,'Position',[150 50 1000 950]);
    hl_faces_mpi=cell(0);
    ht_faces_mpi=[];
    faces_mpi_attrib_info=getfield(load(opts.faces_mpi_inventory_filename),'faces_mpi_attrib_info');
    %generate a list of type names
    iface=0;
    typenames_faces_mpi=cell(0);
    xpos=0;
    for ig=1:faces_mpi_attrib_info.gender.nlevels
        gender_char=faces_mpi_attrib_info.gender.vals(ig);
        for ia=1:faces_mpi_attrib_info.age.nlevels
            age_char=faces_mpi_attrib_info.age.vals(ia);
            ypos=0;
            for is=1:faces_mpi_attrib_info.set.nlevels
                set_char=faces_mpi_attrib_info.set.vals(is);
                    for ie=1:faces_mpi_attrib_info.emo.nlevels
                    emo_char=faces_mpi_attrib_info.emo.vals(ie);
                    ypos=ypos+1;
                    iface=iface+1;
                    typenames_faces_mpi{iface}=cat(2,'000_',age_char,'_',gender_char,'_',emo_char,'_',set_char);
                    [rgb_faces_mpi,symb_faces_mpi]=psg_typenames2colors(typenames_faces_mpi(iface),opts);
                    hp=plot(xpos+[-0.2 0.2],repmat(ypos,1,2),'LineWidth',1);
                    hl_faces_mpi=[hl_faces_mpi;hp];
                    ht_faces_mpi=strvcat(ht_faces_mpi,cat(2,age_char,gender_char,emo_char,set_char));
                    hold on;
                    set(hp,'color',rgb_faces_mpi);
                    set(hp,'Marker',symb_faces_mpi);
                    set(gca,'XLim',[-1 14]);
                    set(gca,'YLim',[0.5 12.5]);
                end %ie
            end %is
            xpos=xpos+1;
        end %ia
    end %ig
    xlabel('gender, age');
    ylabel('set, emo');
    set(gca,'YDir','reverse');
    legend(hl_faces_mpi,ht_faces_mpi,'FontSize',7,'Location','Best');
    %disp(typenames_faces_mpi)
end



