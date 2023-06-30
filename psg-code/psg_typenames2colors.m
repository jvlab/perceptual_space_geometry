function [rgb,symb,vecs,opts_used]=psg_typenames2colors(typenames,opts)
% [rgb,symb,vecs,opts_used]=psg_typenames2colors is a utility that assigns a plotting color and
% plotting symbol to a list of stimulus type names
%
% for btc: if multiple nonzero coord values are present, a sign is assigned only if they are all the same.
% if multiple btc letters are found, the color used is the mean.
%
% if the file in opts.faces_mpi_inventory_filename exists, it uses the
% variable in it, faces_mpi_attrib_info to determine whether typenames are faces, and if so, to parse them.
% if otps.faces_mpi_inventory_filename does not exist, typenames are assumed to be NOT mpi_faces
%
% typenames: cell array of type names, such as {'gp0133','gp0267','gp0400'}
% opts: options: can be omitted
%  opts.colors.[g,b,c,d,e,t,u,v,w,a]: colors to assign to each ray
%  opts.symbs.[z,m,p]: symbols to assign to zero, positive, negative
%
%  rgb: an rgb color triplet
%  symb: a plotting symbol
%  vecs: array of the coordinates found, NaN if unspecified
%    e.g., typenames= {'gp0300','bm0100ap0500','gm0400'} yields 
%   vecs=[...
%         0.30   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;...
%          NaN -0.10   NaN   NaN   NaN   NaN   NaN   NaN   NaN -0.50;...
%        -0.40   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN]
%
%  opts_used: options used
%
% 12Apr23: modify c color; fix bug in assigning symbols to bp0000cm0100 and similar; add symbols for pm and mp
% 29Jun23: begin capability for typenames from mpi faces. opts.type_class='btc' or 'mpi_faces', based on parsing
%
%    See also:  PSG_PLOTCOORDS, PSG_PLOTANGLES, PSG_READCOORD_DATA, PSG_FINDRAYS, BTC_DEFINE, PSG_TYPENAMES2COLORS_TEST,
%    PG_COLORS_LEGACY, FACES_MPI_INVENTORY.
%
if (nargin<2)
    opts=struct;
end
opts=psg_defopts(opts);
opts=filldefault(opts,'colors_nomatch',[0 0 0]);
opts=filldefault(opts,'symbs_nomatch','.');
%
%values if we find no matches
rgb=opts.colors_nomatch;
symb=opts.symbs_nomatch;
vecs=[];
%
ntn=length(typenames);
%
%determine whether typenames are mpi faces, i.e., of form {'132_y_f_n_a'}, and parse
%
if exist(opts.faces_mpi_inventory_filename,'file')
    underscore='_';
    if_mpi_faces=1;
    faces_mpi_attrib_info=getfield(load(opts.faces_mpi_inventory_filename),'faces_mpi_attrib_info');
    %faces_mpi_attrib_info
    attribs=faces_mpi_attrib_info.psg_order;
    nattribs=length(attribs);
    filename_order=cell(1,nattribs+1);
    for ia=1:nattribs
        filename_order{faces_mpi_attrib_info.(attribs{ia}).dim}=attribs{ia}; %order of attributes in file name;
    end
    filename_order{1}='indiv';
    att_table_num=NaN(ntn,nattribs+1);
    att_table_cell=cell(ntn,nattribs+1);
    att_table_avg=NaN(1,nattribs+1);
    for k=1:ntn
        unds=strfind(typenames{k},underscore);
        if length(unds)==nattribs
            %check that each attribute has a legitimate value and parse
            unds=[unds length(typenames{k})+1];
            indiv=str2num(typenames{k}(1:unds(1)-1));
            if ~isempty(indiv) %extract individual
                att_table_num(k,1)=indiv;
                att_table_cell{k,1}=indiv;
            end
            for ia=1:nattribs
                att_char=typenames{k}(unds(ia)+1:unds(ia+1)-1);
%                disp(sprintf('%2.0f %2.0f %s',k,ia,att_char));
%                disp(attribs{ia})
                if length(att_char)==1
                    idx=strmatch(att_char,faces_mpi_attrib_info.(filename_order{ia+1}).vals);
                    if length(idx)==1
                        att_table_num(k,ia+1)=idx;
                        att_table_cell{k,ia+1}=att_char;
                    end
                end
            end
        else
            if_mpi_faces=0;
        end
    end
    if any(isnan(att_table_num(:)))
        if_mpi_faces=0;
    else %now set up symbol and coords
        att_table_avg=mean(att_table_num,1);
    end
    opts.faces_mpi.att_table_avg=att_table_avg;
    opts.faces_mpi.att_table_num=att_table_num;
    opts.faces_mpi.att_table_cell=att_table_cell;
    opts.faces_mpi.att_table_order=filename_order;
    opts.faces_mpi.attrib_info=faces_mpi_attrib_info;
else
    if_mpi_faces=0;
    opts.faces_mpi=struct();
end
if (if_mpi_faces)
    opts.type_class='mpi_faces';
else
    opts.type_class='btc';
end
%determine symbols, colors, and vectors based on typename and type_class
switch opts.type_class
    case 'btc'
        dict=btc_define;
        codel=dict.codel;
        nbtc=length(codel);
        colors_def=struct;
        colors_def.g=[0.50 0.50 0.50];
        colors_def.b=[0.00 0.00 0.75];
        %colors_def.c=[0.25 0.25 1.00]; %modified 12Apr23
        colors_def.c=[0.00 0.80 0.80];
        colors_def.d=[0.00 0.75 0.00];
        colors_def.e=[0.00 1.00 0.10];
        %colors_def.t=[0.75 0.25 0.80];
        %colors_def.u=[0.50 0.25 0.80];
        %colors_def.v=[0.50 0.00 0.80];
        %colors_def.w=[0.75 0.00 0.80];
        colors_def.t=[0.85 0.60 0.30];
        colors_def.u=[0.75 0.65 0.20];
        colors_def.v=[1.00 0.90 0.20];
        colors_def.w=[0.85 0.75 0.20];
        colors_def.a=[1.00 0.00 0.00];
        %
        symbl='zmp';
        symbs_deg=struct;
        symbs_def.z='o';
        symbs_def.m='*';
        symbs_def.p='+';
        symbs_def.pm='v'; %downward triangle
        symbs_def.mp='^'; %upward triangle
        symbvals.z=0;
        symbvals.m=-1;
        symbvals.p=+1;
        opts=filldefault(opts,'colors',colors_def);
        opts=filldefault(opts,'symbs',symbs_def);
        %
        %replace any unspecified values
        for ibtc=1:nbtc
            if ~isfield(opts.colors,codel(ibtc))
                opts.colors.(codel(ibtc))=colors_def.(codel(ibtc));
            end
        end
        for isymb=1:length(symbl)
            if ~isfield(opts.symbs,symbl(isymb))
                opts.symbs.(symbl(isymb))=symbs_def.(symbl(isymb));
            end
        end
        if ~isfield(opts.symbs,'pm')
            opts.symbs.pm=symbs_def.pm;
        end
        if ~isfield(opts.symbs,'mp')
            opts.symbs.mp=symbs_def.mp;
        end
        %
        %values if we find single matches
        rgb=opts.colors_nomatch;
        symb=opts.symbs_nomatch;
        %
        nu=6;
        nc=2; %number of xchars before digits
        signs_found=[];
        lets_found=[]; %letters found for nonzero value
        vecs=[];
        for k=1:ntn %assume typename strings are in sets of nu(=6), like 'ap0300bm0100'
            tn=typenames{k};
            vec_new=NaN(1,nbtc);
            while length(tn)>=nu
                substr=tn(1:nu);
                if ismember(substr(1),codel) & ismember(substr(2),fieldnames(symbs_def))
                    val=symbvals.(substr(2))*str2num(substr(nc+1:nu))/(10^(nu-nc-1));
                    if val~=0
                        signs_found=[signs_found,substr(2)]; %moved after val~=0, 12Apr23
                        lets_found=[lets_found,substr(1)];
                        vec_new(find(dict.codel==substr(1)))=val;
                    end
                end
                tn=tn(nu+1:end);
            end
            if any(~isnan(vec_new))
                vecs=[vecs;vec_new];
            end
        %    vecs
        end
        lets_found=unique(lets_found);
        if length(lets_found)==0
            if ~isempty(signs_found) %an explicit zero was found
                symb=opts.symbs.z;
            end
        else
            signs_found_nz=setdiff(signs_found,'z');
            if length(signs_found_nz)==1
                symb=opts.symbs.(signs_found_nz);
            elseif strcmp(signs_found,'pm')
                symb=opts.symbs.pm;
            elseif strcmp(signs_found,'mp')
                symb=opts.symbs.mp;
            end
            if length(lets_found)==1
                rgb=opts.colors.(lets_found);
            else %average the colors
                rgbs=zeros(length(lets_found),3);
                for ifound=1:length(lets_found)
                    rgbs(ifound,:)=opts.colors.(lets_found(ifound));
                end
                rgb=mean(rgbs);
            end
        end
end %switch
opts_used=opts;
return
