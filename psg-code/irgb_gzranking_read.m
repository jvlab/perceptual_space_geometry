function [s_parsed,s_raw,fn_used]=irgb_gzranking_read(filename)
%[s_parsed,s_raw]=irgb_gzranking_read(filename) reads and parses the Giesel-Zaidi ranking structure
%
% filename: full path and file name, requested if empty
%
% s_parsed: parsed data
% s_raw: raw data
% fn_used: filename used
%
%   See also:  IRGB_GZRANKING_ANALYZE.
%
filename_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat';
n_imgs=261; %total number of images (not all ranked)
if nargin<1
    filename=[];
end
if isempty(filename)
    fn_used=getinp('file name','s',[],filename_def);
else
    fn_used=filename;
end
s_raw=load(fn_used);
%
props=fieldnames(rmfield(s_raw.Ranking,'README'))';
n_props=length(props);
subjs=fieldnames(s_raw.Ranking.(props{1}))';
n_subjs=length(subjs);
%
s_parsed=struct;
s_parsed.README=s_raw.Ranking.README;
s_parsed.props=props;
s_parsed.subjs=subjs;
s_parsed.dims={'imageno','subjs','props'};
s_parsed.rankings_all=zeros(n_imgs,n_subjs,n_props);
ranked_list=[]; %list of image numbers actually ranked
for iprop=1:n_props
    for isubj=1:n_subjs
        ranked_list_check=[];
        lists=s_raw.Ranking.(props{iprop}).(subjs{isubj});
        n_ranks=length(lists);
        total_ranked=0;
        for irank=1:n_ranks
            vals=cell2mat(lists{irank});
            ranked_list_check=unique([ranked_list_check,vals]);
            total_ranked=total_ranked+length(vals);
            if any(s_parsed.rankings_all(vals,isubj,iprop)>0)
                disp(sprintf('property %12s, subject %3s: rank %3.0f duplicates some previous ranks',props{iprop},subjs{isubj},irank));
            else
                s_parsed.rankings_all(vals,isubj,iprop)=irank;
            end
        end %irank
        %check that ranked list is consistent across subjects
        if isempty(ranked_list)
            ranked_list=ranked_list_check;
        else
            if length(ranked_list)~=length(ranked_list_check)
                disp('mismatch of ranked list lengths');
            elseif any(ranked_list~=ranked_list_check)
                disp('mismatch of ranked list');
            end
        end
        disp(sprintf('property %12s, subject %3s: %3.0f images ranked in %3.0f ranks',props{iprop},subjs{isubj},total_ranked,n_ranks));
    end %isubj
end %iprop
s_parsed.ranked_list=ranked_list;
s_parsed.rankings=s_parsed.rankings_all(ranked_list,:,:);
return
