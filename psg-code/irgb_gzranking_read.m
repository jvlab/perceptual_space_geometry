function [s_parsed,s_raw,fn_rank,fn_img]=irgb_gzranking_read(filename_rank,filename_img)
%[s_parsed,s_raw,fn_rank,fn_img]=irgb_gzranking_read(filename_rank,filename_img) reads and parses the Giesel-Zaidi ranking structure
%
% filename_rank: full path and file name for ranking data file, requested if empty
% filename_img: full path and file name for image name file, requested if empty
%
% s_parsed: parsed data
% s_raw: raw data
% fn_rank: filename used for ranking data
% fn_img: filename used for image list
%
%   See also:  IRGB_GZRANKING_ANALYZE.
%
n_imgs=261; %total number of images (not all ranked)
filename_rank_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/RankingData.mat';
filename_img_def='C:/Users/jdvicto/Documents/jv/ENCL/Zaidi/ImageNames.mat';
%
if nargin<1
    filename_rank=[];
end
if isempty(filename_rank)
    fn_rank=getinp('file name for ranking data','s',[],filename_rank_def);
else
    fn_rank=filename_rank;
end
%
if nargin<2
    filename_img=[];
end
if isempty(filename_img)
    fn_img=getinp('file name for image list','s',[],filename_img_def);
else
    fn_img=filename_img;
end
s_raw=load(fn_rank);
s_raw.R=getfield(load(fn_img),'R');
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
s_parsed.image_names_all_orig=s_raw.R.imfiles;
s_parsed.image_names_all_short=strrep(strrep(s_raw.R.imfiles,'Images/',''),'/cut/','/');
s_parsed.image_names_ranked=cell(length(ranked_list),1);
for imgno=1:length(ranked_list) %names of the images actually ranked
    s_parsed.image_names_ranked{imgno}=s_parsed.image_names_all_short{ranked_list(imgno)};
end
return
