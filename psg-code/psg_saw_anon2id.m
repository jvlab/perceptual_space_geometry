%psg_saw_anon2id.m: one-time routine to rename files from SAW study (JN 2024)
%to confirm with subject ID convention
%
%sets creation date to date of execution
%
if ~exist('dir_anon') dir_anon='C:\Users\jdvicto\Downloads\WaraichCoordsAnonymized';end
if ~exist('dir_id') dir_id='C:\Users\jdvicto\Documents\jv\EY7977\psg\psg_data';end
%
%from SAW's email
% MC S1
% BL S2
% EFV  S3
% SJ S4
% SAW  S5
% NK S6
% YCL  S7
% SA S8
% JF S9
% AJ S10
% SN S11
% ZK S12
% CME  S13
disp(sprintf('files from %s will be renamed and copied to %s',dir_anon,dir_id))
if ~exist('X')
    X={'MC','BL','EFV','SJ','SAW','NK','YCL','SA','JF','AJ','SN','ZK','CME'};
end
for ix=1:length(X)
    subj=X{ix};
    disp(' ');
    disp(sprintf(' subject %2.0f: %s',ix,subj));
    dir_anon_cmd=cat(2,'dir ',dir_anon,'\*coords_',sprintf('S%1.0f',ix),'.mat  /w/b');
    [errs,files]=dos(dir_anon_cmd);
    while ~isempty(files)
        [fn_anon,files]=strtok(files);
        if ~isempty(fn_anon)
            fn_id=strrep(fn_anon,sprintf('S%1.0f',ix),subj);
            disp(sprintf('%30s -> %30s',fn_anon,fn_id));
            copy_cmd=cat(2,'copy ',dir_anon,filesep,fn_anon,' ',dir_id,filesep,fn_id);
            errs_copy=dos(copy_cmd);
            if errs_copy>0
                disp(errs_copy);
            end
            setdate_cmd=cat(2,'copy /b ',dir_id,filesep,fn_id,' +,,'); %sets the file date but not the modified date
            errs_setdate=dos(setdate_cmd);
            if errs_setdate>0
                disp(errs_setdate);
            end
        end
    end
end

