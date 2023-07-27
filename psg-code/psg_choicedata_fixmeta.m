%psg_choicedata_fixmeta: script to fix the metadata in choice data files
%
% see psg_read_choicedata_notes.docx
%
fn=getinp('file name to fix (no .mat)','s',[]);
load(fn)
disp(sprintf('responses is %4.0f x %2.0f',size(responses)));
for k=1:size(responses_colnames,1)
    responses_colnames(k,:)=strrep(responses_colnames(k,:),'>','<');
end
fn_new=getinp('file name to write','s',[],cat(2,fn,'_fixmeta'));%
save(fn_new,'desc','responses','responses_colnames','stim_list');
disp(sprintf('%s written',fn_new));
