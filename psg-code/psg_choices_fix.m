%psg_choices_fix: utility to fix a file of coords that has an extra entry in stim_list
%mat-file (e.g., ' tvpm3pt_choices_MC_sess01_10.mat') needs to be loaded.
disp(' fixing stim_list as in tvpm3pt_choices_MC_sess01_10');
whos
stim_list_orig=stim_list
stim_list=[];
for k=1:25;stim_list(k,:)=char(stim_list_orig(k,:));end
stim_list=char(stim_list)
clear k
clear stim_list_orig
whos
disp('now save the file.');

