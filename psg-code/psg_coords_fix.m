%psg_coords_fix: utility to fix a file of coords that has an extra entry in stim_labels
%mat-file (e.g., ' bc6pt_coords_MC_sess01_10.mat') needs to be loaded.
disp(' fixing coordinates as in bc6pt_coords_MC_sess01_10_fix');
whos
stim_labels_orig=stim_labels
stim_labels=[];
for k=1:25;stim_labels(k,:)=char(stim_labels_orig(k,:));end
stim_labels=char(stim_labels)
clear k
clear stim_labels_orig
whos
disp('now save the file.');

