function [ifok,matches]=psg_geomodels_apply_util(y,f1,f2,samediff,tol)
% [ifok,matches]=psg_geomodels_apply_util(y,f1,f2,samediff,tol) is a utility for
% psg_geomodels_apply_test
%
%   See also: PSG_GEOMODELS_APPLY_TEST, PSG_PWAFFINE_APPLY, PERSP_APPLY.
%
%
string_ans={'OK','NOT OK'};
diffs=abs(y.(f1)(:)-y.(f2)(:));
switch samediff
    case 'same'
        ans_ind=1+double(max(diffs)>tol);
        ifok=double(ans_ind==1);
    case 'diff'
        ans_ind=1+double(max(diffs)<=tol);
        ifok=double(ans_ind==2);
end
matches=reshape(double(diffs<=tol),size(y.(f1),1),size(y.(f1),2));
disp(sprintf('comparing %30s with %30s, s.b. %6s: %8s maxdiff=%15.8f, ndiff=%6.0f nsame=%6.0f',f1,f2,samediff,string_ans{ans_ind},max(diffs),sum(diffs>tol),sum(diffs<=tol)));
return
