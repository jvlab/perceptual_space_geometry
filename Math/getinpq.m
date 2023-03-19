function [ansnew,newval]=getinpq(answers,varnam,varprompt,vartype,varlim,vardef)
% [ansnew,newval]=getinpq(answers,varnam,varprompt,vartype,varlim,vardef] is a
% utility routine to add a field "varnam" to the structure ansnew.
%    if the field occurs in answers with a non-blank value, it is used.
%    if the field does not occur in answers, the default value vardef is used.
%    if the field occurs with a value [], then getinp is used to ask its value
%        with prompt varprompt, type (s,d,f) vartype, limits varlin, and default vardef
%
%   newval: new value added
%
%    See also: VBMG_GETDEF, GETINP.
%
ansnew=answers;
if (isfield(answers,varnam))
    if (isempty(getfield(answers,varnam)))
        ansnew=setfield(answers,varnam,getinp(varprompt,vartype,varlim,vardef));
    end
else
    ansnew=setfield(answers,varnam,vardef);
end
newval=getfield(ansnew,varnam);
return
