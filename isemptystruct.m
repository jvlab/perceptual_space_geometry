function tf=isemptystruct(x)
% function tf=isemptystruct(x) returns true if x is empty or a structure with no fields
% and false otherwise.  Note matlab's isempty(struct())=false.
if isempty(x)
    tf=true;
elseif isstruct(x)
    if length(fieldnames(x))==0
        tf=true;
    else
        tf=false;
    end
else
    tf=false;
end
return

    
