function varchar=tablecol2char(t,varname)
% varchar=tablecol2char(t,varname) is a utility that convers one column of
% a table to a character array.  
% The column must be singleton cells, and all of the same type.
% If the cells are non-numeric, then an empty cell (i.e., {''} is first converted to {'[]'}
%   This ensures that varchar has the same number of rows as the table.
% if the cells are numeric, then num2str is used
%
% t: a table
% varname: a variable nname
%
% varchar: a character array, size(varchar,1)=size(t,1)
%
v=t{:,varname};
if ~isnumeric(v{1}) 
    v(strmatch('',v,'exact'))={'[]'};
    varchar=strvcat(v);
else
    varchar=strvcat(num2str(cell2mat(table2array(t(:,varname)))));
end
return
