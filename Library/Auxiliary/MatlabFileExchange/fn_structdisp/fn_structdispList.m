function [names,values] = fn_structdispList(Xname)
no = 1;
% function fn_structdisp Xname
% function fn_structdisp(X)
%---
% Recursively display the content of a structure and its sub-structures
%
% Input:
% - Xname/X     one can give as argument either the structure to display or
%               or a string (the name in the current workspace of the
%               structure to display)
%
% A few parameters can be adjusted inside the m file to determine when
% arrays and cell should be displayed completely or not

% Thomas Deneux
% Copyright 2005-2012

if ischar(Xname)
    X = evalin('caller',Xname);
else
    X = Xname;
    Xname = inputname(1);
end

if ~isstruct(X), error('argument should be a structure or the name of a structure'), end
rec_structdisp(Xname,X);

%---------------------------------
function rec_structdisp(Xname,X)
%---
    

%-- PARAMETERS (Edit this) --%

ARRAYMAXROWS = 10;
ARRAYMAXCOLS = 100;
ARRAYMAXELEMS = 100;
CELLMAXROWS = 10;
CELLMAXCOLS = 10;
CELLMAXELEMS = 30;
CELLRECURSIVE = false;

%----- PARAMETERS END -------%

% if(isnumeric(X))    
%     disp([Xname ':',num2str(X)])
% elseif(ischar(X))
%     disp([Xname ':',X]);
% elseif(isstruct(X))
%     dispStruct(X);
% else
%    disp(X);
%end

%fprintf('\b')

if isstruct(X) || isobject(X)
    F = fieldnames(X);
    nsub = length(F);
    Y = cell(1,nsub);
    subnames = cell(1,nsub);
    for i=1:nsub
        f = F{i};
        Y{i} = X.(f);
        subnames{i} = [Xname '.' f];
    end
elseif CELLRECURSIVE && iscell(X)
    nsub = numel(X);
    s = size(X);
    Y = X(:);
    subnames = cell(1,nsub);
    for i=1:nsub
        inds = s;
        globind = i-1;
        for k=1:length(s)
            inds(k) = 1+mod(globind,s(k));
            globind = floor(globind/s(k));
        end
        subnames{i} = [Xname '{' num2str(inds,'%i,')];
        subnames{i}(end) = '}';
    end
else
    return
end

for i=1:nsub
    a = Y{i};
    if isstruct(a) || isobject(a)
        if length(a)==1
            rec_structdisp(subnames{i},a);
        else
            for k=1:length(a)                
                ak = a(k);
                rec_structdisp([subnames{i} '(' num2str(k) ')'],a(k));
            end
        end
    elseif iscell(a)
        if size(a,1)<=CELLMAXROWS && size(a,2)<=CELLMAXCOLS && numel(a)<=CELLMAXELEMS
            rec_structdisp(subnames{i},a)
        end
    elseif size(a,1)<=ARRAYMAXROWS && size(a,2)<=ARRAYMAXCOLS && numel(a)<=ARRAYMAXELEMS                        
        if(isnumeric(a) && ~isempty(a))    
            names{no} = subnames{i};                   
            values{no} = num2str(a);
            no = no + 1;
            %disp([subnames{i} ':' num2str(a)])
        elseif(ischar(a))
            names{no} = subnames{i};        
            values{no} = a;
            no = no + 1;
            %disp([subnames{i} ':' a]);        
        elseif(~isempty(a))
            disp('ERROR');
        %elseif(~isempty(a))
%            disp(a);
        end        
    else
        disp(a);
%        disp([subnames{i} ':'])  disp(a)
    end
end
end

end
% 
% function dispStruct(st)
%     
%     fn = fieldnames(st);
%     
%     for j = 1:length(fn)
%         a = st.(fn{j});
%         if(isnumeric(a) && isscalar(a) && ~isempty(a))
%             disp([fn{j},' : ',num2str(a)]);
%         elseif(isnumeric(a) && ~isscalar(a)&& ~isempty(a))
%             disp([fn{j},' : ']); disp(a);
%         elseif(ischar(a))
%             disp([fn{j},' : ',a]);
%         end
%     end
