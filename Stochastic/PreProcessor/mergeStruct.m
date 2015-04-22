function outStruct = mergeStruct(varargin)
%outStruct = mergeStruct(varargin)
%  merges together the structures given in vargin.  Note there's no error
%  checking so make sure e.g. the structures have disjoint element names.

% get input info
structs=varargin;
nStruct=nargin;

% concatenate field names -- should really do some error checking here
% and data from each structure
fNames = [];
cells = [];
for iStruct = 1:nStruct
    
    fieldnamesi=fieldnames(structs{iStruct});
    
    if(~isempty(fieldnamesi))
        fNames = [fNames ; fieldnamesi]; %#ok        
        cells = [cells ; struct2cell(structs{iStruct})]; %#ok
    end
    
end

% convert from cells back to a structure
if(isempty(cells) && isempty(fNames))
    outStruct = struct();
else
    outStruct = cell2struct(cells, fNames, 1);
end
