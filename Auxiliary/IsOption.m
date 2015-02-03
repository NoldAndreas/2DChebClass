function bool = IsOption(opt,field)
    bool = ~isempty(find(ismember(opt,field)));
end