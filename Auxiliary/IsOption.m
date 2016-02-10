function bool = IsOption(opt,field)
    if(ischar(opt))
        opt = {opt};
    end
    bool = ~isempty(find(ismember(opt,field)));
end