function cols = GetGreyShades(nocols)
    cols = {}; %{'g','b','c','k','r'};      
    for iC = 1:nocols
        cols{end+1} = (nocols-iC)/nocols*[1 1 1];
    end
end