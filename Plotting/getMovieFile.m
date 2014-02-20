function str = getMovieFile(filename,path)

    if(ispc())
        addPath = '\Movies\';
    else
        addPath = '/Movies/';
    end
    
    if(nargin == 2 && path == true)
        str = [getResultsPath(),addPath];
    else
        time = clock; % Gets the current time as a 6 element vector
        str = ([getResultsPath(),addPath,filename,'_' ...
            num2str(time(1)),'_'... % Returns year as character
            num2str(time(2)),'_'... % Returns month as character
            num2str(time(3)),'_'... % Returns day as char
            num2str(time(4)),'_'... % returns hour as char..
            num2str(time(5)),... %returns minute as char
            '.gif']);
    end
        
end