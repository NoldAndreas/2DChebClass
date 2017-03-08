function str = getMovieFile(filename,path)

    global dirData

    if(nargin == 2 && path == true)
        str = dirData;
    else
        time = clock; % Gets the current time as a 6 element vector
        str = ([dirData filesep filename,'_' ...
            num2str(time(1)),'_'... % Returns year as character
            num2str(time(2)),'_'... % Returns month as character
            num2str(time(3)),'_'... % Returns day as char
            num2str(time(4)),'_'... % returns hour as char..
            num2str(time(5)),'_'... %returns minute as char
            num2str(round(time(6))),... %returns second as char
            '.gif']);
    end
        
end