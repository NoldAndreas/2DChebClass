function str = SaveToFile(name,output,path)

    time = clock; % Gets the current time as a 6 element vector
    if(nargin <= 2)
        path = [pwd,'\'];
    end

    str = ([path,name,'_' ...
        num2str(time(1)),'_'... % Returns year as character
        num2str(time(2)),'_'... % Returns month as character
        num2str(time(3)),'_'... % Returns day as char
        num2str(time(4)),'_'... % returns hour as char..
        num2str(time(5)),'_'... %returns minute as char
        num2str(time(6)),... % returns seconds as char
        '.mat']);
    
    save(str,'output');
    disp(['Data saved in ',str]);

end