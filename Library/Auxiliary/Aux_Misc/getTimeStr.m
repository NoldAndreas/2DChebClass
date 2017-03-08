function str = getTimeStr()
    time  = clock;
    str   = ([num2str(time(1)),'_'... % Returns year as character
              num2str(time(2)),'_'... % Returns month as character
              num2str(time(3)),'_'... % Returns day as char
              num2str(time(4)),'_'... % Returns hour as char..
              num2str(time(5)),'_'... % Returns minute as char          
              num2str(round(time(6))),...    % Returns second as char          
              ]);        
end
