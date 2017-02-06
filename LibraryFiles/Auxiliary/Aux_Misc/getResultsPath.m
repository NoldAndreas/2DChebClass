function pathResults = getResultsPath()

    global dirData   
    pathResults = dirData;
    
%     if(exist('D:\','dir'))
%         pathResults = 'D:\Results\2DCode\';
%     elseif(exist('/Users/NoldAndreas/','dir'))
%         pathResults = '/Users/NoldAndreas/Documents/Results/2DCode/';    
%     elseif(exist('/home/bgoddard/','dir'))
%         pathResults = '/home/bgoddard/work/MATLAB/Fluids/NComponent2D/Data/2DChebData';
%     elseif(exist('/Users/Ben/','dir'))
%         pathResults = '/Users/Ben/work/MATLAB/Fluids/NComponent2D/Data/2DChebData';
%     else
%         disp('Unknown computer; using current directory for save data');
%         pathResults = pwd;
%     end

end