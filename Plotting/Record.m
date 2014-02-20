function Record(iPlot,gifFile,delayTime)
    
    %***************************** 
    %Initialization
    if(nargin == 2)
        delayTime = 0.1;
    end
    
    hRPf = gcf;
    
    %***************************** 
    
    % plot the current frame
    f = getframe(hRPf);

    % get image data
    % need new map each frame as colours can change
    [im,map] = rgb2ind(f.cdata,256,'nodither');

    if(iPlot==1)
        %make poster file for beamer presentations
        %posterFile=[movieFile filename];
        %posterFile= filename;
        %imwrite(im,map,posterFile,'png');

        % write first frame
        imwrite(im,map,gifFile,'DelayTime',delayTime,'LoopCount',inf);
    else
        % append this frame
        [im,map] = rgb2ind(f.cdata,256,'nodither'); 
        imwrite(im,map,gifFile,'gif','WriteMode','append','DelayTime',delayTime);
    end

end