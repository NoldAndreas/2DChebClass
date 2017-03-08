function scatterPlot3D(stocStruct,nBins)

    if(nargin<2)
        nBins = [20,20];
    end

    x = stocStruct.xInitial;
    x = x.';
    x = x(:);
    
    x1 = x(1:3:end);
    x2 = x(2:3:end);
    x3 = x(3:3:end);
    
    %figure
    %scatter3(x1,x2,x3);

    figure
    hist3([x1,x2],nBins);
    
    figure
    hist3([x1,x3],nBins);
    
    figure
    hist3([x2,x3],nBins);



end
