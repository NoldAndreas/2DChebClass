function pairDensity(X)
    [X1,X2] = meshgrid(X);
        
    mask = eye(length(X));
    

    X1 = X1(:);
    X2 = X2(:);

    X1 = X1(~mask);
    X2 = X2(~mask);
        
    %histogram2(X1,X2);
    
    nBins = 10;
    
    edges = linspace(0,1,nBins);
    
    [N2,Xedges2,Yedges2] = histcounts2(X1,X2,edges,edges)
    
    [N,edges] = histcounts(X,edges)
    
    NN = kron(N,N);
    
    size(NN)
    size(X)
    
    NN = reshape(NN,nBins-1,nBins-1);
    
    figure
    
    subplot(1,2,1)
    
    bar3(N2)
    
    subplot(1,2,2)
    
    bar3(NN)


end