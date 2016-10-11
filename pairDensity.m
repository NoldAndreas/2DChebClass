function [rho2,rhorho,rho] = pairDensity(X,opts)
    
    [X1,X2] = meshgrid(X);
    
    % prevent counting the same particle twice
    mask = eye(length(X));
    X1 = X1(:);
    X2 = X2(:);
    X1 = X1(~mask);
    X2 = X2(~mask);
        
    %histogram2(X1,X2);
    
    edges = opts.edges;
    
    nBins = length(edges);
    
    %edges = linspace(0,1,nBins);
    
    %rho2 = histcounts2(X1,X2,edges,edges);
    
    edges2 = {edges,edges};
    
    rho2 = hist3([X1 X2],'Edges',edges2);
    
    rho = histc(X,edges);
    
    %rho = histcounts(X,edges);
    
    rhorho = kron(rho,rho);
        
    rhorho = reshape(rhorho,nBins,nBins);
    
    if(isfield(opts,'doPlots') && opts.doPlots)
    
        figure
        subplot(1,2,1)
        bar3(rho2)
        subplot(1,2,2)
        bar3(rhorho)

    end

end