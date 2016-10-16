function pairDensityAll(stocStruct)

    x = stocStruct.x;
    
    nTimes = size(x,3);
    nSamples = size(x,1);
    
    edges = linspace(-5,5,20)';
    opts.edges = edges;
    
    rho2 = zeros(length(edges),length(edges),nTimes);
    rhorho = rho2;
    rho = zeros(length(edges),nTimes);
    
    ht = waitbar(0,'Times');
    hs = waitbar(0,'Samples');
    
    %hf = figure;
    
    for iTime = 1:20:nTimes
        
        waitbar(iTime/nTimes,ht);
        
        xt = x(:,:,iTime);
        
        for iSample = 1:nSamples
        
            waitbar(iSample/nSamples,hs);
            
            xs = xt(iSample,:);
            xs = xs(:);
            
            [rho2ts,rhorhots,rhots] = pairDensity(xs,opts);
            
            rho2(:,:,iTime) = rho2(:,:,iTime) + rho2ts;
            rhorho(:,:,iTime) = rhorho(:,:,iTime) + rhorhots;
            rho(:,iTime) = rho(:,iTime) + rhots;
            
%             rrts = kron(rhots,rhots);
%             rrts = reshape(rrts,length(edges),length(edges));
%         
%             isequal(rrts,rhorhots)
 
%             figure(hf)
%             subplot(1,2,1)
%             bar3(rhorhots)
%             subplot(1,2,2)
%             bar3(rhorho(:,:,iTime))
           
            
        end
        
        figure
        subplot(1,3,1)
        %bar3(rho2(:,:,iTime));
        pcolor(rho2(:,:,iTime))
        subplot(1,3,2)
        %bar3(rhorho(:,:,iTime));
        pcolor(rhorho(:,:,iTime))
        subplot(1,3,3)
        pcolor(rho2(:,:,iTime)./rhorho(:,:,iTime))
        
%        rr = kron(rho(:,iTime),rho(:,iTime));
%        rr = reshape(rr,length(edges),length(edges));
        
%        isequal(rr,rhorho(:,:,iTime))
        
%        figure
%        bar3(rr)
        
%        return
        
    end
    
    

end