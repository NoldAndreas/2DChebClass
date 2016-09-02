function [f1,f1leq,rho,meanP] = oneBodyDistribution(X,P,opts)
    
    nBinsX = opts.nBinsX;
    nBinsP = opts.nBinsP;
    nParticles = opts.nParticles;
    nSamples = length(X)/nParticles;
    
    binEndsX=linspace(min(X),max(X),nBinsX+1);
    binWidthX = binEndsX(2)-binEndsX(1);
    binMidsX = (binEndsX(1:end-1)+binWidthX/2).';
    Xboxes = floor( nBinsX* ( X-min(X) ) / (max(X)-min(X)) ) +1;
    Xboxes(Xboxes==nBinsX+1)=nBinsX;
    
    binEndsP=linspace(min(P),max(P),nBinsP+1);
    binWidthP = binEndsP(2)-binEndsP(1);
    binMidsP = (binEndsP(1:end-1)+binWidthP/2).';
    Pboxes = floor( nBinsP* ( P-min(P) ) / (max(P)-min(P)) ) +1;
    Pboxes(Pboxes==nBinsP+1)=nBinsP;
    
    midsXP  = {binMidsX,binMidsP};
    
    f1 = hist3([X P],midsXP);
    f1 = f1'; % so X is now the x coordinate
    f1 = f1/nSamples/binWidthX/binWidthP;
    
    size(f1);
    
    rho = ( accumarray(Xboxes, ones(size(X)) ) );
    
    meanP = accumarray(Xboxes, P);    

    rhoTemp = rho;
    rhoTemp(rhoTemp==0)=1;
    
    meanP = meanP./rhoTemp;
    
    % normalize rho
    rho = rho/nSamples/binWidthX;
    
    sum(rho)*binWidthX;
    
    MB = zeros(nBinsP,nBinsX);
    for iBin = 1:nBinsX
        % normalized MB
        MBtemp = exp(-(binMidsP-meanP(iBin)).^2/2);
        N = sum(MBtemp)*binWidthP;
        MBtemp = MBtemp/N;
        %sum(MBtemp)*binWidthP;
        
        MB(:,iBin) = MBtemp';
    end
    
    
    
    rhoFull = repmat(rho',nBinsP,1);
    
    f1leq = rhoFull.*MB;
        
    sum(sum(f1leq))*binWidthX*binWidthP;
    
    sum(sum(f1))*binWidthX*binWidthP;
    
    f1diff = f1-f1leq;
    
    f1diffabs = abs(f1diff);
    
    cBottom = min( [ min(min(f1)) , min(min(f1leq)), min(min(f1diff)) ] );
    cTop = max( [ max(max(f1)) , max(max(f1leq)), max(min(f1diff)) ] );
    
    figure('Position',[0,0,1200,800]);
    subplot(1,3,1)
    pcolor(binMidsX,binMidsP,f1)
    caxis manual
    caxis([cBottom cTop]);
    colorbar
    
    subplot(1,3,2)
    pcolor(binMidsX,binMidsP,f1leq)
    caxis manual
    caxis([cBottom cTop]);
    colorbar

    subplot(1,3,3)
    pcolor(binMidsX,binMidsP,f1diffabs)
    caxis manual
    caxis([cBottom cTop]);
    colorbar

    % calculate Talyor expansion
    dv = diff(meanP);
    dv = dv(1:end-1);
    dv = repmat(dv',nBinsP,1);
    
    f1cut = f1(:,2:end-1);
    binMidsXCut = binMidsX(2:end-1);
    
    figure('Position',[0,0,1200,800]);
    subplot(1,3,1)
    pcolor(binMidsX,binMidsP,f1diff)
    colorbar
    
    subplot(1,3,2)
    pcolor(binMidsXCut,binMidsP,dv)
    colorbar
    
    subplot(1,3,3)
    pcolor(binMidsXCut,binMidsP,f1cut./dv)
    colorbar
    
    

end