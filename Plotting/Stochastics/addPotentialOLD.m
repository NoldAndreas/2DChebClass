function quiverNorm=addPotential(hRa,ddft,optsPlot,optsPhys)

    contourWidth=get(0,'defaultlinelinewidth');

    quiverSkip=4;

    hold(hRa,'on');
    
    lineColour=optsPlot.lineColour;
    
    V1DV1=str2func(optsPhys.V1DV1);
    t=optsPlot.time;
    
    x1=ddft.Pts.y1_kv;
    x2=ddft.Pts.y2_kv;
    N1=ddft.Interp.N1;
    N2=ddft.Interp.N2;
    
    X1=reshape(x1,N2,N1);
    X2=reshape(x2,N2,N1); 
    
    X1Plot=ddft.Interp.pts1;
    X2Plot=ddft.Interp.pts1;
    N1Plot=ddft.Interp.Nplot1;
    N2Plot=ddft.Interp.Nplot2;
    
    switch optsPlot.geom

    case 'polar2D'
        % x1=r, x2=theta
        Y1 = X1.*cos(X2); 
        Y2 = X1.*sin(X2);
     
        Y1Plot = X1Plot.*cos(X2Plot); 
        Y2Plot = X1Plot.*sin(X2Plot);

    case 'planar2D'
        Y1=X1;
        Y2=X2;
        
        Y1Plot=X1Plot;
        Y2Plot=X2Plot;
        
    end
     
    qMask1=1:quiverSkip:N2;
    qMask2=1:quiverSkip:N1;
    
    nSpecies=length(optsPhys.nParticlesS);
    
    x1S=repmat(x1,1,nSpecies);
    x2S=repmat(x2,1,nSpecies);
    
    [VBack_S,VAdd_S]=V1DV1(x1S,x2S,t,optsPhys);
    
    VS=VBack_S.V+VAdd_S.V;
    

    DV1S = VBack_S.dy1+VAdd_S.dy1;
    
    DV2S = VBack_S.dy2+VAdd_S.dy2;

    quiverNorm=optsPlot.quiverNorm;
        
    y1pad=optsPlot.rMin{1}-1;
    y2pad=optsPlot.rMin{2}-1;

    Y1cut=Y1(qMask1,qMask2);
    Y2cut=Y2(qMask1,qMask2);
    Y1cut=Y1cut(:);
    Y2cut=Y2cut(:);
    
    y1Min=optsPlot.rMin{1};
    y2Min=optsPlot.rMin{2};

    y1Max=optsPlot.rMax{1};
    y2Max=optsPlot.rMax{2};
    
    maskC= (Y1cut >= y1Min) & (Y1cut <=  y1Max) & (Y2cut >= y2Min) & (Y2cut <=  y2Max);
    mask= (Y1Plot >= y1Min) & (Y1Plot <=  y1Max) & (Y2Plot >= y2Min) & (Y2Plot <=  y2Max);
      
    for iSpecies=1:nSpecies
    
        DV1=reshape(DV1S(:,iSpecies),N2,N1);
        DV2=reshape(DV2S(:,iSpecies),N2,N1);
     

        
        switch optsPlot.geom
            case 'polar2D'

                %u = ur.*cos(theta)-utheta.*sin(theta);
                %v = ur.*sin(theta)+utheta.*cos(theta);

                DV1 = DV1.*cos(X2) - DV2.*sin(X2);
                DV2 = DV1.*sin(X2) + DV2.*cos(X2); 
                
                V=reshape(fft(reshape(V,N2,N1)),N1*N2,1);
        
        end
 
        V=real(ddft.Interp.InterPol*V);


        
        DV1cut=DV1(qMask1,qMask2);
        DV2cut=DV2(qMask1,qMask2);
        
        DV1cut=DV1cut(:);
        DV2cut=DV2cut(:);
        
        if(isempty(quiverNorm))
               quiverNorm=0.1*max(max(max(abs([DV1cut;DV2cut])))); 
        end
               
        [~,h]=contour(Y1Plot(mask),Y2Plot(mask),V(mask));
        set(h,'color',lineColour{iSpecies},'linewidth',contourWidth);
        
        h=quiver([y1pad; Y1cut(maskC)],[y2pad; Y2cut(maskC)],-[quiverNorm; DV1cut(maskC)],-[0; DV2cut(maskC)]);
        set(h,'color',lineColour{iSpecies});
    end
      
    hold(hRa,'off');
    
end