function addPotential(hRa,optsPlot)

    NPlot1=optsPlot.NPlot1;
    NPlot2=optsPlot.NPlot2;
    
    contourWidth=get(0,'defaultlinelinewidth');

    quiverSkip=optsPlot.quiverSkip;
    % note swapped for meshgrid
    qMask1=1:quiverSkip:NPlot1;
    qMask2=1:quiverSkip:NPlot2;


    hold(hRa,'on');
    
    lineColour=optsPlot.lineColour;
    
    V1DV1=str2func(optsPlot.V1DV1);
    t=optsPlot.time;
    
    Y1=optsPlot.YPlot1;
    Y2=optsPlot.YPlot2;
    X1=optsPlot.YPlot1V;
    X2=optsPlot.YPlot2V;
            
    [VBack_S,VAdd_S]=V1DV1(X1,X2,t,optsPlot);
    
    VS=VBack_S.V+VAdd_S.V;
    DV1S = VBack_S.dy1+VAdd_S.dy1;
    DV2S = VBack_S.dy2+VAdd_S.dy2;
        
    quiverNorm=optsPlot.quiverNorm;
        
    y1pad=optsPlot.rMin{1}-1;
    y2pad=optsPlot.rMin{2}-1;

    
    Y1cut=Y1(qMask2,qMask1);
    Y2cut=Y2(qMask2,qMask1);
    Y1cut=Y1cut(:);
    Y2cut=Y2cut(:);
    
    nSpecies=size(VS,2);
    
    switch optsPlot.geom
        case 'polar2D'

            %u = ur.*cos(theta)-utheta.*sin(theta);
            %v = ur.*sin(theta)+utheta.*cos(theta);

            DV1Stemp=DV1S;
            DV2Stemp=DV2S;

            DV1S = DV1Stemp.*cos(X2) - DV2Stemp.*sin(X2);
            DV2S = DV1Stemp.*sin(X2) + DV2Stemp.*cos(X2); 
    
    end
        
    for iSpecies=1:nSpecies
    
        V=VS(:,iSpecies);
        
        V=reshape(V,NPlot2,NPlot1);
        
        DV1=DV1S(:,iSpecies);
        DV2=DV2S(:,iSpecies);
               
        DV1=reshape(DV1,NPlot2,NPlot1);
        DV2=reshape(DV2,NPlot2,NPlot1);
        
        
        DV1cut=DV1(qMask2,qMask1);
        DV2cut=DV2(qMask2,qMask1);
        
        DV1cut=DV1cut(:);
        DV2cut=DV2cut(:);
        
%         if(isempty(quiverNorm))
%                quiverNorm=0.1*max(max(max(abs([DV1(:);DV2(:)])))); 
%         end

        [~,h]=contour(Y1,Y2,V);
        set(h,'color',lineColour{iSpecies},'linewidth',contourWidth);
        
        %surf(Y1,Y2,V);
        
        h=quiver([y1pad; Y1cut],[y2pad; Y2cut],-[quiverNorm; DV1cut],-[0; DV2cut],2);
        set(h,'color',lineColour{iSpecies});
    end
          
    hold(hRa,'off');
    
    
    
end