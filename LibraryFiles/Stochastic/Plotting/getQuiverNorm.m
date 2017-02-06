function quiverNorm=getQuiverNorm(optsPlot)

    X1=optsPlot.YPlot1V;
    X2=optsPlot.YPlot2V;

    plotTimes=optsPlot.plotTimes;
    nPlots=length(plotTimes);
    
    quiverNorm=0;
    
    V1DV1=str2func(optsPlot.V1DV1);
    
    for iPlot=1:nPlots
    
        [VBack_S,VAdd_S]=V1DV1(X1,X2,plotTimes(iPlot),optsPlot);
        DV1S = VBack_S.dy1+VAdd_S.dy1;
        DV2S = VBack_S.dy2+VAdd_S.dy2;
        
        qNorm=sqrt(DV1S(:).^2+DV2S(:).^2);
        
        quiverNormTemp=max(qNorm); 
        
        if(quiverNormTemp>quiverNorm)
            quiverNorm=quiverNormTemp;
        end
        
    end
    
end