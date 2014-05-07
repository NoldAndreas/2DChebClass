function averagesStruct=getAveragesDDFT(opts,DDFTStruct)

    % check what happens with the 2D case
    %[rMean,vMean]=getRVmeansDDFT(DDFTStruct);
    
    [rMean,fluxMean]=getRFluxMeansDDFT(DDFTStruct);
    
    averagesStruct.rMean    = rMean;
    averagesStruct.fluxMean = fluxMean;

end

