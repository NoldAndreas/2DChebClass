function averagesStruct=getAveragesDDFT(opts,DDFTStruct)

    % check what happens with the 2D case
    %[rMean,vMean]=getRVmeansDDFT(DDFTStruct);

    if(isprop(DDFTStruct.shape,'N1'))
        dim = 2;
    else
        dim = 1;
    end
    
    if(dim==1)
        [rMean,fluxMean]=getRFluxMeansDDFT(DDFTStruct);
    else
        [rMean,fluxMean]=getRFluxMeansDDFT2D(DDFTStruct);
    end
    
    averagesStruct.rMean    = rMean;
    averagesStruct.fluxMean = fluxMean;

end

