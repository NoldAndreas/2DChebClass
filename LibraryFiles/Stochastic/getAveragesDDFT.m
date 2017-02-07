function averagesStruct=getAveragesDDFT(opts,DDFTStruct)

    % check what happens with the 2D case
    %[rMean,vMean]=getRVmeansDDFT(DDFTStruct);

    if(isfield(DDFTStruct,'shape'))
        dim = 1;
    else
        dim = 2;
    end
    
    if(dim==1)
        [rMean,fluxMean,vMean]=getRFluxMeansDDFT(DDFTStruct);
    else
       % [rMean,fluxMean]=getRFluxMeansDDFT2D(DDFTStruct);
        rMean = 0;
        fluxMean = 0;
        vMean = 0;
        averagesStruct.error = 'Not implemented in 2D';
    end
    
    averagesStruct.rMean    = rMean;
    averagesStruct.fluxMean = fluxMean;
    averagesStruct.vMean = vMean;
    
end

