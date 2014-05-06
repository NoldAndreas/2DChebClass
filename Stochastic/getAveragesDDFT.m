function averagesStruct=getAveragesDDFT(opts,DDFTStruct)

    % check what happens with the 2D case
    [rMean,vMean]=getRVmeansDDFT(DDFTStruct);
    
    averagesStruct.rMean = rMean;
    averagesStruct.vMean = vMean;

end

