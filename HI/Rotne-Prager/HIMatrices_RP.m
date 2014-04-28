function HIStruct = HIMatrices_RP(optsPhys,domain)

    sigmaHS = optsPhys.HI.sigmaHS;
    sigmaS  = optsPhys.HI.sigmaS;
    
    sigmaPairs = [sigmaHS(:) sigmaS(:)];
    
    nSpecies = size(sigmaHS,1);
    
    uniquePairs = unique(sigmaPairs,'rows');
    
    nPairs = size(uniquePairs,1);

    HIStruct(nSpecies,nSpecies).HIInt11 = 0;
    
    for iPair = 1:nPairs
         
        optsPhys.sigmaHij = uniquePairs(iPair,1);
        optsPhys.HIShapeParams.RMin = uniquePairs(iPair,2);
        
        HITemp12 = RP12(optsPhys,domain);
        HITemp11 = zeros(size(HITemp12));
        
        for iSpecies = 1:nSpecies
            for jSpecies = 1:nSpecies;
                
                pairIJ = [sigmaHS(iSpecies,jSpecies) sigmaS(iSpecies,jSpecies)];
                   
                if(isequal(pairIJ,uniquePairs(iPair,:)))
                    HIStruct(iSpecies,jSpecies).HIInt11 = HITemp11;
                    HIStruct(iSpecies,jSpecies).HIInt12 = HITemp12;
                end
            end
        end
        
    end
     
end

