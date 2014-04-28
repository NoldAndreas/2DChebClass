function W = HIMatrices_RP(optsPhys,domain)

    sigmaHS = optsPhys.HI.sigmaHS;
    sigmaS  = optsPhys.HI.sigmaS;
    
    sigmaPairs = [sigmaHS(:) sigmaS(:)];
    
    nSpecies = size(sigmaHS,1);
    
    uniquePairs = unique(sigmaPairs,'rows');
    
    nPairs = size(uniquePairs,1);

    W(nSpecies,nSpecies).W1 = 0;
    
    for iPair = 1:nPairs
         
        optsPhys.sigmaHij = uniquePairs(iPair,1);
        optsPhys.HIShapeParams.RMin = uniquePairs(iPair,2);
        
        W2 = RP12(optsPhys,domain);
        W1 = zeros(size(W2));
        
        for iSpecies = 1:nSpecies
            for jSpecies = 1:nSpecies;
                
                pairIJ = [sigmaHS(iSpecies,jSpecies) sigmaS(iSpecies,jSpecies)];
                   
                if(isequal(pairIJ,uniquePairs(iPair,:)))
                    W(iSpecies,jSpecies).W2 = W2;
                    W(iSpecies,jSpecies).W1 = W1;
                end
            end
        end
        
    end
     
end

