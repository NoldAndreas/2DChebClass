function W = HIMatrices_RP(optsPhys,domain)

    sigmaHS = optsPhys.HI.sigmaHS;
    sigmaS  = optsPhys.V2.sigmaS;
    
    nSpecies = size(sigmaHS,1);
    
    uniqueSigmaHS = unique(sigmaHS);
    
    nSigmaHS = length(uniqueSigmaHS);

    W(nSpecies,nSpecies).W1 = 0;
    
    for iSigmaHS = 1:nSigmaHS
        optsPhys.sigmaHij = uniqueSigmaHS(iSigmaHS);
        optsPhys.HIShapeParams.RMin = optsPhys.sigmaHij;
        
        W2 = RP12(optsPhys,domain);
        W1 = zeros(size(W2));
        
        for iSpecies = 1:nSpecies
            for jSpecies = 1:nSpecies;
                if(sigmaHS(iSpecies,jSpecies)==uniqueSigmaHS(iSigmaHS))
                    W(iSpecies,jSpecies).W2 = W2;
                    W(iSpecies,jSpecies).W1 = W1;
                end
            end
        end
        
    end
    
%     for iSpecies = 1:nSpecies
%         for jSpecies = iSpecies:nSpecies
%             optsPhys.sigmaHij = sigmaHS(iSpecies,jSpecies);
%             optsPhys.HIShapeParams.RMin = optsPhys.HIShapeParams.RMinS(iSpecies,jSpecies);
%             W(iSpecies,jSpecies).W2 = RP12(optsPhys,domain);
%             W(iSpecies,jSpecies).W1 = zeros(size(W(iSpecies,jSpecies).W2));
%             W(jSpecies,iSpecies) = W(iSpecies,jSpecies);
%         end
%     end
       
end

