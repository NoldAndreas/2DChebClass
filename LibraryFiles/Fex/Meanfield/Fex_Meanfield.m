function [MeanField,J] = Fex_Meanfield(rho,IntMatrFex,kBT)
    % get number of species
    N         = size(rho,1);
    nSpecies  = size(rho,2);    
    
    MeanField = zeros(size(rho));
    J         = zeros(N*nSpecies);
    
    if(isempty(IntMatrFex)) % default value for no mean field part
                             % (see DDFT_2D/Preprocess_MeanfieldContribution)        
        return;
    end
    
    if((nSpecies > 1) || isstruct(IntMatrFex))
        for iSpecies=1:nSpecies
            for jSpecies=1:nSpecies

                MeanField(:,iSpecies) = MeanField(:,iSpecies) ...
                  + IntMatrFex(iSpecies,jSpecies).Conv*rho(:,jSpecies);
                            
                J((1+(iSpecies-1)*N):(iSpecies*N),...
                  (1+(jSpecies-1)*N):(jSpecies*N)) = IntMatrFex(iSpecies,jSpecies).Conv;
            end
        end    
        
    else
        MeanField = IntMatrFex*rho;
        J         = IntMatrFex;
    end
    
end
