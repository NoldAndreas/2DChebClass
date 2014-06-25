function MeanField = Fex_Meanfield(rho,IntMatrFex,kBT)
    % get number of species

    nSpecies=size(rho,2);

    MeanField=zeros(size(rho));

    if((nSpecies > 1) || isstruct(IntMatrFex))
        for iSpecies=1:nSpecies
            for jSpecies=1:nSpecies
                MeanField(:,iSpecies) = MeanField(:,iSpecies) ...
                  + IntMatrFex(iSpecies,jSpecies).Conv*rho(:,jSpecies);
            end
        end    
    else
        MeanField = IntMatrFex*rho;
    end
end
