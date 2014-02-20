function HI = ComputeHI(rho,Dmu,IntMatrHI)

% note that Dmu is v in the inertial case, but the computation of the term
% needed to give the HI takes exactly the same form.

    nSpecies=size(rho,2);
    
    rhoDmu = rho.*Dmu;

    HI_11=zeros(size(rho));
    HI_12=zeros(size(rho));
    
    for iSpecies=1:nSpecies
        for jSpecies=1:nSpecies
            HI_11(:,iSpecies) =HI_11(:,iSpecies) ...
              + IntMatrHI(iSpecies,jSpecies).HIInt11*rho(:,jSpecies);
            HI_12(:,iSpecies) =HI_12(:,iSpecies) ...
              + IntMatrHI(iSpecies,jSpecies).HIInt12*rhoDmu(:,jSpecies);
        end
    end    
    
    HI_11 = HI_11.*Dmu;
    
    HI = HI_11 + HI_12;
    
end
