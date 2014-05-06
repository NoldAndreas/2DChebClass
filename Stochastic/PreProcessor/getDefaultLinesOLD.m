function D = getDefaultLines(D)

nParticlesS=D.nParticlesS;
stocName=D.stocName;
DDFTName=D.DDFTName;

%--------------------------------------------------------------------------

nSpecies=length(nParticlesS);
nStoc=length(stocName);
nDDFT=length(DDFTName);

purple=[0.8 0.2 1];
green=[0 0.8 0];

colours={'r','b',green,purple};
markers={'o','x','p','^'};
styles={'-','--',':','-.'};

%--------------------------------------------------------------------------

D.stocMarker=cell(1,nStoc);
D.stocColour=cell(1,nStoc);
D.stocStyle=cell(1,nStoc);
D.stocTitle=cell(1,nStoc);

iC=5;

for iStoc=1:nStoc
    
    mTemp=cell(1,nSpecies);
    cTemp=cell(1,nSpecies);
    sTemp=cell(1,nSpecies);
    tTemp=cell(1,nSpecies);
    
    iC=iC-1;
    if(iC<1)
        iC=4;
    end
     
    iM=0;
    
    for iSpecies=1:nSpecies
        iM=iM+1;
        if(iM>4)
            iM=1;
        end
        
        mTemp{iSpecies}=markers{iM};
        cTemp{iSpecies}=colours{iC};
        sTemp{iSpecies}='none';
        
        if(nSpecies==1)
            tTemp{iSpecies}=[stocName{iStoc}];
        else
            tTemp{iSpecies}=[stocName{iStoc} ' Species' num2str(iSpecies)];
        end
        
    end
    
    D.stocMarker{iStoc}=mTemp;
    D.stocColour{iStoc}=cTemp;
    D.stocStyle{iStoc}=sTemp;
    D.stocTitle{iStoc}=tTemp;
        
    
end

%--------------------------------------------------------------------------

D.DDFTMarker=cell(1,nDDFT);
D.DDFTColour=cell(1,nDDFT);
D.DDFTStyle=cell(1,nDDFT);
D.DDFTTitle=cell(1,nDDFT);

iC=0;

for iDDFT=1:nDDFT
    
    mTemp=cell(1,nSpecies);
    cTemp=cell(1,nSpecies);
    sTemp=cell(1,nSpecies);
    tTemp=cell(1,nSpecies);
    
    iC=iC+1;
    if(iC>4)
        iC=1;
    end
     
    iS=0;
    
    for iSpecies=1:nSpecies
        iS=iS+1;
        if(iS>4)
            iS=1;
        end
        
        mTemp{iSpecies}='none';
        cTemp{iSpecies}=colours{iC};
        sTemp{iSpecies}=styles{iS};
        
        if(nSpecies==1)
            tTemp{iSpecies}=[DDFTName{iDDFT}];
        else
            tTemp{iSpecies}=[DDFTName{iDDFT} ' Species ' num2str(iSpecies)];
        end
        
    end
    
    D.DDFTMarker{iDDFT}=mTemp;
    D.DDFTColour{iDDFT}=cTemp;
    D.DDFTStyle{iDDFT}=sTemp;
    D.DDFTTitle{iDDFT}=tTemp;
    
end

%--------------------------------------------------------------------------

eqS=cell(1,nSpecies);
eqM=cell(1,nSpecies);
eqC=cell(1,nSpecies);
eqT=cell(1,nSpecies);

for iSpecies=1:nSpecies
    eqS{iSpecies}=':';
    eqM{iSpecies}='none';
    eqC{iSpecies}='k';
    eqT{iSpecies}='Equilibrium';
end

D.eqStyle= { eqS };
D.eqMarker={ eqM };
D.eqColour={ eqC };
D.eqText = { eqT };