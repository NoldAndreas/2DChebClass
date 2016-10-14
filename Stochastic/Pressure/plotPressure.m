function plotPressure(stocStruct,optsPhys)

    R = stocStruct.x;
    P = stocStruct.p;

    R = permute(R,[2,1,3]);  %nParticles x nSamples x nTimes
    P = permute(P,[2,1,3]);

    nBins = 100;
    optsPhys.rMin = -6;
    optsPhys.rMax = 0;

    for tPos = 61:20:101

        Rt = R(:,:,tPos);
        Pt = P(:,:,tPos);

        Rt = reshape(Rt,size(Rt,1),size(Rt,2)*size(Rt,3));
        Pt = reshape(Pt,size(Pt,1),size(Pt,2)*size(Pt,3));

        [PK,PV,binMids,nR,meanP,PbyBin] = getPressure1D(Rt,Pt,nBins,optsPhys);

        % Check local equilibrium approximation
        figure

        DmeanP = diff(meanP);
        
        subplot(1,3,1)
        hold off
        scatter(nR,PK,'r');
        hold on
        scatter(nR,PV,'b');
        title('density')
        
        subplot(1,3,2)
        hold off
        scatter(meanP,PK,'r');
        hold on
        scatter(meanP,PV,'b');
        title('velocity')
        
        subplot(1,3,3)
        hold off
        scatter(DmeanP,PK(2:end),'r');
        hold on
        scatter(DmeanP,PV(2:end),'b');
        title('Dv');
        
        figure
        scatter3(nR(2:end),DmeanP,PK(2:end),'r');
        hold on
        scatter3(nR(2:end),DmeanP,PV(2:end),'b');
        xlim([0,3])
        ylim([0,0.1])
        pause
        
    end

end