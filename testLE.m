function testLE(stocStruct,optsPhys)

    R = stocStruct.x;
    P = stocStruct.p;

    R = permute(R,[2,1,3]);  %nParticles x nSamples x nTimes
    P = permute(P,[2,1,3]);

    nBins = 25;
    optsPhys.rMin = -6;
    optsPhys.rMax = 0;

    for tPos = 101:10:101

        Rt = R(:,:,tPos);
        Pt = P(:,:,tPos);

        Rt = reshape(Rt,size(Rt,1),size(Rt,2)*size(Rt,3));
        Pt = reshape(Pt,size(Pt,1),size(Pt,2)*size(Pt,3));

        [PK,PV,binMids,nR,meanP,PbyBin] = getPressure1D(Rt,Pt,nBins,optsPhys);

        % Check local equilibrium approximation
        figure

        for iPlot = 1:nBins
            hold off
            [N,Pmid] = hist(PbyBin{iPlot},20);
            bar(Pmid,N);
            hold on
            maxN = max(N);
            plot(Pmid, maxN*exp(-(Pmid-meanP(iPlot)).^2/2) );
            title(num2str(binMids(iPlot)));
            xlim([-4,4])
            pause
        end

    end

end