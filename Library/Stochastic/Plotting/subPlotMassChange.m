figure

plotTimes = optsPlot.plotTimes;
dt = plotTimes(2)-plotTimes(1);

plot(plotTimes,-cumsum(DDFTStruct(1).dynamicsResult.Subspace.accFlux)*dt,'r')
hold on
plot(plotTimes,-cumsum(DDFTStruct(2).dynamicsResult.Subspace.accFlux)*dt,'b')
plot(plotTimes,-cumsum(DDFTStruct(3).dynamicsResult.Subspace.accFlux)*dt,'g')
plot(plotTimes,-cumsum(DDFTStruct(4).dynamicsResult.Subspace.accFlux)*dt,'m')
plot(plotTimes,-cumsum(DDFTStruct(5).dynamicsResult.Subspace.accFlux)*dt,'c')