function PlotIntermediateRegionStokes

    PhysArea = struct('N',[50,30],'y2Min',10,'y2Max',20,'L1',10,... %12,80,50
                      'NBorder',200);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',100,...
                      'y2Min',PhysArea.y2Min,'y2Max',PhysArea.y2Max,'N2',100);   
	       

    IC = InfCapillaryQuad(PhysArea); 
    IC.ComputeAll(PlotArea);        

    [u,Psi,p] =  GetHuhScriven_Solution(Cart2PolPts(IC.GetCartPts),pi/2);

    figure; IC.plotFlux(u);
    figure; IC.plotStreamlines(u,IC.Ind.finite)
    figure; IC.plot(p);
end