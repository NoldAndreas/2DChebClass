function data = CheckHardRod()

    Phys_Area = struct('N',100,'yMin',0,'yMax',13.29);
    Plot_Area = struct('N',400,'yMin',0,'yMax',13.29);
    Fex_Num   = struct('Fex','Percus','N',100);
    
    optsNum   = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','CheckHardRod',...
                         'Tmax',7,'TN',50,...  
                         'name','default');
                     
    V1       = struct('V1DV1','gravity1D','g',0);              
    V2       = struct('V2DV2','zeroPotential1D');         

    optsPhys = struct('V1',V1,'V2',V2,...                                            
                  'kBT',1,'nParticlesS',10,'sigmaS',1); 
              
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    
    close all;  
    disp(['** ',optsNum.DDFTCode,' **']);
    
    %************************************************
    %***************  Initialization ****************
    %************************************************    
    PhysArea  = optsNum.PhysArea;
    
    aLine = SpectralLine(PhysArea);
    Int   = aLine.ComputeIntegrationVector;
    
    [Interp,yPlot] = aLine.InterpolationPlot(optsNum.PlotArea);
    InterPol = Interp.InterPol;
    
    optsFMT = Fex_Num;
    optsFMT.sigma = optsPhys.sigmaS;
    
    kBT = optsPhys.kBT;
    nSpecies = length(optsPhys.nParticlesS);
    
    IntMatrFex = aLine.ComputeFMTMatrices(optsFMT);

    getFex = str2func([Fex_Num.Fex]);

    Conv = zeros(size(IntMatrFex.F));
    
    getV = str2func([V1.V1DV1]);
    
    [VBack_S,VAdd_S] = getV(aLine.Pts.y,0,optsPhys.V1);

    VAdd = VAdd_S.V;
    
    %fsolveOpts=optimset('MaxFunEvals',2000000,'MaxIter',200000,'Display','off');    
    
    fsolveOpts = optimset('Display','off','TolFun',1e-10);
    
    rho0 = optsPhys.nParticlesS/(PhysArea.yMax-PhysArea.yMin) * ones(size(aLine.Pts.y));
    
    y0 = kBT*log(rho0);   
    
    y0mu0 = [y0;0];
    
    [xOut,h1,flag] = fsolve(@f,y0mu0,fsolveOpts);    
    
    rhoOut = exp(xOut(1:end-1)/optsPhys.kBT);
    
    Int*rhoOut;
    
    figure
    background = imread('HardRods.gif');
    
    temp = 14.29;
    xmin = - 1.85;
    xmax = temp + 0.55;
    ymin = -0.4;
    ymax = 3.05;
    
    imagesc([xmin xmax], [ymin ymax], flipud(background));
    colormap(gray);
    
    set(gca,'ydir','normal');
    hold on;
    
    plot([0, 0],[ymin,ymax],'g');
    plot([temp, temp],[ymin,ymax],'g');
    
    sigma = optsPhys.sigmaS;
    
    plot(aLine.Pts.y+sigma/2, rhoOut, 'or');
    hold on
    plot(yPlot+sigma/2, InterPol*rhoOut, '-b','LineWidth',2);
    
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************             
    
    function y = f(x)
        %solves for T*log*rho + Vext                         
        
        mu = x(end);
        x = x(1:end-1);
        
        y            = GetExcessChemPotential(x,0,mu);
        y            = y(:);
        
        rho_s = exp(x/kBT);  
        
        y = [y; optsPhys.nParticlesS - Int*rho_s];
        
    end
    
    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        
        mu_s = mu_s + x + Conv*rho_s + VAdd;
    end

end