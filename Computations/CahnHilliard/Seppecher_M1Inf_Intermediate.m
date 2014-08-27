function Seppecher_M1Inf_Intermediate()

    close all;
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    
    %% Parameters    
    PhysArea = struct('N',[60,40],'y2Min',0,'y2Max',15,'L1',10,... %12,80,50
                      'NBorder',200);

    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',100,...
                      'y2Min',0,'y2Max',PhysArea.y2Max,'N2',100);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('theta',pi/2,'g',0,...
                        'D_A',0,...
                        'zeta',10+2/3,'eta',1,...
                        'Cak',0.1,'Cn',4/3,...
                        'UWall',1,...
                        'rho_m',4,...
                        'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
    
    surfaceTension = 4/3;
        
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    
    DI = DiffuseInterface(config);
    DI.Preprocess();                  
                             
    %***********************************************
    %******** Start iterative procedure ************
    %***********************************************    
    
    rho       = DI.InitialGuessRho();    
    theta     = 90*pi/180;%DI.FindInterfaceAngle(rho);        
    mu        = 0;
    
    
    %eps = 10^(-5);        
    %DI.IC.doPlotFLine([2,100],[PhysArea.y2Max,PhysArea.y2Max],rho+1,'CART'); ylim([-eps,eps]);        
    
    for j = 1:20             
        close all;        
        if(mod(j,3)==2 || j > 5)
            [rho,theta]  = DI.GetEquilibriumDensity(mu,theta,rho,'findTheta');
        else
            [rho,theta]  = DI.GetEquilibriumDensity(mu,theta,rho);
        end
        [mu,uv] = DI.GetVelocityAndChemPot(rho,0,theta);        
        
  %      mu      = DoublewellPotential(rho,Cn) - Cn*(DI.IC.Diff.Lap*rho);
                       
        DI.PlotResultsMu(mu,uv);
        DI.PlotResultsRho(uv,rho,theta);
        
        
        figure; L_ana = 10;
        DI.IC.doPlotFLine([-L_ana L_ana],...
                         [PhysArea.y2Max,PhysArea.y2Max],mu,'CART'); hold on;
        DI.IC.doPlotFLine([-L_ana L_ana],...
                         [PhysArea.y2Max,PhysArea.y2Max]/2,mu,'CART'); hold on;
        muM = mean(mu(DI.IC.Ind.left));   
        muP = mean(mu(DI.IC.Ind.right));
        plot([0 2*L_ana],[muM muM],'k:');
        plot([0 2*L_ana],[muP muP],'k:');
                
        
        % 2a
        % get pressure difference between +/- infinty
        % compute origin
        pM = DI.GetPressure_from_ChemPotential(muM,-1);
        pP = DI.GetPressure_from_ChemPotential(muP,1);
        
        R = surfaceTension/(pM-pP);
        
        ptC.y1 = PhysArea.y2Max/tan(theta) - sin(theta)*R;
        ptC.y2 = PhysArea.y2Max + cos(theta)*R;
                
        
        DI.DisplayFullError(rho,uv);        
    end
    
    
    DI.SavePlotResults(uv,rho,theta,mu);
  
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************                               
    function [h2,dh2,ddh2] = Interface(y)
        
        InterpIF   = Spectral_Interpolation(y,PtsIFY2,MapsIF);        
        
        h2         = InterpIF.InterPol*IF.y1;
        
        if(nargout > 1)
            dh2   = InterpIF.InterPol*(DiffIF.Dy*IF.y1);
            ddh2  = InterpIF.InterPol*(DiffIF.DDy*IF.y1);
        end
    end    
    function IPUpdate = UpdateInterfaceAndMap()
        %1st step: Compute new Interface
        IFOld = IF;
        IFNew = FindInterface(InterpFunc,Pts,Maps,rho);
        
        %2nd step: Compute new Points 
        IF            = IFNew;
        [y1New,y2New] = Comp_to_Phys(kron(Pts.x1,ones(N2,1)),kron(ones(N1,1),Pts.x2));
        IF            = IFOld;
        [x1New,x2New] = Phys_to_Comp(y1New,y2New);        
        
        %3rd step: Interpolate V onto new points
        IPUpdate = zeros(M);
        for ii=1:M
            hInt            = M1SpectralSpectral_Interpolation(x1New(ii),x2New(ii),Pts,Maps);
            IPUpdate(ii,:)  = hInt.InterPol;
        end        
        
        %4th step: Compute new Matrices etc
        IF = IFNew;
        [Pts,Diff,Int,Ind,Interp]   = M1SpectralSpectral(Maps,N1,N2,PlotArea.x1Plot,PlotArea.x2Plot);                  
        [PtsIFY2,DiffIF,IntIF]      = Spectral(MapsIF,N2);
    
        [Path,InterpPath,Int_of_path,Int_SubOnFull] = SubSpace(SubArea,...
               @M1SpectralSpectral_Interpolation,Pts,Maps,'normal','cart');        

        [InterpPathUpper,Int_of_pathUpper] = Path2DVec(InterpFunc,Pts,Maps,@f_pathUpperLimit,N1*10,'normal');    
        [InterpPathLower,Int_of_pathLower] = Path2DVec(InterpFunc,Pts,Maps,@f_pathLowerLimit,N1*10,'normal');
        
        %5th Step: Update other related values
        %slope      = (IF.y1(1)-IF.y1(end))/PhysArea.y2Max;
        %slope      = (-IF.y1(end))/PhysArea.y2Max;
        theta      = pi/2  - atan(DiffIF.Dy(end,:)*IF.y1);
        y10        = IF.y1(end) + cos(pi-theta)*PhysArea.y2Max;
        nParticles = IntIF*IF.y1;% (IF.y1(end)-IF.y1(1))*PhysArea.y2Max/2;
        disp(['theta:' , num2str(theta)]);
        disp(['y10:' , num2str(y10)]);
        disp(['nParticles:' , num2str(nParticles)]);
    end

    %***************************************************************
    %Auxiliary functions:
    %***************************************************************    
    function [y1 , y2 , dy1_dt , dy2_dt ] = f_pathUpperLimit(t) 
        [y1,y2,J]   = Comp_to_Phys(-t,ones(size(t)));        
        dy1_dt      = -J(:,1,1);
        dy2_dt      = -J(:,2,1);                
    end    
    function [y1 , y2 , dy1_dt , dy2_dt ] = f_pathLowerLimit(t) 
        [y1,y2,J]   = Comp_to_Phys(t,-ones(size(t)));        
        dy1_dt      = J(:,1,1);
        dy2_dt      = J(:,2,1);                
    end        
    function [InterpPlotUV,PtsPlotSepp] = GetUVInterpPlotting()%y1Int,y2Int,n1,n2)
        
        n1 = PlotArea.N1Vecs;  n2 = PlotArea.N2Vecs;
        
        y1Int(1) = PlotArea.y1Min; y1Int(2) = PlotArea.y1Max; 
        y2Int(1) = PlotArea.y2Min; y2Int(2) = PlotArea.y2Max;         
        
        y1s = y1Int(1) + (y1Int(2) - y1Int(1))*(0:n1)'/n1;
        y2s = y2Int(1) + (y2Int(2) - y2Int(1))*(0:n2)'/n2;
        [interpPts1x,~] = Phys_to_Comp(y1s,0);
        [~,interpPts2x] = Phys_to_Comp(0,y2s);    
        PtsPlot.y1_kv = kron(y1s,ones(size(y2s)));
        PtsPlot.y2_kv = kron(ones(size(y1s)),y2s);         
        PtsPlot.N1 = length(y1s);   PtsPlot.N2 = length(y2s);    
        InterpPlotUV   = M1SpectralSpectral_Interpolation(interpPts1x,interpPts2x,Pts,Maps);            
        
        
        %Initialization of Plotting Sets
        n2Sepp = ceil(2/(PlotArea.y2Max-PlotArea.y2Min)*n2);
        y1s = PlotArea.y1Min + (PlotArea.y1Max - PlotArea.y1Min)*(0:n1)'/n1;
        y2s = PlotArea.y2Max + PlotArea.Hy2*(0:n2Sepp)'/n2Sepp;                
        %y1s = (Plot_Area.y1Min:0.5:Plot_Area.y1Max)';
        %y2s = (PhysArea.y2Max:0.5:(PhysArea.y2Max+2))';    
        PtsPlotSepp.y1_kv = kron(y1s,ones(size(y2s)));
        PtsPlotSepp.y2_kv = kron(ones(size(y1s)),y2s);         
        PtsPlotSepp.N1 = length(y1s);   PtsPlotSepp.N2 = length(y2s);
        
    end
    function PlotPtsSepp()
        uBCSepp = GetSeppecherSolutionCart([PtsPlotSepp.y1_kv-y10,PtsPlotSepp.y2_kv],UWall,D_A,D_B,theta);
        %GetSeppecherSolutionCart(PtsPlotSepp,UWall,D_A,D_B,theta);

        NormQuiverPlot_NSpecies(PtsPlotSepp,uBCSepp,true(size(PtsPlotSepp.y1_kv)),uwall,[0 -1],1,{'r'}); hold on;
        plot([IF.y1(end);IF.y1(end)+cos(theta)*PlotArea.Hy2],...
             [IF.y2(end),IF.y2(end)+sin(theta)*PlotArea.Hy2],'linewidth',2,'color','k'); hold on;
        ylim([PhysArea.y2Min (PhysArea.y2Max+PlotArea.Hy2)]);             
    end    

end
