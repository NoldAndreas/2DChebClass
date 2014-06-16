function data = Seppecher_M1Inf()
%************************************************************************* 
%data = Seppecher(optsPhys,optsNum)
%
%            (Eq)  0 = 
% with   (BC Wall) 0 = Cn*normal*grad(rho_s) + cos(theta)*(rho-1)*rho
%(BC outer Radius) 0 = normal*grad(rho_s);
%
%*************************************************************************   

    %% Parameters    
    PhysArea = struct('N',[100,20],'y2Min',0,'y2Max',20,'L1',12);

    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',120,...
                       'y2Min',0,'y2Max',PhysArea.y2Max,'N2',40,...
                       'N1Vecs',40,'N2Vecs',6,'Hy2',3);    

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('theta',pi/2,'g',0,...
                        'D_A',0,...
                        'zeta',10+2/3,'eta',1,...
                        'Cak',0.02,'Cn',4/3,...
                        'UWall',-1,...
                        'rho_m',2);
                    
    config = v2struct(optsPhys,optsNum);
                               
    %************************************************
    %***************  Initialization ****************
    %************************************************
    
%    M  = PhysArea.N(1)*PhysArea.N(2);
%    FF = false(2*M,1); F = false(M,1); 
%    TT = true(2*M,1);  T = true(M,1); 
%    OO = ones(2*M,1);  O = ones(M,1);
%    ZZ = zeros(2*M,1); Z = zeros(M,1);

    %uwall      = [UWall*O ; 0*O];       
        
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    
    DI = DiffuseInterface(config);
    DI.Preprocess();                  
         
    %BC at wall and left and right boundaries (at +/- infinity)
%    uvBound  = Z;
%    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));               
                
    %****************************************************************
    %******** Solve for equilibrium density distribution ************
    %****************************************************************        
    
    nParticles = 0; D_B = 0;
    
    rho    = DI.InitialGuessRho();    
    theta  = DI.FindInterfaceAngle(rho);
    rho    = DI.GetEquilibriumDensity(theta,nParticles,rho);
    
    eps = 10^(-5);    
    figure('Name','Check accuracy of map');
    DI.IC.doPlotFLine([2,100],[PhysArea.y2Max,PhysArea.y2Max],rho+1,'CART'); ylim([-eps,eps]);    
    
    for k = 1:3
        %nParticles = (pi/2 - theta)*(PhysArea.y1Max)^2;        
        for j = 1:10             
            
            %*** 1st step ***
            D_B    = DI.SetD_B(theta,rho,D_B);               
            DI.PlotSeppecherSolution(D_B,theta,rho);
            
            %*** 2nd step ***
            %Solve continuity and momentum eq. for chem. potential and
            %velocities            
            [mu,uv,A,b] = DI.GetVelocityAndChemPot(rho,D_B,theta);
            DI.PlotMu_and_U(mu,uv);
            
            %*** 3rd step ***
            %Plot intermediate step
            %subplot(2,1,1);            
                        
            figure;
            StreamlinePlotInterp(InterpPlotUV,Pts,uv); hold on;
            NormQuiverPlotInterp_NSpecies(InterpPlotUV,Pts,uv,uwall,[0 -1],1,{'b'});
            PlotPtsSepp();
            doPlots_IP_Contour(Interp,rho);  hold off;                 
            %subplot(2,1,2);
            figure;
            doPlots_SC(Interp,Pts,mu); 
            
            %*** 4th Step ***
            %Solve for density, for given chemical potential field            
            for i=1:2
                theta = pi/2 + atan((IF.y1(1)-IF.y1(end))/PhysArea.y2Max);
                %rhoInf = NewtonMethod(rhoInf,@f_eq_inf);
                y = fsolve(@f_eqFull,[0;rho]);
                
                %y     = NewtonMethod([0;rho],@f_eqFull);  
                rho   = y(2:end);
                %y     = NewtonMethod([0;rho(bulkSolve)],@f_eq);                 
                %y     = fsolve(@f_eq,rho(bulkSolve));
                %rho   = GetFullRho(y(2:end));                    

                %*** 2nd Step ***
                %Determine angle theta from density field            
                IPUpdate  = UpdateInterfaceAndMap();
                rho = IPUpdate*rho;   mu  = IPUpdate*mu;
                
                subplot(2,1,1); doPlots_IP_Contour(Interp,rho); hold on;
                plot(IF.y1,IF.y2,'linewidth',2); hold on;    PlotPtsSepp(); 
                subplot(2,1,2); doPlots_SC(Interp,Pts,rho);                            
            end

        end
    end
  
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************                            
    function [y] = f_eqFull(x)
        %solves for T*log*rho + Vext          
        mu_s        = x(1);        
        rho_s       = x(2:end);
        [~,y,J]     = CahnHilliardFreeEnergy(rho_s,Cn,Diff,mu+mu_s);
        
        %Boundary Condition for the density
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - g;
        %J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),Ind.normalBottom*Diff.grad];
        
        E               = eye(M);
        ETop            = E(Ind.top,:);
        topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
        y(Ind.top)      = topDirection*(Diff.grad*rho_s); %Ind.normalTop*(Diff.grad*rho_s); %TODO!!!
        %J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad];
        
        [~,y(Ind.right)]     = DoubleWellPotential(rho_s(Ind.right),Cn,mu_s);
        [~,y(Ind.left)]     = DoubleWellPotential(rho_s(Ind.left),Cn,mu_s);        
        
        y           = [Int_SubOnFull*rho_s - nParticles;y];
        
        
        %J           = [0,Int_SubOnFull;J];
        %y           = [Int_SubOnFull*rho_s - nParticles;y(bulkSolve)];
        %J           = [Int_SubOnFull(bulkSolve);J(bulkSolve,[false;bulkSolve])];
        %y           = y(bulkSolve);
        %J           = J(bulkSolve,[false;bulkSolve]);
        %y           = [Int_SubOnFull*rho_s - nParticles;y(bulkSolve)];
        %J           = [Int_SubOnFull(bulkSolve);J(bulkSolve,[false;bulkSolve])];
        
    end

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

    function aT = transposeVec(a)              
        aT                = zeros(M,2*M);
        aT(:,1:M)     = diag(a(1:M));
        aT(:,1+M:end) = diag(a(1+M:end));
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

%     function [w,dw,J] = DoubleWellPotential(rho,Cn,mu)
%         %rho in [-1,1]
%         %w = 1/(2*Cn^2)*(1-rho.^2).^2 + 0.5 |grad(rho)|^2 
%         w       = 1/(2*Cn^2)*(1-rho.^2).^2;
%         dw      = - 2/(Cn^2)*(1-rho.^2).*rho ;  
%         J       = - 2/(Cn^2)*diag(1-3*rho.^2);    
% 
%         if(nargin == 3) %if mu is passed as a variable
%             w      = w - mu.*rho;
%             dw     = dw - mu;  
%             J      = [-ones(length(rho),1),J];  
%         end
%     end

end
