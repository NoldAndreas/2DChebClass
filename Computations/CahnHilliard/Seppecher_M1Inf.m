function data = Seppecher_M1Inf(optsPhys,optsNum)
%************************************************************************* 
%data = Seppecher(optsPhys,optsNum)
%
%            (Eq)  0 = 
% with   (BC Wall) 0 = Cn*normal*grad(rho_s) + cos(theta)*(rho-1)*rho
%(BC outer Radius) 0 = normal*grad(rho_s);
%
%*************************************************************************   
    if(nargin == 0)        
        %Numerical Parameters    
        Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N1',50,...
                           'y2Min',0,'y2Max',20,'N2',20,'L1',12); %y1Max = 20

        Plot_Area = struct('y1Min',-20,'y1Max',20,'N1',120,...
                           'y2Min',0,'y2Max',Phys_Area.y2Max,'N2',40,...
                           'N1Vecs',40,'N2Vecs',6,'Hy2',3);

        Sub_Area  = struct('y1Min',-5,'y1Max',5,'N1',50,...
                           'y2Min',0,'y2Max',Phys_Area.y2Max,'N2',30);

        optsNum   = struct('PhysArea',Phys_Area,...
                           'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                           'DDFTCode','Seppecher_M1Inf',...
                           'plotTimes',0:0.1:5);                            

        optsPhys  = struct('theta',pi/2,'g',0,...
                           'D_A',0.,'rho_m',2,'nu',10,'Ca',0.02,... 
                           'nParticles',0,'UWall',-1);

        set(0,'defaultaxesfontsize',20);
        set(0,'defaultlinelinewidth',1.);
        set(0,'defaulttextfontsize',15.);
        %set(0,'defaultaxeswidth',1);
    end                     

    disp(['** ',optsNum.DDFTCode,' **']);
    
    %************************************************
    %***************  Initialization ****************
    %************************************************
    
    %Maps   = struct('PhysSpace',@Comp_to_Phys,'CompSpace',@Phys_to_Comp);
    %MapsIF = struct('PhysSpace',@Comp_to_PhysIF,'CompSpace',@Phys_to_CompIF);
    IF.y1  = zeros(optsNum.PhysArea.N2,1);    
    
    [nParticles,g,theta,D_A,rho_m,nu,Ca,Cn,UWall]  = LoadPhysData_Sepp();                   
    [N1,N2,PhysArea,SubArea,PlotArea,~]            = LoadNumDataM1_Sepp();

    FF = false(2*N1*N2,1); F = false(N1*N2,1); 
    TT = true(2*N1*N2,1);  T = true(N1*N2,1); 
    OO = ones(2*N1*N2,1);  O = ones(N1*N2,1);
    ZZ = zeros(2*N1*N2,1); Z = zeros(N1*N2,1);

    uwall      = [UWall*O ; 0*O];   
    InterpFunc = @M1SpectralSpectral_Interpolation;
        
    %************************************************
    %****************  Preprocess  ******************
    %************************************************       
    [PtsIFY2,DiffIF,IntIF]      = Spectral(MapsIF,N2);
    [Pts,Diff,Int,Ind,Interp]   = M1SpectralSpectral(Maps,N1,N2,PlotArea.x1Plot,PlotArea.x2Plot);
    
    [Path,InterpPath,Int_of_path,Int_SubOnFull] = SubSpace(SubArea,...
               @M1SpectralSpectral_Interpolation,Pts,Maps,'normal','cart');        

    [InterpPathUpper,Int_of_pathUpper] = Path2DVec(InterpFunc,Pts,Maps,@f_pathUpperLimit,N1*5,'normal');    
    [InterpPathLower,Int_of_pathLower] = Path2DVec(InterpFunc,Pts,Maps,@f_pathLowerLimit,N1*5,'normal');   
    
    [InterpPlotUV,PtsPlotSepp]         = GetUVInterpPlotting();%[-15 15],[1 5],40,30);        
            
    IBB       = [Ind.bound;Ind.bound];
    bulkSolve = (~Ind.right & ~Ind.left);
    [~,y2]    = Maps.PhysSpace(zeros(N2,1),Pts.x2);
         
    %BC at wall and left and right boundaries (at +/- infinity)
    uvBound  = Z;
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));               
                
    %****************************************************************
    %******** Solve for equilibrium density distribution ************
    %****************************************************************        
    mu    = zeros(N1*N2,1);               
    rho   = - tanh(Pts.y1_kv*3/4);        
    IF    = FindInterface(InterpFunc,Pts,Maps,rho);
    theta = pi/2; nParticles = 0; y10 = 0; D_B = 0; rhoInf = -1;
    
    
    theta = pi/2 + atan((IF.y1(1)-IF.y1(end))/PhysArea.y2Max);
    rhoInf = NewtonMethod(rhoInf,@f_eq_inf);

    y     = NewtonMethod([0;rho(bulkSolve)],@f_eq);                 
    rho   = GetFullRho(y(2:end));

    eps = 10^(-5);
    doPlots_SC_Path(InterpPathUpper,rho); xlim([2 100]); ylim([-1-eps,-1+eps]);
    
    for k = 1:3
        %nParticles = (pi/2 - theta)*(PhysArea.y1Max)^2;        
        for j = 1:10             
                                    
            %*** 1st Step ***
            %Solve for parameter D_B to ensure Mass Balance
            D_B  = fsolve(@GetMassInflux,D_B);     
            
            u_flow = GetSeppecherSolutionCart(Pts,UWall,D_A,D_B,theta);
            figure;
            doPlots_SC_Path(InterpPathUpper,u_flow.*(repmat(rho,2,1)+rho_m));
            xlim([-40 40]);            

            %*** 2nd step ***
            %Solve continuity and momentum eq. for chem. potential and
            %velocities
            uvBound([Ind.top;Ind.top]) = ...
                GetSeppecherSolutionCart([Pts.y1_kv(Ind.top)- y10,Pts.y2_kv(Ind.top)],UWall,D_A,D_B,theta);
            [mu,uv,A,b] = GetVelocityAndChemPot(rho);
            
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
    function [y,dy] = f_eq_inf(rhoInf)
        muInf       = mu(Ind.right);
        [~,y,dy]     = DoubleWellPotential(rhoInf,Cn,muInf(1));
        dy           = dy(2);
    end
    function [y,J] = f_eq(x)
        %solves for T*log*rho + Vext          
        mu_s        = x(1);        
        rho_s       = GetFullRho(x(2:end));
        [~,y,J]     = CahnHilliardFreeEnergy(rho_s,Cn,Diff,mu+mu_s);
        
        %Boundary Condition for the density
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - g;
        J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),Ind.normalBottom*Diff.grad];
        
        E               = eye(N1*N2);
        ETop            = E(Ind.top,:);
        topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
        y(Ind.top)      = topDirection*(Diff.grad*rho_s); %Ind.normalTop*(Diff.grad*rho_s); %TODO!!!
        J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad];
        
        y           = [Int*rho_s - nParticles;y(bulkSolve)];
        J           = [0,Int(bulkSolve);J(bulkSolve,[true;bulkSolve])];
        %y           = [Int_SubOnFull*rho_s - nParticles;y(bulkSolve)];
        %J           = [Int_SubOnFull(bulkSolve);J(bulkSolve,[false;bulkSolve])];
        %y           = y(bulkSolve);
        %J           = J(bulkSolve,[false;bulkSolve]);
        %y           = [Int_SubOnFull*rho_s - nParticles;y(bulkSolve)];
        %J           = [Int_SubOnFull(bulkSolve);J(bulkSolve,[false;bulkSolve])];
        
    end
    function [y] = f_eqFull(x)
        %solves for T*log*rho + Vext          
        mu_s        = x(1);        
        rho_s       = x(2:end);
        [~,y,J]     = CahnHilliardFreeEnergy(rho_s,Cn,Diff,mu+mu_s);
        
        %Boundary Condition for the density
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - g;
        %J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),Ind.normalBottom*Diff.grad];
        
        E               = eye(N1*N2);
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
    function [mu,uv,A,b] = GetVelocityAndChemPot(rho)
        
        A               = eye(3*N1*N2);  b = zeros(3*N1*N2,1);
        
        [Af,bf]         = ContMom_DiffuseInterfaceSingleFluid(rho,optsPhys,Diff);
        A([~Ind.bound;~IBB],:)   = Af([~Ind.bound;~IBB],:);   
        b([~Ind.bound;~IBB])     = bf([~Ind.bound;~IBB]);
        
        b([F;IBB])               = uvBound(IBB);        
         
        %******************************************************
        %BC for chemical Potential
        %
        %1. Top and Bottom Boundary:
        %Boundary condition for chemical potential      
        
        %******** Take momentum eq normal to boundary ********
        %A([Ind.top;FF],:)    = Ind.normalTop*Af([F;TT],:);
        %b([Ind.top;FF])      = Ind.normalTop*bf([F;TT],:);
        
        %A([Ind.bottom;FF],:) = Ind.normalBottom*Af([F;TT],:);
        %b([Ind.bottom;FF])   = Ind.normalBottom*bf([F;TT],:);        
        
        %*****************************************************
        %*** Take divergence of momentum eq ***
        logGradRho_T   = transposeVec(Diff.grad*log(rho + rho_m));        
        [~,ys,~]       = CahnHilliardFreeEnergy(rho,Cn,Diff);
        
        Cmu            = -Diff.Lap*diag(rho + rho_m);           
        Cuv            = -Ca*(2+nu)*Diff.Lap*logGradRho_T;
        C              = [Cmu,Cuv];               
        bBound         = - repmat(ys,2,1).*(Diff.grad*rho); 
        
        A([Ind.top|Ind.bottom;FF],:)  = C(Ind.top|Ind.bottom,:);                        
        b([Ind.top|Ind.bottom;FF])    = Diff.div(Ind.top|Ind.bottom,:)*bBound;
        %*****************************************************
               
        %2. Left boundary:
        %Set Chemical potential at -infinty to zero        
        A([Ind.left;FF],:) = 0;        
        A([Ind.left;FF],[Ind.left;FF]) = eye(sum(Ind.left)); 
        b([Ind.left;FF])               = 0;
        
        %3. Right boundary:
        OR                 = ones(sum(Ind.right),1);        
        IntPathUpLow       = (Int_of_pathUpper.Vec + Int_of_pathLower.Vec);        
        [Tt11,Tb11]        = CahnHilliard_StressTensorIJ(rho,optsPhys,Diff,1,1);
        [Tt12,Tb12]        = CahnHilliard_StressTensorIJ(rho,optsPhys,Diff,1,2);            
        TtHor  = [Tt11;Tt12];        
        TbHor = [Tb11;Tb12];        
                
        A(Ind.right,:)        = (OR*IntPathUpLow)*TtHor;
        A(Ind.right,[Ind.right;FF])   = A(Ind.right,[Ind.right;FF]) - ...
                                diag(rho(Ind.right)+rho_m)*PhysArea.y2Max;
        b(Ind.right)          = -IntPathUpLow*TbHor;        
        %******************************************************
              
        x               = A\b;
        
        disp(['Error: ',num2str(max(abs(A*x-b)))]);

        mu              = x(1:N1*N2);
        uv              = x(1+N1*N2:end);
                
        [At,bt] = CahnHilliard_DivergenceOfStressTensor(rho,optsPhys,Diff);
        disp(['Error of divergence of stress tensor: ',num2str(max(abs(At*[mu;uv] + bt)))]);       
        
        figure;
        doPlots_SC_Path(InterpPathUpper,TtHor*[mu;uv]+TbHor);
        xlim([-40 40]);
        
    end
    function m = GetMassInflux(dB)
        %Density at infinity : -1,
        %Density at -infinity: 1.
        dmCL   = ((-1+rho_m)-(1+rho_m))*(PhysArea.y2Max-0);
        
        u_flow = GetSeppecherSolutionCart(Pts,UWall,D_A,dB,theta);
        m      = Int_of_pathUpper.Vec*(u_flow.*(repmat(rho,2,1)+rho_m)) + dmCL*UWall;
    end
    %***************************************************************
    %Mapping functions:
    %***************************************************************
    function [y1_kv,y2_kv,J,dH1,dH2] = Comp_to_Phys(x1,x2)
        n  = length(x1);            
        
        [y2_kv,dy2dx2] =  LinearMap(x2,0,PhysArea.y2Max);
        [h2,dh2,ddh2]  =  Interface(y2_kv);
        [y1_kv,Diffy1] =  M1QuadMap(x1,PhysArea.L1,inf);
        y1_kv          =  y1_kv + h2;
        
        if(nargout >= 3)
            J        = zeros(n,2,2);
            J(:,1,1) = Diffy1.dydx; 
            J(:,1,2) = dh2;
            J(:,2,1) = zeros(n,1);
            J(:,2,2) = dy2dx2;            
        end
        
        if(nargout >= 4) %d^2(y1)/d(..)
            dH1        = zeros(n,2,2);
            dH1(:,1,1) = Diffy1.dyddx; 
            dH1(:,2,2) = ddh2;
        end
        
        if(nargout >= 4) %d^2(y2)/d(..)
            dH2        = zeros(n,2,2);            
        end
    end
    function [x1,x2] = Phys_to_Comp(y1,y2,phys_Area)
        if(nargin == 3)
            x2 = InvLinearMap(y2,0,phys_Area.y2Max);
            x1 = InvQuadMap(y1,phys_Area.L1,inf);
        else
            x2 = InvLinearMap(y2,0,PhysArea.y2Max);
            x1 = InvQuadMap(y1-Interface(y2),PhysArea.L1,inf);
        end
    end

    function [y2,dy2,dx,ddx,dddx,ddddx] = Comp_to_PhysIF(x2)
        [y2,dy2,dx,ddx,dddx,ddddx] =  LinearMap(x2,0,PhysArea.y2Max);        
    end
    function [x2] = Phys_to_CompIF(y2)
        x2 = InvLinearMap(y2,0,PhysArea.y2Max);
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
        IPUpdate = zeros(N1*N2);
        for ii=1:N1*N2
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
    function rho_Full = GetFullRho(rho_s)
        rho_Full             = Z;
        rho_Full(bulkSolve)  = rho_s;
        rho_Full(Ind.right)  = rhoInf;
        rho_Full(Ind.left)   = 1; 
    end

    function [nParticles,g,theta,D_A,rho_m,nu,Ca,Cn,uwall] = LoadPhysData_Sepp()
        nParticles  = optsPhys.nParticles;
        g           = optsPhys.g;
        theta       = optsPhys.theta;
                
        D_A         = optsPhys.D_A;    
        rho_m       = optsPhys.rho_m;
        nu          = optsPhys.nu;    
        Ca          = optsPhys.Ca;                 
        Cn          = 4/3;
        optsPhys.Cn = Cn;        
        uwall       = optsPhys.UWall;
    end
    function [N1,N2,PhysArea,SubArea,PlotArea,plotTimes] = LoadNumDataM1_Sepp()
        N1        = optsNum.PhysArea.N1; 
        N2        = optsNum.PhysArea.N2;
        PhysArea  = optsNum.PhysArea;
        SubArea   = optsNum.SubArea;   
        PlotArea  = optsNum.PlotArea;           
        
        CompSpace = Maps.CompSpace;
        
        CheckAreaBoundaries(PhysArea,PlotArea);
        CheckAreaBoundaries(PhysArea,SubArea);
        
        x1Min = min(CompSpace(PlotArea.y1Min,PlotArea.y2Min,PhysArea),...
                    CompSpace(PlotArea.y1Min,PlotArea.y2Max,PhysArea));
        x1Max = max(CompSpace(PlotArea.y1Max,PlotArea.y2Min,PhysArea),...
                    CompSpace(PlotArea.y1Max,PlotArea.y2Max,PhysArea));                
        x2Min = min(CompSpace(PlotArea.y1Min,PlotArea.y2Min,PhysArea),...
                    CompSpace(PlotArea.y1Max,PlotArea.y2Min,PhysArea));
        x2Max = max(CompSpace(PlotArea.y1Min,PlotArea.y2Max,PhysArea),...
                    CompSpace(PlotArea.y1Max,PlotArea.y2Max,PhysArea));
        
        PlotArea.x1Plot = GetArray(x1Min,x1Max,PlotArea.N1);
        PlotArea.x2Plot = GetArray(x2Min,x2Max,PlotArea.N2);
                        
        plotTimes  = optsNum.plotTimes;        
    end

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
        aT                = zeros(N1*N2,2*N1*N2);
        aT(:,1:N1*N2)     = diag(a(1:N1*N2));
        aT(:,1+N1*N2:end) = diag(a(1+N1*N2:end));
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