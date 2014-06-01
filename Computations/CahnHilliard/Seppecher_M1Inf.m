function data = Seppecher_M1Inf()
%************************************************************************* 
%data = Seppecher(optsPhys,optsNum)
%
%            (Eq)  0 = 
% with   (BC Wall) 0 = Cn*normal*grad(rho_s) + cos(theta)*(rho-1)*rho
%(BC outer Radius) 0 = normal*grad(rho_s);
%
%*************************************************************************   

    % Numerical Parameters    
    PhysArea = struct('N',[100,20],'y2Min',0,'y2Max',20,'L1',12);

    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',120,...
                       'y2Min',0,'y2Max',PhysArea.y2Max,'N2',40,...
                       'N1Vecs',40,'N2Vecs',6,'Hy2',3);    
	nParticles  = 0;
    
    % Physical Parameters
    g           = 0;
    theta       = pi/2;
    D_A         = 0.;    
    rho_m       = 2;
    nu          = 10;
    zeta        = 10 + 2/3;
    Ca          = 0.02; 
    Cn          = 4/3;    
    UWall       = -1;        
    
    %************************************************
    %***************  Initialization ****************
    %************************************************
    
    M  = PhysArea.N(1)*PhysArea.N(2);
    FF = false(2*M,1); F = false(M,1); 
    TT = true(2*M,1);  T = true(M,1); 
    OO = ones(2*M,1);  O = ones(M,1);
    ZZ = zeros(2*M,1); Z = zeros(M,1);

    uwall      = [UWall*O ; 0*O];       
        
    %************************************************
    %****************  Preprocess  ******************
    %************************************************    
    IC                        = InfCapillaryQuad(PhysArea);    
    [Pts,Diff,Int,Ind,Interp] = IC.ComputeAll(PlotArea);   
    PtsCart                   = IC.GetCartPts();
    
    IC.SetUpBorders(300);
            
    IBB       = [Ind.bound;Ind.bound];
    bulkSolve = (~Ind.right & ~Ind.left);
         
    %BC at wall and left and right boundaries (at +/- infinity)
    uvBound  = Z;
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));               
                
    %****************************************************************
    %******** Solve for equilibrium density distribution ************
    %****************************************************************        
    mu     = zeros(M,1);               
    rho    = - tanh(Pts.y1_kv/Cn);
    
    alpha  = FindInterfaceAngle();
    theta  = pi/2; nParticles = 0; y10 = 0; D_B = 0; rhoInf = -1;
        
    theta  = pi/2 + alpha;
    rhoInf = NewtonMethod(rhoInf,@f_eq_inf);

    y     = NewtonMethod([0;rho(bulkSolve)],@f_eq);                 
    rho   = GetFullRho(y(2:end));

    eps = 10^(-5);    
    figure('Name','Check accuracy of map');
    subplot(2,2,1); IC.doPlotFLine([2,100],[PhysArea.y2Max,PhysArea.y2Max],rho+1,'CART'); ylim([-eps,eps]);    
    
    for k = 1:3
        %nParticles = (pi/2 - theta)*(PhysArea.y1Max)^2;        
        for j = 1:10             
                                    
            %*** 1st Step ***
            %Solve for parameter D_B to ensure Mass Balance
            D_B  = fsolve(@GetMassInflux,D_B);                 
            
            u_flow = GetSeppecherSolutionCart(Pts,UWall,D_A,D_B,theta);
            
            subplot(2,2,2); IC.doPlotFLine([-100,100],[PhysArea.y2Max,PhysArea.y2Max],rho.*u_flow(end/2+1:end),'CART'); 
            title('Check accuracy of map');

            IC.doPlotsStreamlines(u_flow,Ind.top); %IC.doPlotsFlux(u_flow)
            
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
        [~,y,dy]    = W(rhoInf);
        y           = y - muInf(1);        
    end

    function [mu_s,J] = GetExcessChemPotential(rho_s,mu_offset)    
        [WE,dWE,ddWE]    = W(rho_s);
        mu_s             = dWE - Cn*(Diff.Lap*rho_s) - mu_offset;                                   
        J                = diag(ddWE) - Cn*Diff.Lap;
        J                = [-ones(length(rho),1),J];  
    end

    function [z,dz,ddz] = W(rho)        
        z   = (1-rho.^2).^2/(2*Cn);
        dz  = -2*rho.*(1-rho.^2)/Cn;
        ddz = 2*(3*rho.^2-1)/Cn;
    end

    function [y,J] = f_eq(x)
        %solves for T*log*rho + Vext          
        mu_s        = x(1);        
        rho_s       = GetFullRho(x(2:end));
                
        [y,J] = GetExcessChemPotential(rho_s,mu+mu_s);        
        
        
        %% Boundary condition for the density at wall
        % 
        % $${\bf n}\cdot {\nabla \rho} = g$$
        % 
        
        
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - g;
        J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),Ind.normalBottom*Diff.grad];
        
        %% Boundary condition for the density at y2max
        % 
        % $$(\cos(\theta),\sin \theta)^T \cdot {\nabla \rho} = 0$$
        % 
        E               = eye(M);
        ETop            = E(Ind.top,:);
        topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
        y(Ind.top)      = topDirection*(Diff.grad*rho_s); %Ind.normalTop*(Diff.grad*rho_s); %TODO!!!
        J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad];
        
        %% Mass condition
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

    function [mu,uv,A,b] = GetVelocityAndChemPot(rho)
        
        A               = eye(3*M);  b = zeros(3*M,1);
        
        [Af,bf]         = ContMom_DiffuseInterfaceSingleFluid(rho);
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
        % ******** check
        logGradRho_T   = transposeVec(Diff.grad*log(rho + rho_m));        
        [~,ys,~]       = W(rho);
        ys             = ys - Cn*Diff.Lap*rho;
        
        Cmu            = -Diff.Lap*diag(rho + rho_m);           
        Cuv            = -Ca*(2+nu)*Diff.Lap*logGradRho_T;
        C              = [Cmu,Cuv];               
        bBound         = - repmat(ys,2,1).*(Diff.grad*rho); 
        
        A([Ind.top|Ind.bottom;FF],:)  = C(Ind.top|Ind.bottom,:);                        
        b([Ind.top|Ind.bottom;FF])    = Diff.div(Ind.top|Ind.bottom,:)*bBound;
        %*****************************************************
               
        %2. Left boundary:
        %Set Chemical potential at -infinty to zero        
        A([Ind.left;FF],:)             = 0;        
        A([Ind.left;FF],[Ind.left;FF]) = eye(sum(Ind.left)); 
        b([Ind.left;FF])               = 0;
        
        %3. Right boundary:
        OR                 = ones(sum(Ind.right),1);        
        IntPathUpLow       = IC.borderTop.IntNormal + IC.borderBottom.IntNormal; %?? => to check!!
%        IntPathUpLow       = (Int_of_pathUpper.Vec + Int_of_pathLower.Vec);        
        [Tt11,Tb11]        = CahnHilliard_StressTensorIJ(rho,1,1);
        [Tt12,Tb12]        = CahnHilliard_StressTensorIJ(rho,1,2);            
        TtHor              = [Tt11;Tt12];        
        TbHor              = [Tb11;Tb12];        
                
        A(Ind.right,:)        = (OR*IntPathUpLow)*TtHor;
        A(Ind.right,[Ind.right;FF])   = A(Ind.right,[Ind.right;FF]) - ...
                                        diag(rho(Ind.right)+rho_m)*PhysArea.y2Max;
        b(Ind.right)          = -IntPathUpLow*TbHor;        
        %******************************************************
              
        x               = A\b;
        
        disp(['Error: ',num2str(max(abs(A*x-b)))]);

        mu              = x(1:M);
        uv              = x(1+M:end);
                
        [At,bt] = CahnHilliard_DivergenceOfStressTensor(rho);
        disp(['Error of divergence of stress tensor: ',num2str(max(abs(At*[mu;uv] + bt)))]);       
        
      %  figure;
      %  doPlots_SC_Path(InterpPathUpper,TtHor*[mu;uv]+TbHor);
      %  xlim([-40 40]);
        
    end

    function m = GetMassInflux(dB)
        %Density at infinity : -1,
        %Density at -infinity: 1.
        massFlux   = ((-1+rho_m)-(1+rho_m))*(PhysArea.y2Max-0)*UWall;      %due to mapping to infinity
        u_flow     = GetSeppecherSolutionCart(Pts,UWall,D_A,dB,theta);
        m          = IC.borderTop.IntNormal*(u_flow.*(repmat(rho,2,1)+rho_m)) + massFlux;
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
    function rho_Full = GetFullRho(rho_s)
        rho_Full             = Z;
        rho_Full(bulkSolve)  = rho_s;
        rho_Full(Ind.right)  = rhoInf;
        rho_Full(Ind.left)   = 1; 
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
    function alpha = FindInterfaceAngle()

        fsolveOpts   = optimset('Display','off');                
        
        pt.y2_kv  =  min(PtsCart.y2_kv);
        y1CartStart = fsolve(@rhoX1,0,fsolveOpts);
                
        pt.y2_kv  =  max(PtsCart.y2_kv);
        y1CartEnd = fsolve(@rhoX1,0,fsolveOpts);                
        
        alpha = atan((y1CartStart-y1CartEnd)/...
                        (max(PtsCart.y2_kv)-  min(PtsCart.y2_kv)));
    
            
        function z = rhoX1(y1)
            pt.y1_kv = y1;
            IP       = IC.SubShapePtsCart(pt);
            z        = IP*rho;
        end    
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

    function [A,b] = ContMom_DiffuseInterfaceSingleFluid(rho)
    %Continuitiy: div(rho*uv) = 0
    %Momentum: - rho*grad(mu) +Ca*(Lap(uv) + (zeta + 1/3)*grad(div(uv)) ) +...
    %          ... +grad(rho)*(W'(rho) - mu - Cn*Lap(rho) )
    %
    % A*[mu;uv] = b corresponds to momentum and continuity Eq. for given rho               

        rho_f          = rho + rho_m;
        rho_f2         = repmat(rho_f,2,1);
        gradRho_T      = Diff.grad*rho_f;

        aT             = zeros(M,2*M);
        aT(:,1:M)      = diag(gradRho_T(1:M));
        aT(:,1+M:end)  = diag(gradRho_T(1+M:end));

        A_cont_mu      = zeros(M);
        A_cont_uv      = aT + diag(rho_f)*Diff.div; 

        A_mom_mu       = -diag(rho_f2)*Diff.grad - [diag(Diff.Dy1*rho);diag(Diff.Dy2*rho)];
        A_mom_uv       = Ca*(Diff.LapVec + (zeta + 1/3)*Diff.gradDiv);

        %A_mom_mu       = -Diff.grad;
        %A_mom_uv       = Diff.LapVec;

        A_cont         = [A_cont_mu,A_cont_uv];
        A_mom          = [A_mom_mu, A_mom_uv];
        A              = [A_cont;A_mom];  

        b              = zeros(3*M,1);
        [~,ys,~]       = W(rho);
        ys             = ys - Cn*(Diff.Lap*rho);
                
        b(1+M:end)     = - repmat(ys,2,1).*(Diff.grad*rho); 

    end
    function [A,b] = CahnHilliard_StressTensorIJ(rho,i,j)
    % get matrices for
    % T = Ca*( grad(u) + grad(u)^T + (zeta - 2/3) div(u)*I ) +...
    %   + (W(rho) + Cn/2*|grad(rho)|^2 - mu*(rho+rho_m))*I - Cn*(grad(rho) X grad(rho))

        bDiag   = W(rho) + Cn/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);
        %bDiag   = W(rho) + 1/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);
        if(i == 1 && j == 1)
            Amu     = -diag(rho+rho_m);
            Auv     = Ca*[2*Diff.Dy1 , zeros(M)] + Ca*(zeta - 2/3)*Diff.div;     
            b       = bDiag - Cn*(Diff.Dy1*rho).*(Diff.Dy1*rho);
        elseif((i==2 && j == 1)  || (i==1 && j == 2))
            Amu     = zeros(M);
            Auv     = Ca*[Diff.Dy2 , Diff.Dy1] ;
            b       = - Cn*(Diff.Dy1*rho).*(Diff.Dy2*rho);            
        elseif(i==2 && j == 2)
            Amu     = -diag(rho+rho_m);
            Auv     = Ca*[zeros(M) , 2*Diff.Dy2] + Ca*(zeta - 2/3)*Diff.div;
            b       = bDiag - Cn*(Diff.Dy2*rho).*(Diff.Dy2*rho);                 
        end

        A = [Amu Auv];

    end
    function [A,b] = CahnHilliard_DivergenceOfStressTensor(rho)

        [A11,b11] = CahnHilliard_StressTensorIJ(rho,1,1); 
        [A12,b12] = CahnHilliard_StressTensorIJ(rho,1,2); 
        [A21,b21] = CahnHilliard_StressTensorIJ(rho,2,1); 
        [A22,b22] = CahnHilliard_StressTensorIJ(rho,2,2); 

        A1        = Diff.Dy1 * A11  + Diff.Dy2*A21;
        A2        = Diff.Dy1 * A12  + Diff.Dy2*A22;

        b1        = Diff.Dy1 * b11  + Diff.Dy2*b21;
        b2        = Diff.Dy1 * b12  + Diff.Dy2*b22;

        A   = [A1;A2];
        b   = [b1;b2];      

    end

end
