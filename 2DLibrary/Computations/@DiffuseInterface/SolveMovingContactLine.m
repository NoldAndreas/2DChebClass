function SolveMovingContactLine(this,noIterations)    

    if(nargin == 1)
        params.noIterations = 20;
    else
        params.noIterations = noIterations;
    end
    
    params.optsNum      = this.optsNum;
    params.optsPhys     = this.optsPhys;
    
    otherInput.thetaInitialGuess = 90*pi/180;

    res = DataStorage('CahnHilliardSolver',@IterativeSolverCahnHilliard,params,otherInput);
    
    this.rho   = res.rho;
    this.uv    = res.uv;
    this.theta = res.thetaIter(end);
    
    disp('*** Results ***');
    disp(['theta = ',num2str(this.theta*180/pi),' [deg]']);
    disp(['a = ',num2str(res.aIter(end))]);
    disp(['Max error of equations excluding boundaries: ',num2str(res.errorIterations(end))]);
    
    figure('color','white');
    semilogy(res.errorIterations,'ro','MarkerFaceColor','r'); hold on;
    semilogy(res.errorAverage,'ko','MarkerFaceColor','k');
    legend({'Maximal error','Average error'});

    function res = IterativeSolverCahnHilliard(params,otherInput)
        
         noIter = params.noIterations;
         
         theta  = otherInput.thetaInitialGuess;
         rho    = InitialGuessRho(this);
         mu     = 0;
         
         error          = zeros(noIter,1);
         errorAverage   = zeros(noIter,1);
         thetaIter      = zeros(noIter,1);
         aIter          = zeros(noIter,1);
         
         %eps = 10^(-5);        
        %this.IC.doPlotFLine([2,100],[this.optsNum.PhysArea.y2Max,...
        %   this.optsNum.PhysArea.y2Max],rho+1,'CART'); ylim([-eps,eps]);        
         
         for j = 1:noIter
            disp(['** Iteration no ',num2str(j),' **']) 
            close all;        
            if(mod(j,3)==2 || j > 5)
                [rho,theta]  = GetEquilibriumDensity(this,mu,theta,rho,'findTheta');
            else
                [rho,theta]  = GetEquilibriumDensity(this,mu,theta,rho);
            end
            [mu,uv,~,~,a] = GetVelocityAndChemPot(this,rho,theta);
            
      %      mu      = DoublewellPotential(rho,Cn) - Cn*(DI.IC.Diff.Lap*rho);

         %   DI.PlotResultsMu(mu,uv);
         %   DI.PlotResultsRho(uv,rho,theta);


%             figure; L_ana = 10;
%             DI.IC.doPlotFLine([-L_ana L_ana],...
%                              [PhysArea.y2Max,PhysArea.y2Max],mu,'CART'); hold on;
%             DI.IC.doPlotFLine([-L_ana L_ana],...
%                              [PhysArea.y2Max,PhysArea.y2Max]/2,mu,'CART'); hold on;
%             muM = mean(mu(DI.IC.Ind.left));   
%             muP = mean(mu(DI.IC.Ind.right));
%             plot([0 2*L_ana],[muM muM],'k:');
%             plot([0 2*L_ana],[muP muP],'k:');
%             pM = DI.GetPressure_from_ChemPotential(muM,-1);
%             pP = DI.GetPressure_from_ChemPotential(muP,1);
%             R = surfaceTension/(pM-pP);            

            [error(j),errorAverage(j)] = DisplayFullError(this,rho,uv);      
            thetaIter(j) = theta;
            aIter(j)     = a;            
            
            hold on;
            plot(j,error(j),'ro','MarkerFaceColor','r'); hold on;
            plot(j,errorAverage(j),'ko','MarkerFaceColor','k'); drawnow;           
            
         end
         
         res.errorIterations = error;
         res.errorAverage    = errorAverage;
         res.rho             = rho;
         res.uv              = uv;                
         res.thetaIter       = thetaIter;         
         res.aIter           = aIter;
    end

end