function SolveMovingContactLine(this,noIterations)    

    M = this.IC.M;

    if(nargin == 1)
        params.noIterations = 20;
    else
        params.noIterations = noIterations;
    end
    
    params.optsNum      = this.optsNum;
    params.optsPhys     = this.optsPhys;
    
    if(isempty(this.theta))
        otherInput.thetaInitialGuess = pi/2;
    else
        otherInput.thetaInitialGuess = this.theta;
    end
        
    IterativeSolverCahnHilliard_Full(params,otherInput);
    %[res,~,Parameters] = DataStorage('CahnHilliardSolver',@IterativeSolverCahnHilliard,params,otherInput);
        
    %this.theta = res.thetaIter(end);    
    
%     this.errors.errorIterations = res.errorIterations;
%     this.errors.errorAverage    = res.errorAverage;
%     this.errors.aIter           = res.aIter;
%     this.errors.thetaIter       = res.thetaIter;
%     
%     this.filename               = Parameters.Filename;
%     
% 	FindStagnationPoint(this);
	
    function res = IterativeSolverCahnHilliard_Full(params,otherInput)
         noIter = params.noIterations;
         
         if(otherInput.thetaInitialGuess == pi/2)
            theta  = pi/2;
            phi    = InitialGuessRho(this);
            mu     = 0;
         else
            theta  = otherInput.thetaInitialGuess;
            phi    = this.phi;
            mu     = this.GetMu();
         end
         
         this.uv  = GetBoundaryCondition(this);%,theta,phi);            
         this.p   = zeros(M,1);
         this.mu  = zeros(M,1);
         this.phi = phi;                   
         
         G    = this.mu - this.s*this.phi;             
         vec  = [this.uv;this.phi;G;this.p];
         IterationStepFullProblem(this,vec);                 
       
         figure;
         this.PlotU();
         hold on;
         this.IC.doPlots(this.phi,'contour');
         figure;
         this.IC.doPlots(this.mu,'SC');
         figure;
         this.IC.doPlots(this.p,'SC');
    end
    
	function res = IterativeSolverCahnHilliard2(params,otherInput)
        noIter = params.noIterations;
         
         if(otherInput.thetaInitialGuess == pi/2)
            theta  = pi/2;
            phi    = InitialGuessRho(this);
            mu     = 0;
         else
            theta  = otherInput.thetaInitialGuess;
            phi    = this.phi;
            mu     = this.GetMu();
         end
         uv = zeros(2*M,1);
         
         error          = zeros(noIter,1);
         errorAverage   = zeros(noIter,1);
         thetaIter      = zeros(noIter,1);
         aIter          = zeros(noIter,1);
         
         %eps = 10^(-5);        
        %this.IC.doPlotFLine([2,100],[this.optsNum.PhysArea.y2Max,...
        %   this.optsNum.PhysArea.y2Max],phi+1,'CART'); ylim([-eps,eps]);
        mu                  = SolvePhasefieldForChemPot(this,uv,phi);
        for j = 1:100
           % for j = 1:5
            [phi,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,phi);               
            mu                  = SolvePhasefieldForChemPot(this,uv,phi);
            [p,uv,~,~,a]        = GetVelocityAndChemPot(this,phi,mu,theta);
            mu                  = SolvePhasefieldForChemPot(this,uv,phi);
            %end
            
            
            DisplayFullError(this);
            
            f1 = figure;
            subplot(2,2,1)
            PlotResultsPhi(this,phi,uv,theta)
            subplot(2,2,2)
            this.IC.doPlots(p,'SC');
            subplot(2,2,3)
            this.IC.doPlots(mu,'SC');
            close(f1);
        end

        DisplayFullError(this);
    end
    
    function res = IterativeSolverCahnHilliard(params,otherInput)
        
         noIter = params.noIterations;
         
         if(otherInput.thetaInitialGuess == pi/2)
            theta  = otherInput.thetaInitialGuess;
            phi    = InitialGuessRho(this);
            mu     = 0;
         else
            theta  = otherInput.thetaInitialGuess;
            phi    = this.phi;
            mu     = this.GetMu();
         end
         uv = zeros(2*M,1);
         
         error          = zeros(noIter,1);
         errorAverage   = zeros(noIter,1);
         thetaIter      = zeros(noIter,1);
         aIter          = zeros(noIter,1);
         
         %eps = 10^(-5);        
        %this.IC.doPlotFLine([2,100],[this.optsNum.PhysArea.y2Max,...
        %   this.optsNum.PhysArea.y2Max],phi+1,'CART'); ylim([-eps,eps]);        
        [phi,theta]  = GetEquilibriumDensity_Flux(this,theta,phi,uv);
        [p,uv,~,~,a] = GetVelocityAndChemPot(this,phi,theta);
        DisplayFullError(this);
         
        %2nd Iteration
        mob = this.optsPhys.mobility; 
        for m = 100:-5:55
            this.optsPhys.mobility = m;
            [phi,theta]  = GetEquilibriumDensity_Flux(this,theta,phi,uv);
            [p,uv,~,~,a] = GetVelocityAndChemPot(this,phi,theta);
            [phi,theta]  = GetEquilibriumDensity_Flux(this,theta,phi,uv);
            [p,uv,~,~,a] = GetVelocityAndChemPot(this,phi,theta);
            DisplayFullError(this);      
        end
        
        %3rd Iteration
        for j = 1:10            
            [phi,theta]  = GetEquilibriumDensity_Flux(this,theta,phi,uv);
            [p,uv,~,~,a] = GetVelocityAndChemPot(this,phi,theta);
            DisplayFullError(this);   
        end
         
         
         res.errorIterations = error;
         res.errorAverage    = errorAverage;
         res.phi             = phi;
         res.uv              = uv;                
         res.thetaIter       = thetaIter;         
         res.aIter           = aIter;
    end

end