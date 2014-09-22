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
        
    IterativeSolverCahnHilliard2(params,otherInput);
    %[res,~,Parameters] = DataStorage('CahnHilliardSolver',@IterativeSolverCahnHilliard,params,otherInput);
        
    this.phi   = res.phi;
    this.uv    = res.uv;
    this.theta = res.thetaIter(end);    
    
    this.errors.errorIterations = res.errorIterations;
    this.errors.errorAverage    = res.errorAverage;
    this.errors.aIter           = res.aIter;
    this.errors.thetaIter       = res.thetaIter;
    
    this.filename               = Parameters.Filename;
    
	FindStagnationPoint(this);
    
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
        
        [phi,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,phi);   
        [p,uv,~,~,a]        = GetVelocityAndChemPot(this,phi,theta);
        mu                  = SolvePhasefieldForChemPot(this,uv,phi);
        for j = 1:10
            for j = 1:5
                [phi,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,phi);   
                [p,uv,~,~,a]        = GetVelocityAndChemPot(this,phi,theta);
            end
            mu                  = SolvePhasefieldForChemPot(this,uv,phi);
            
            subplot(2,2,1)
            PlotResultsPhi(this,phi,uv,theta)
            subplot(2,2,2)
            this.IC.doPlots(p,'SC');
            subplot(2,2,3)
            this.IC.doPlots(mu,'SC');
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