function SolveMovingContactLine(this,noIterations)    

    if(nargin == 1)
        params.noIterations = 20;
    else
        params.noIterations = noIterations;
    end
    
    params.optsNum      = this.optsNum;
    params.optsPhys     = this.optsPhys;
    
    if(isempty(this.theta))
        otherInput.thetaInitialGuess = 90*pi/180;
    else
        otherInput.thetaInitialGuess = this.theta;
    end
        

    [res,~,Parameters] = DataStorage('CahnHilliardSolver',@IterativeSolverCahnHilliard,params,otherInput);
        
    this.phi   = res.phi;
    this.uv    = res.uv;
    this.theta = res.thetaIter(end);    
    
    this.errors.errorIterations = res.errorIterations;
    this.errors.errorAverage    = res.errorAverage;
    this.errors.aIter           = res.aIter;
    this.errors.thetaIter       = res.thetaIter;
    
    this.filename               = Parameters.Filename;
    
	FindStagnationPoint(this);
    
    function res = IterativeSolverCahnHilliard(params,otherInput)
        
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
         
         
         error          = zeros(noIter,1);
         errorAverage   = zeros(noIter,1);
         thetaIter      = zeros(noIter,1);
         aIter          = zeros(noIter,1);
         
         %eps = 10^(-5);        
        %this.IC.doPlotFLine([2,100],[this.optsNum.PhysArea.y2Max,...
        %   this.optsNum.PhysArea.y2Max],phi+1,'CART'); ylim([-eps,eps]);        
         
         for j = 1:noIter
            disp(['** Iteration no ',num2str(j),' **']) 
            close all;        
            if((j==2 || j > 4) && ~isempty(a))
                [phi,theta]  = GetEquilibriumDensity(this,mu,theta,phi,'findTheta');
            else
                [phi,theta]  = GetEquilibriumDensity(this,mu,theta,phi);
            end
            [mu,uv,~,~,a] = GetVelocityAndChemPot(this,phi,theta);                  

            [error(j),errorAverage(j)] = DisplayFullError(this,phi,uv);      
            
            if(~isempty(a))
                thetaIter(j) = theta;
                aIter(j)     = a;            
            end
            
            hold on;
            plot(j,error(j),'ro','MarkerFaceColor','r'); hold on;
            plot(j,errorAverage(j),'ko','MarkerFaceColor','k'); drawnow;           
            
         end
         
         res.errorIterations = error;
         res.errorAverage    = errorAverage;
         res.phi             = phi;
         res.uv              = uv;                
         res.thetaIter       = thetaIter;         
         res.aIter           = aIter;
    end

end