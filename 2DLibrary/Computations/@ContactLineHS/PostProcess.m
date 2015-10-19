function PostProcess(this,opts)

	if(isfield(this.optsNum,'PlotAreaCart'))
        shapeSL = struct('yMin',this.optsNum.PlotAreaCart.y1Min,...
                         'yMax',this.optsNum.PlotAreaCart.y1Max,...
                         'N',150);                       
        this.y1_SpectralLine = SpectralLine(shapeSL);
        this.y1_SpectralLine.ComputeAll();
    end
    
    opts.shape_y1_Line = shapeSL;
    opts.FilenameEq    = this.FilenameEq;

    res = DataStorage(['Equilibrium' filesep 'PostProcess'],@ComputePostProcess,opts,[]);%[]);

    this.AdsorptionIsotherm    = res.AdsorptionIsotherm;
    this.disjoiningPressure_II = res.disjoiningPressure_II;
    this.disjoiningPressure_IV = res.disjoiningPressure_IV;    
    
    this.hII       = res.hII;
    this.hIII      = res.hIII;
    this.hIV       = res.hIV;
    this.hContour  = res.hContour;
    
    if(isfield(res,'hI') && isfield(res,'y1_I'))
        this.hI   = res.hI;
        this.y1_I = res.y1_I;
    else
        Compute_hI(this);        
    end
         
    function res = ComputePostProcess(opts,misc)
        
         %ComputeAdsorptionIsotherm(this,'load'); %epw = 0.7: '\2DChebData\POF_FMT_ContactLine\deg90\IterativeContinuationPostProcess\2014_8_13_16_55_32.496'
         if(isfield(opts,'AdsorptionIsotherm_file'))
             ComputeAdsorptionIsotherm(this,opts.AdsorptionIsotherm_file);
         else
            ComputeAdsorptionIsotherm(this,200);
         end
        
         ix   = find(abs(this.AdsorptionIsotherm.FT) > this.optsNum.PlotAreaCart.y2Max,1);
         mark = (1:ix);
         this.AdsorptionIsotherm.FT       = this.AdsorptionIsotherm.FT(mark);
         this.AdsorptionIsotherm.mu       = this.AdsorptionIsotherm.mu(mark);
         this.AdsorptionIsotherm.rho      = this.AdsorptionIsotherm.rho(mark,:);
         this.AdsorptionIsotherm.OmEx     = this.AdsorptionIsotherm.OmEx(mark);
         this.AdsorptionIsotherm.dmuCheck = this.AdsorptionIsotherm.dmuCheck(mark);
           
         Compute_DisjoiningPressure_II(this);
         Compute_DisjoiningPressure_IV(this);

    %     %Compute height profiles    
        
         Compute_hII(this,'II');
         Compute_hIII(this);
         Compute_hI(this);
         Compute_hII(this,'IV');
         
         Compute_hContour(this,0.5);
            
         res.AdsorptionIsotherm    = this.AdsorptionIsotherm;
         res.disjoiningPressure_II = this.disjoiningPressure_II;
         res.disjoiningPressure_IV = this.disjoiningPressure_IV;
        
         
         res.hI        = this.hI;
         res.y1_I      = this.y1_I;
         res.hII       = this.hII;
         res.hIII      = this.hIII;
         res.hIV       = this.hIV;
         res.hContour  = this.hContour;        
         
    end
% 
%     InitAnalysisGrid(this,y1Int,y2Int);
%     ComputeInterfaceContour(this);
%     
%     ComputeAdsorptionIsotherm(this);	    
%     disjoingPressure1DCheck();    
%     
%     PostProcess_2DDisjoiningPressure(this);   
%     SumRule_DisjoiningPotential(this);
    
    %f = f2;
    %f = fB;            
%    this.grandPot      = f;   
%    this.grandPot2     = f2;           
    
  %  this.wIntLoc       = wIntLoc;
   % this.wIntHS        = wIntHS;
    %Compute Film thickness as contour line       
    
    
end