function [err,epsRel,pi_II] = SumRuleIIError(this,interval_y1)

     shapeSL = struct('yMin',interval_y1(1),...
                      'yMax',interval_y1(2),...
                      'N',150);                       
                  
	 y1_SL_old = this.y1_SpectralLine;
	 this.y1_SpectralLine = SpectralLine(shapeSL);
	 this.y1_SpectralLine.ComputeAll();
     
     %(1)
	 Compute_DisjoiningPressure_II(this);
     pi_II.DP = this.disjoiningPressure_II;
     pi_II.y1 = this.y1_SpectralLine;
     
     %(2)Find correct interval
     eps    = abs(this.disjoiningPressure_II(1) -this.disjoiningPressure_II(end));     
     epsRel = eps/max(abs(this.disjoiningPressure_II));  
     
     i1 = find((abs(this.disjoiningPressure_II)- 2*eps) > 0,1);
     i2 = find((abs(this.disjoiningPressure_II)- 2*eps) > 0,1,'last');
     
     if(isempty(i1) || isempty(i2))
         err    = 1;     
         disp('Error of disjoining pressure computation too large to allow prediction')
         this.y1_SpectralLine = y1_SL_old;
         return;
     end
     
     shapeSL.yMin = this.y1_SpectralLine.Pts.y(i1);
     shapeSL.yMax = this.y1_SpectralLine.Pts.y(i2);
    
     this.y1_SpectralLine = SpectralLine(shapeSL);
	 this.y1_SpectralLine.ComputeAll();

     %(3) Recompute disjoining pressure in that interval    
     Compute_DisjoiningPressure_II(this);
     
     %(4) Compute Sum rule error
     err    = SumRule_DisjoiningPressure(this,'II');                 

    % this.y1_SpectralLine = y1_SL_old;

end  