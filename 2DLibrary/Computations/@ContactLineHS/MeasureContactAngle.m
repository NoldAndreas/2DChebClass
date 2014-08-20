function  [CA_measured,err] = MeasureContactAngle(this,type,yInt) %,pt1,pt2]

%     theta_CS          = this.optsNum.PhysArea.alpha_deg*pi/180;
%     upperY2_angleComp = min(this.optsNum.maxComp_y2,15);        
%     rhoV              = (rhoLiq_sat+rhoGas_sat)/2; 
%     x2 = 0;
%     
%     fsolveOpts=optimset('Display','off');    
%     
%     %x2_1 = HS.CompSpace2(upperY2_angleComp/(2*sin(theta_CS)));
%     x2_1 = this.HS.CompSpace2((upperY2_angleComp-3)/(sin(theta_CS)));
%     x2   = x2_1;
%     x1_1 = fsolve(@rhoX1,-0.5,fsolveOpts);
%     x2_2 = this.HS.CompSpace2(upperY2_angleComp/sin(theta_CS));
%     x2   = x2_2;
%     x1_2 = fsolve(@rhoX1,-0.5,fsolveOpts);
% 
%     [y1_1,y2_1] = this.HS.PhysSpace(x1_1,x2_1);
%     [y1_2,y2_2] = this.HS.PhysSpace(x1_2,x2_2);
%     pt1         = this.HS.GetCartPts(y1_1,y2_1);
%     pt2         = this.HS.GetCartPts(y1_2,y2_2);
% 
%     slope  = (pt2.y2_kv-pt1.y2_kv)/(pt2.y1_kv-pt1.y1_kv);
%  %   b      = pt2.y2_kv - slope*pt2.y1_kv;
%     alphaM = mod(atan(slope),pi);
% %    y2I    = (PlotArea.y2Min:0.1:PlotArea.y2Max)';
%     fprintf(['Measured Contact Angle: ',num2str(alphaM*180/pi),' [deg]\n']);
%     
%     function z = rhoX1(x1)
%         IP = this.HS.ComputeInterpolationMatrix(x1,x2);
%         z  = IP.InterPol*rho-rhoV;
%     end    

    %Measure Contact Angle:
    
    switch type
        case 1
            %take slope of iso-interface at a certain distance from the
            %wall
            if(nargin>2)
                InitAnalysisGrid(this,[],yInt);
                ComputeInterfaceContour(this,0.1);
            end
            
            measureAtLevel = 6;
            isoInt         = this.IsolineInterface;
            [h,i]          = min(abs(isoInt - measureAtLevel));
            Diff          = barychebdiff(this.y1,1);        
            CA_measured   = atan(Diff.Dx(i,:)*isoInt)*180/pi;    
            err           = 0;
        case 2
            %average slope in a certain interval from the wall
            
            if(nargin>2)
                [this.y2,this.Int_y2,this.DiffY2] = InitAnalysisGridY(this,yInt,100);                   
                ComputeInterfaceContourY2(this,0.5);
            end
            
            D             = this.DiffY2;
            theta         = atan(1./(D*this.IsolineInterfaceY2))*180/pi;
            theta(theta<0) = theta(theta<0) + 180;
            
            CA_measured   = (this.Int_y2*theta)/sum(this.Int_y2);
            err           = std(theta);
            %err           = max(abs(theta-CA_measured));
    end
    
    disp(['Measured Contact Angle:', num2str(CA_measured),' [deg] +/- ',num2str(err),' [deg]']);    
	this.CA_deg_measured        = CA_measured;
    this.CA_deg_measured_error  = err;
    
end