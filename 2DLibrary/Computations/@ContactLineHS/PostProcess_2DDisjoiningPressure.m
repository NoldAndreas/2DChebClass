function PostProcess_2DDisjoiningPressure(this)    
       
    %m       = 200;
    %y1      = y1Int(1) + (0:(m-1))'/(m-1)*(y1Int(2)-y1Int(1));        
    %Int_y1  = ([y1(2:end)-y1(1:end-1);0]/2 + [0;y1(2:end)-y1(1:end-1)]/2)';    
    
    %[x,w]   = ClenCurtFlip(200);
    %y1      = y1Int(1) + (x+1)/2*(y1Int(2)-y1Int(1));
    %Int_y1  = w*(y1Int(2)-y1Int(1))/2;
    
    %this.Int_y1 = Int_y1;
    y1  = this.y1;

    mu_sat   = this.optsPhys.mu_sat;        
    %Disjoining Pressure based on 2D density profile:    
    rho      = this.rho_eq;    
    rho_wl   = kron(ones(this.HS.Pts.N1,1),this.rho1D_wl);
    rho_wg   = kron(ones(this.HS.Pts.N1,1),this.rho1D_wg);       
    
    fB            = zeros(size(y1));   
    
    [rhoGas_eq,rhoLiq_eq,p] = BulkValues(mu_sat,this.optsPhys,[],false);
    %checkContactDensity = (p + Int_1D*(rho_ic1D.*dVAdd.dy2) )/kBT;
    
    y2Max = this.optsNum.maxComp_y2 + 10;
    
    hw = waitbar(0,'Analyzing disjoining pressure');
    for iy1 = 1:length(y1)
        y1i = y1(iy1);
        
%        [floc,wIntLoc(iy1,:)]                              = this.HS.doIntFLine([y1i y1i],[0.5 y2Max],f_loc-floc_Bulk,'CHEB');
%        [fhsStrip,wIntHS(iy1,1:this.HS.AD.Mstrip)]         = this.HS.AD.Sub_Strip.doIntFLine([y1i y1i],[0 1],f_hs(this.HS.AD.mark_id(:,1))-f_hs_Bulk,'CHEB');
%        [fhsHalfSpace,wIntHS(iy1,1+this.HS.AD.Mstrip:end)] = this.HS.AD.Sub_HalfSpace.doIntFLine([y1i y1i],[1 y2Max],f_hs(this.HS.AD.mark_id(:,2))-f_hs_Bulk,'CHEB');

%        f(iy1) = floc + fhsStrip + fhsHalfSpace + R*f_hs_Bulk; 
        
        %f2(iy1)  = doIntNormalLine(this,y2Max,y1i,f_loc-floc_Bulk,f_hs-f_hs_Bulk) ;
        pts.y1_kv = y1i; pts.y2_kv = 0.5;
        IP0 = this.HS.SubShapePtsCart(pts);

        [I,w,weights,IP,pts] = this.HS.doIntFLine([y1i y1i],[0.5 y2Max],[],'CHEB');
        [h1,dVadd_i]         = getVAdd(pts.y1_kv,pts.y2_kv,0,this.optsPhys.V1);
        fB(iy1)              = - weights*(dVadd_i.dy2.*(IP*(this.rho_eq-rho_wl))) + this.optsPhys.kBT*(IP0*(this.rho_eq - rho_wl));        
        %fB(iy1)              = - weights*(dVadd_i.dy2.*(IP*this.rho_eq)) + this.optsPhys.kBT*(IP0*this.rho_eq) - p;        
                
      %  dmuCheck(i) = -(kBT*(rho_i(1) - rho_wl(1)) - Int_1D*( (rho_i - rho_wl).*dVAdd.dy2))/(rhoLiq_eq-rhoGas_eq);
        
        %fB(iy1)  = -this.HS.doIntFLine([y1i y1i],[0.5 (this.optsNum.maxComp_y2+5)],dVAdd.dy2.*(this.rho_eq-rho_wl),'CHEB') + this.optsPhys.kBT*(IP*(this.rho_eq - rho_wl));
        %[I,w] = doIntFLine(this,y1P,y2P,f,TRAP_CHEB)
        
        %if(~isempty(this.IntMatrFex))
%            f2(iy1) = f2(iy1) + R*f_hs_Bulk;
%        end
        
               %... %this.HS.doIntFLine([y1i y1i],[0.5 y2Max],f_loc-floc_Bulk,'CHEB') +...            
               %this.HS.AD.Sub_Strip.doIntFLine([y1i y1i],[0 1],f_hs(this.HS.AD.mark_id(:,1))-f_hs_Bulk,'CHEB') +...
               %this.HS.AD.Sub_HalfSpace.doIntFLine([y1i y1i],[1 y2Max],f_hs(this.HS.AD.mark_id(:,2))-f_hs_Bulk,'CHEB') + ...
               %R*f_hs_Bulk;
           
        %filmThickness(iy1) = this.HS.doIntFLine([y1i y1i],[0.5 y2Max],rho-rhoGas_sat,'CHEB')/(rhoLiq_sat-rhoGas_sat);
        
        hw = waitbar(iy1/length(y1));
    end
    close(hw);        
        
    this.disjoiningPressure = fB;            
    
    PostProcess_FilmThickness(this);    
            
end