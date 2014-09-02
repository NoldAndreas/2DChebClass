 function [uvBound_Corr,a] = GetFluidInterfaceVelocityBC(this,theta,rho)
    InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
    IntNormal_Path   = this.IC.borderTop.IntNormal_Path;   
    ptsBorderTop     = this.IC.borderTop.Pts;
    UWall            = this.optsPhys.UWall;
    rho_m            = this.optsPhys.rho_m;
    y2Max            = this.optsNum.PhysArea.y2Max;
    PtsCart          = this.IC.GetCartPts();
    Ind              = this.IC.Ind;
    rhoL             = rho(Ind.top & Ind.left);
    rhoR             = rho(Ind.top & Ind.right);

    y1Delta = GetDeltaX(this,rho,theta);

    ptsBorderTop.y1_kv  = ptsBorderTop.y1_kv - y1Delta;                        
    u_flow     = GetSeppecherSolutionCart(ptsBorderTop,UWall,0,0,theta);                                    

    rhoBorder  = InterpOntoBorder*rho;
    rhoBorder2 = repmat(rhoBorder,2,1);

    massFlux   = ((rhoR+rho_m)-(rhoL+rho_m))*(y2Max-0)*UWall;      %due to mapping to infinity

    fsolveOpts = optimset('Display','off');
    [a,~,flag] = fsolve(@massInflux,0,fsolveOpts);            
    if(flag < 1)
        cprintf('*r','CorrectVelocityProfile: Finding parameter a: No solution found.\n');                                
    else
        disp(['a = ',num2str(a)]);
    end                        

    u_flow     = GetSeppecherSolutionCart([PtsCart.y1_kv(Ind.top) - y1Delta,...
                                 PtsCart.y2_kv(Ind.top)],UWall,0,0,theta);          
    rhoBorder2 = repmat(rho(Ind.top),2,1);

    uvBound_Corr = u_flow .*(1 + a*(rhoL-rhoBorder2).^2.*(rhoR-rhoBorder2).^2);

    function m = massInflux(a)
         u_corrected = u_flow .*(1 + a*(rhoL-rhoBorder2).^2.*(rhoR-rhoBorder2).^2);
         m          = IntNormal_Path*(u_corrected.*(rhoBorder2+rho_m)) + massFlux;    
    end             

end  