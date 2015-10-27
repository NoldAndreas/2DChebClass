function ThesisNanoscale_Fig11_ComputeDynamicContactLines()

    AddPaths('ThesisNanoscale');            
    close all;
    
%     config = ThesisNanoscale_GetStandardConfig(90,[]);
%     epw = FindEpwFromContactAngle(config,180);
%     disp(epw);    

    DoDynamicComputation(90,1.071,400); %Eq: 60 degrees
    DoDynamicComputation(90,0.594,400); %Eq: 120 degrees
    
    DoDynamicComputation(45,0.856,400); %Eq: 90 degrees
    DoDynamicComputation(60,1.27,400);  %Eq: 0  degrees
    
    DoDynamicComputation(135,0.856,400); %Eq: 90 degrees
    DoDynamicComputation(120,0.0,400); %Eq: 180 degrees
        
    
%     DoDynamicComputation(90,1.154,400); %Eq: 45 degrees
%     DoDynamicComputation(90,0.453,400); %Eq: 135 degrees
%     
%     DoDynamicComputation(45,0.856,400); %Eq: 90 degrees
%     DoDynamicComputation(60,1.154,400); %Eq: 45 degrees
%     
%     DoDynamicComputation(120,0.856,400); %Eq: 90 degrees
%     DoDynamicComputation(135,0.453,400); %Eq: 135 degrees
          
    function DoDynamicComputation(alpha_deg,epw,maxT)        
        try 
            config = ThesisNanoscale_GetStandardConfig(alpha_deg,epw,maxT);
            config.optsNum.PlotAreaCart       = struct('y1Min',-5,'y1Max',10,...
                                                   'y2Min',0.5,'y2Max',15.5,...
                                                   'N1',100,'N2',100,'NFlux',20);

            CL = ContactLineHS(config);
            CL.Preprocess(); 
            
            %     %**********************************************
            %     % Equilibration from off-equilibrium IC
            %     %**********************************************
            rhoGas = CL.optsPhys.rhoGas_sat;
            rhoLiq = CL.optsPhys.rhoLiq_sat;

            [om,rho1D_wl,params] = CL.Compute1D('WL');            
            [om,rho1D_wg,params] = CL.Compute1D('WG');
            [om,rho1D_lg,params] = CL.Compute1D('LG');

            rho1D_wl = repmat(rho1D_wl,CL.IDC.N1,1);
            rho1D_wg = repmat(rho1D_wg,CL.IDC.N1,1);
            rho1D_lg = kronecker(rho1D_lg,ones(CL.IDC.N2,1));

            rho_ic = rho1D_wg + (rho1D_wl - rho1D_wg).*(rho1D_lg-rhoGas)/(rhoLiq - rhoGas);
            CL.x_eq = CL.optsPhys.kBT*log(rho_ic) + CL.Vext;    
            
            
            CL.ComputeDynamics();
            CL.PostprocessDynamics([4,5.5]);           
            CL.PlotDynamicValue({'entropy','rho_t','fittedInterface','UV_t','contactangle_0'},{'save','MovingFrameOfReference'});
        catch err
            disp('ERROR')
            rethrow(err);        
        end
    end
    

end