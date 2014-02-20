function data = DDFT_FMT_InfCapillary(optsPhys,optsNum,optsPlot)
%************************************************************************* 
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics:
%   (DYN i) rho_i/dt = div(rho_i*grad(mu_s_i))
%*************************************************************************   
    if(nargin == 0)
        [optsNum,optsPhys,optsPlot] = Test_DDFT_FMT_InfCapillary();
    end

    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************
        
    [kBT,nParticlesS,HS_f,eta]  = LoadPhysData_2Species(optsPhys);
    PhysArea                    = optsNum.PhysArea;   
    nSpecies                    = length(nParticlesS);
    plotTimes = optsNum.plotTimes;
    Rs        = diag(optsPhys.V2.sigmaS)/2;
    fBulk     = str2func(['FexBulk_',optsNum.FexNum.Fex]);
    
    N1        = PhysArea.N(1);   N2 = PhysArea.N(2);
        
    getFex = str2func(['Fex_',optsNum.FexNum.Fex]);
    
    if(nargin<3)
        optsPlot.lineColourDDFT={'r','b','g','k','m'};
        optsPlot.doDDFTPlots=true;
    end
    
    
    HS                        = InfCapillary_FMT(PhysArea,Rs);
    [Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(optsNum.PlotArea);
    Interp1D                  = HS.ComputeInterpolationMatrix(1,(-1:0.01:1)',true);
    
    %HS.AD.PlotGrid();
      
    subPts.y2_kv              = (0.:0.01:PhysArea.y2Max)';
    subPts.y1_kv              = inf*ones(size(subPts.y2_kv));           
    IP                        = HS.AD.SubShapePts(subPts);
    Interp1D_AD.InterPol      = IP(:,HS.AD.Pts.y1_kv == inf);
    Interp1D_AD.pts1          = subPts.y1_kv; 
    Interp1D_AD.pts2          = subPts.y2_kv;
    
	rhoBulk   = eta*6/pi;
    mu        = kBT*log(rhoBulk) + fBulk(rhoBulk,kBT);
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    
    tic
    fprintf(1,'Computing Fex matrices ...');   
    
    params         = optsPhys;
    params.FexNum  = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;
    params.Pts     = HS.Pts;     
    params.Polar   = 'cart';           
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);       
    IntMatrFex     = DataStorage('FexData',func,params,HS); %true
    IntMatrFex     = Get1DMatrices(IntMatrFex,HS);
    
    fprintf(1,'done.\n');
    tfex = toc;
    disp(tfex);        
    
    y1S     = repmat(Pts.y1_kv,1,nSpecies); 
    y2S     = repmat(Pts.y2_kv,1,nSpecies);
    [Vext,Vext_grad]  = getVBackDVBack(y1S,y2S,optsPhys.V1);       

    %I  = eye(N1*N2);
    %eyes=repmat(I,1,2);
    
    t_preprocess = toc;
  
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************        
    tic
    
    VAdd0=getVAdd(y1S,y2S,0,optsPhys.V1);
    y0 = zeros(N1*N2,nSpecies);%getInitialGuess(VAdd0);        
    y0 = y0(Pts.y1_kv==inf);
    
    PlotRosenfeldFMT_AverageDensitiesInf(HS,IntMatrFex(1),ones(size(y0)));

    %PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),ones(size(y0)));        
           
    fprintf(1,'Computing initial condition ...');
    fsolveOpts=optimset('MaxFunEvals',2000000,'MaxIter',200000);%'Display','off');   
    %[x_ic,flag]    = fsolve(@f,y0,fsolveOpts); 
    [x_ic,flag]    = fsolve(@f,y0,fsolveOpts); 
    rho_ic         = exp((x_ic-Vext(Pts.y1_kv==inf))/kBT);
    fprintf(1,'done.\n');
    
    if(flag<0)
        fprintf(1,'fsolve failed to converge\n');
        pause
    else
        fprintf(1,'Found initial equilibrium\n');
    end
    
    if(true) 
        figure;
        shift = struct('xmin',0.235,'xmax',3.31,'ymin',-1.75,'ymax',12.25,'yref',0,'xref',0.5);
        PlotBackgroundImage(['Fex' filesep 'FMT' filesep 'SpecialPlotting' filesep 'PackingFraction0_4783_RothFMTReview.gif'],shift);
            
        plot(Pts.y2_kv(Pts.y1_kv == inf),rho_ic,'o'); hold on
        plot(Interp1D.pts2,Interp1D.InterPol*[rho_ic;rho_ic;rho_ic],'linewidth',1.5);
        xlim([0.5 4.5]);    
        save2pdf([optsNum.name,'.pdf'],gcf);
        
        figure;
        PlotRosenfeldFMT_AverageDensitiesInf(HS,IntMatrFex(1),rho_ic);
        
    else
        if(optsPlot.doDDFTPlots)
            figure        
            doPlots_SC(Interp,Pts,rho_ic);
            PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),rho_ic);        
        end
    end

    t_eqSol = toc;

	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);

    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function y = f(x)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(x,0,mu);         
        y            = y(:);
    end
    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT);
        
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        
        mu_s = mu_s + x;                   
    end   
%     function mu_s = GetExcessChemPotential(x,t,mu_offset)
%         rho_s = exp((x-Vext)/kBT);                
%         mu_s  = getFex(rho_s,IntMatrFex,kBT);
%                        
%         for iSpecies=1:nSpecies
%            mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
%         end
% 
%         mu_s = mu_s + x + HS_f(rho_s,kBT) + getVAdd(y1S,y2S,t,optsPhys.V1);
%                    
%     end

    %***************************************************************
    %Auxiliary functions:
    %***************************************************************
    function [kBT,nParticles,HS_f,eta] = LoadPhysData_2Species(optsPhys)        
        kBT        = optsPhys.kBT;   
        nParticles = optsPhys.nParticlesS;     
        if(isfield(optsPhys,'HSBulk'))
            HS_f       = str2func(optsPhys.HSBulk);  
        elseif(isfield(optsNum,'HSBulk'))
            HS_f       = str2func(optsNum.HSBulk); 
        else
            HS_f       = @ZeroMap;
        end     
        eta = optsPhys.eta;
    end
    function [muSC,fnCS] = ZeroMap(h1s,h2s)
        muSC = 0;
        fnCS = 0;
    end

    function PlotRosenfeldFMT_AverageDensitiesInf(FMTShape,FMTMatrices,rho)

        figure('name','Average densities');

        nStruct = FMTMatrices.AD;    
        fields  = fieldnames(nStruct);

        noRows  = ceil(sqrt(numel(fields)));
        noCols  = ceil(numel(fields)/noRows);

        for i=1:numel(fields)
          %fields(i)
          %nStruct.(fields(i))      
            subplot(noRows,noCols,i);
            
            val = nStruct.(fields{i})*rho;
            
            plot(FMTShape.AD.Pts.y2_kv(FMTShape.AD.Pts.y1_kv == inf),val,'o'); hold on
            plot(Interp1D_AD.pts2,Interp1D_AD.InterPol*val);
            xlim([min(Interp1D_AD.pts2) max(Interp1D_AD.pts2)]);    
                        
            title(fields(i));    
        end

    %     subplot(2,2,1);
    %     FMTShape.PlotFull(FMTMatrices.AD.n2*rho,1);    
    %     title('n_2');
    %      
    % 	subplot(2,2,2);
    %     FMTShape.PlotFull(FMTMatrices.AD.n1*rho,1);
    %     title('n_1');
    %     
    %     subplot(2,2,3);
    %     FMTShape.PlotFull(FMTMatrices.AD.n1_v_1*rho,1);    
    %     title('n1_v_1');    
    %     
    %     subplot(2,2,4);
    %     FMTShape.PlotFull(FMTMatrices.AD.n1_v_2*rho,1);       
    %     title('n1_v_2');                

    end
    function y0 = getInitialGuess(VAdd)        

      if(~isfield(optsPhys,'HSBulk') && ~isfield(optsNum,'HSBulk'))
        % initial guess for mu doesn't really matter
        muInit=zeros(1,nSpecies);

        % however, require relatively good guess for rho
        % rho without interaction ie exp(-(VBack+VAdd )/kBT)
        rhoInit = exp(-(Vext+VAdd)/kBT);
        
        % inverse of normalization coefficient
        normalization = repmat( Int*rhoInit./nParticlesS' , size(rhoInit,1) ,1);

        %rho     = exp((x-Vext)/kBT) = N exp((-Vadd - Vext)/kBT
        %x/kBT - Vext/kBT = log(N) - Vadd/kBT - Vext/kBT
        %x = kBT log(N) - Vadd
        
        y0 = -kBT*log(normalization) -VAdd;
        
        y0 = [muInit; y0];
      else
          y0 = zeros(1+N1*N2,nSpecies);
      end
    
    end
    function dxdt = dx_dt(t,x)
        x       = reshape(x,N1*N2,nSpecies);
        
        mu_s     = GetExcessChemPotential(x,t,mu);
        %mu_s(Pts.y1_kv==inf,:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        h_s([Pts.y1_kv==inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions at infinity
        if(PhysArea.y1Max == Inf)        
            dxdt(Ind.outR,:)  = x(Ind.outR,:) - x_ic(Ind.outR,:);   
        else
            flux_dir        = Diff.grad*mu_s;
            dxdt(Ind.outR,:)  = Ind.normalOutR*flux_dir;             
        end
        
        
        dxdt = dxdt(:);
    end
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

end