classdef HalfSpace_FMT < HalfSpace & ConvolutionFiniteSupport
    properties
        R
        AD %Check how to handle multiple Species!        
        %InterpFull
        y2wall
        Accuracy        
    end
    
    methods
        function this = HalfSpace_FMT(Geometry,R,accuracy)   
            
            shapeHS       = Geometry;
            shapeHS.y2Min = Geometry.y2wall+R;
            
            this@HalfSpace(shapeHS);
            
            this.R        = R;
            this.y2wall   = Geometry.y2wall;
            if(nargin == 3)
                this.Accuracy = accuracy;
            else
                this.Accuracy = 1;
            end
            
            %Geometry includes:
           % y2Min,h,L1,L2,N
            %Check how to handle multiple Species!
            %this.Boundary(1,length(R)) = InfCapillary();
            %for i = 1:length(R)
            shapeAD = struct('L1',Geometry.L1,'L2',Geometry.L2_AD,...
                            'y2Min',Geometry.y2wall,...
                            'h',Geometry.h,...
                            'N',Geometry.N,'N2bound',Geometry.N2bound);

            this.AD = HalfSpace_Composed(shapeAD);                        
            
             %There are four numberings: 
             % (here, there is only one domain,simplifying the numbering, as
             % all identifiers = 1)
             %(1) the numberings of the domains. The points with identifyer
             %     X are given by mark_id(:,X)       
             this.mark_id = true(this.M,1);
             %(2) the numbering of the domains in the first dimension. Points
             %in the first dimension belonging to the part with identifyer Y
             %are given by mark_id_1(:,Y)
             this.mark_id_1 = true(this.N1,1);
             %(3) the numbering of the domains in the second dimension. Points
             %in the second dimension belonging to the part with identifyer Z
             %are given by mark_id_2(:,Z)                      
             this.mark_id_2 = true(this.N2,1);
             %mark12 makes a link between identifyers of dimensions 1,2 and 
             % the global identifyer. Points which are in the first dimension
             % located in a domain with identifyer Y and in the second
             % direction in a domain with identifyer Z, belong to the global
             % domain with the identifyer X = mark_12(Y,Z).                 
             this.mark_12   = ones(1,1); 
        end                  
        function do1Dplot(this,val,y1)  
            
            %TODO: Do proper interpolation
            if(nargin == 2)
                y1 = inf;
            end
            
            %subPtsAAD.y2_kv           = (0.5:0.01:3.5)';
            %subPtsAAD.y1_kv           = inf*ones(size(subPtsAAD.y2_kv));           
            %IP                        = SubShapePts(subPtsAAD);
            
            y2_kv_1D                  = this.Pts.y2_kv(this.Pts.y1_kv == y1);
            y2_1D                     = this.Interp.pts2(Interp.pts1 == y1);
            Interp1D_IP               = this.Interp.InterPol(Interp.pts1 == y1,this.Pts.y1_kv == y1);
            
            %Interp1D_AAD.InterPol     = IP(:,this.Pts.y1_kv == inf);
            %Interp1D_AAD.pts1         = subPtsAAD.y1_kv; 
            %Interp1D_AAD.pts2         = subPtsAAD.y2_kv;
            
            plot(y2_kv_1D,val,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); hold on
            plot(y2_1D,Interp1D_IP.InterPol*val,'linewidth',1.5);
            xlim([min(Interp1D_IP.pts2) max(Interp1D_IP.pts2)]);
            h = xlabel('$y$');     set(h,'Interpreter','Latex'); set(h,'fontsize',25);
        	%h = ylabel('$\rho$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);                        
            %title(sel(i));    
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                                    
        end
        
        function AD = ComputeConvolutionFiniteSupport(this,area,weights,pts)
            %only tested if pts are in a shape of a HalfSpace                   
            %**********************************************************
            %INPUT:
            %   - area: an object with properties: ...
            %   - weights: a list of functions which get a structure with 
            %         two arrays [y1_kv,y2_kv] as input and returns one
            %         array of same size. Here, [y1_kv,y2_kv] are in polar
            %         coordinates, representing the radial and angular
            %         component, respectively.
            %   - pts : structure with 'y1_kv','y2_kv','y1','y2'. 
            %           (y1_kv,y2_kv) is in the cartesian coordinate
            %           system, and is a grid defined through [y1 (X) y2]            
            %
            %OUTPUT:
            %AD  - Average densities to get average densities
            %
            %AD(i,:,l)*rho  = int(rho(r_i+rd)*weights{k}(rd), 
            %            rd in area and (r_i + rd) in Halfspace)   
            %           where r_i is defined through the input pts
            %            
            %**********************************************************            
                        
            %*********************************
            %Initialization:
            fprintf('Computing interpolation for matrices for averaged densities..\n');
            tic
            
            AD  = zeros(this.AD.M,this.M,numel(weights)+1);%always include unity weight
            %*********************************
            
            markY2  = (pts.y2 < this.y2wall + 2*this.R);
            markYkv = (pts.y2_kv < this.y2wall + 2*this.R);
            
            ptsStrip.y1_kv = pts.y1_kv(markYkv);
            ptsStrip.y2_kv = pts.y2_kv(markYkv);
            ptsStrip.y2    = pts.y2(markY2);
            ptsStrip.y1    = pts.y1;
            
            ptsHS.y1_kv = pts.y1_kv(~markYkv);
            ptsHS.y2_kv = pts.y2_kv(~markYkv);
            ptsHS.y2    = pts.y2(~markY2);
            ptsHS.y1    = pts.y1;
    
            ptsy2 = ptsStrip.y2; 
            for iPts = 1:length(ptsy2)
                dataAD(iPts) = Intersect(this,area,struct('offset_y2',ptsy2(iPts)));
            end
            
            AD(markYkv,:,:)   = Conv_LinearGridX(this,ptsStrip,dataAD,weights);        
            AD(~markYkv,:,:)  = Conv_LinearGridXY(this,ptsHS,area,weights);    

            t = toc;
            disp([num2str(t),'s']);  
        end        
        
        function AD = ComputeConvolutionFiniteSupport2(this,area,weights,pts,params)
            %only tested if pts are in a shape of a HalfSpace                   
            %**********************************************************
            %INPUT:
            %   - area: an object with properties: ...
            %   - weights: a list of functions which get a structure with 
            %         two arrays [y1_kv,y2_kv] as input and returns one
            %         array of same size. Here, [y1_kv,y2_kv] are in polar
            %         coordinates, representing the radial and angular
            %         component, respectively.
            %   - pts : structure with 'y1_kv','y2_kv','y1','y2'. 
            %           (y1_kv,y2_kv) is in the cartesian coordinate
            %           system, and is a grid defined through [y1 (X) y2]            
            %
            %OUTPUT:
            %AD  - Average densities to get average densities
            %
            %AD(i,:,l)*rho  = int(rho(r_i+rd)*weights{k}(rd), 
            %            rd in area and (r_i + rd) in Halfspace)   
            %           where r_i is defined through the input pts
            %            
            %**********************************************************            
                        
            %*********************************
            %Initialization:
            fprintf('Computing interpolation for matrices for averaged densities..\n');
            tic
            
            AD  = zeros(this.M,this.M,numel(weights)+1);%always include unity weight
            %*********************************
            
            markY2  = (pts.y2 < this.y2wall + 2*this.R);
            markYkv = (pts.y2_kv < this.y2wall + 2*this.R);
            
            ptsStrip.y1_kv = pts.y1_kv(markYkv);
            ptsStrip.y2_kv = pts.y2_kv(markYkv);
            ptsStrip.y2    = pts.y2(markY2);
            ptsStrip.y1    = pts.y1;
            
            ptsHS.y1_kv = pts.y1_kv(~markYkv);
            ptsHS.y2_kv = pts.y2_kv(~markYkv);
            ptsHS.y2    = pts.y2(~markY2);
            ptsHS.y1    = pts.y1;
    
            ptsy2 = ptsStrip.y2;
            
            for iPts = 1:length(ptsy2)
                dataAD(iPts) = Intersect(this,area,struct('offset_y2',ptsy2(iPts)));
            end
            
            if(nargin==5)
                AD(markYkv,:,:)   = Conv_LinearGridX(this,ptsStrip,dataAD,weights,params);
                AD(~markYkv,:,:)  = Conv_LinearGridXY(this,ptsHS,area,weights,params);    
            else
                AD(markYkv,:,:)   = Conv_LinearGridX(this,ptsStrip,dataAD,weights);        
                AD(~markYkv,:,:)  = Conv_LinearGridXY(this,ptsHS,area,weights);    
            end
            
            t = toc;
            disp([num2str(t),'s']);  
        end        
        
        function [AD,AAD] = GetAverageDensities(this,area,weights)
            %**********************************************************
            %INPUT:
            %   - area: an object with properties: ...
            %   - weights: a list of functions which get a structure with 
            %         two arrays [y1_kv,y2_kv] as input and returns one
            %         array of same size. Here, [y1_kv,y2_kv] are in polar
            %         coordinates, representing the radial and angular
            %         component, respectively.
            %
            %OUTPUT:
            %AD  - Average densities to get average densities
            %AAD - average the average densities to compute free energy
            %
            %AD(i,:,l)*rho  = int(rho(r_i+rd)*weights{k}(rd), 
            %            rd in area and (r_i + rd) in Halfspace)   
            %           where r_i are the points defined in this.AD.Pts. 
            %           This is the grid on which the weighted densities
            %           are defined.
            %
            %AAD(j,:)*rho_AD = int(rho(r_j + rd)*weights{k}(rd),rd in area)
            %           where r_j are the points defined in this.Pts. 
            %           This is the grid on which the density
            %           is defined.
            %
            %NOTE:             
            %**********************************************************
            AD   = ComputeConvolutionFiniteSupport(this,area,weights,this.AD.Pts);
            AAD  = this.AD.ComputeConvolutionFiniteSupport(area,weights,this.Pts);
        end   
    end
    
    
end