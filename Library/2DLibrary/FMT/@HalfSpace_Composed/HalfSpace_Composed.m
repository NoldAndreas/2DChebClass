 classdef HalfSpace_Composed < ConvolutionFiniteSupport
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Always Check: The y1-Map (PhysSpace1) has to be the same for SubStrip
% and SubHalfSpace
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   properties 
       Sub_Strip        %Area_ID: 0
       Sub_HalfSpace    %Area_ID: 1
       %mark_id_1,mark_id_2,mark_id,mark_12
       Pts,N1,N2
       Mstrip,MHS,M     %no of points (in both dimensions)
       alpha           
   end
   
   methods
       function this = HalfSpace_Composed(Geometry)
           %Geometry includes:
           % y2Min,h,L1,L2,N,N2bound
           if(~isfield(Geometry,'alpha'))
               Geometry.alpha = pi/2;
           end
           
           shapeStrip      = struct('y2Min',Geometry.y2Min,...
                                    'y2Max',Geometry.y2Min + Geometry.h,...
                                    'L1',Geometry.L1,...
                                    'N',[Geometry.N(1),Geometry.N2bound],...
                                    'alpha',Geometry.alpha);
                                
           this.Sub_Strip  = InfCapillarySkewed(shapeStrip); 
           
           shapeHalfSpace  = struct('y2Min',Geometry.y2Min + Geometry.h,...
                                    'L2',Geometry.L2,'L1',Geometry.L1,...
                                    'N',Geometry.N,...
                                    'alpha',Geometry.alpha);
           this.Sub_HalfSpace = HalfSpaceSkewed(shapeHalfSpace);
           
           this.alpha       = Geometry.alpha;
           
           this.Pts.y1_kv   = [this.Sub_Strip.Pts.y1_kv;...
                               this.Sub_HalfSpace.Pts.y1_kv];
                           
           this.Pts.y2_kv   = [this.Sub_Strip.Pts.y2_kv;...
                               this.Sub_HalfSpace.Pts.y2_kv];
                           
           this.Pts.y1      = this.Pts.y1_kv(this.Pts.y2_kv == 0);
           this.Pts.y2      = this.Pts.y2_kv(this.Pts.y1_kv == inf);
           
           this.Pts.x1      = this.Sub_Strip.Pts.x1;
           
           this.N1          = Geometry.N(1);
           this.N2          = Geometry.N(2) + Geometry.N2bound;
                           
           this.Mstrip      = length(this.Sub_Strip.Pts.y1_kv);
           this.MHS         = length(this.Sub_HalfSpace.Pts.y1_kv);
           this.M           = this.Mstrip + this.MHS;
           
           
           %There are four numberings:
           %(1) the numberings of the domains. The points with identifyer
           %     X are given by mark_id(:,X)
           this.mark_id      = true(this.M,2);
           this.mark_id(:,1) = [true(this.Mstrip,1);false(this.MHS,1)];
           this.mark_id(:,2) = [false(this.Mstrip,1);true(this.MHS,1)];           
           
           %(2) the numbering of the domains in the first dimension. Points
           %in the first dimension belonging to the part with identifyer Y
           %are given by mark_id_1(:,Y)
           this.mark_id_1   = true(this.N1,1);
           
           %(3) the numbering of the domains in the second dimension. Points
           %in the second dimension belonging to the part with identifyer Z
           %are given by mark_id_2(:,Z)           
           this.mark_id_2       = true(this.N2,2);                     
           this.mark_id_2(:,1)  = [true(this.Sub_Strip.N2,1);false(this.Sub_HalfSpace.N2,1)];
           this.mark_id_2(:,2)  = [false(this.Sub_Strip.N2,1);true(this.Sub_HalfSpace.N2,1)];
           
           %mark12 makes a link between identifyers of dimensions 1,2 and 
           % the global identifyer. Points which are in the first dimension
           % located in a domain with identifyer Y and in the second
           % direction in a domain with identifyer Z, belong to the global
           % domain with the identifyer X = mark_12(Y,Z).           
           this.mark_12      = zeros(1,2);
           this.mark_12(1,1) = 1;
           this.mark_12(1,2) = 2;
       end       
       function ptsCart = GetCartPts(this,pts_y1,pts_y2)            
            if(nargin == 1)
                pts_y1 = this.Pts.y1_kv;
                pts_y2 = this.Pts.y2_kv;
            end                                              
            [ptsCart.y1_kv,ptsCart.y2_kv] = SkewedGrid(pts_y1,pts_y2,this.alpha);
       end    
       function pts = GetInvCartPts(this,ptsCart_y1,ptsCart_y2)                                              
            [pts.y1_kv,pts.y2_kv] = InvSkewedGrid(ptsCart_y1,ptsCart_y2,this.alpha);            
        end              
       function plot(this,V,opts)                      
           xl = [-5 5];
           yl = [0 5];
           N  = 50;
           
           PlotArea = struct('y1Min',xl(1),'y1Max',xl(2),...
                             'y2Min',this.Sub_Strip.y2Min,'y2Max',this.Sub_Strip.y2Max,...
                             'N1',N,'N2',N);
           this.Sub_Strip.InterpolationPlot(PlotArea,true);           
           
           PlotArea.y2Min = this.Sub_HalfSpace.y2Min;
           PlotArea.y2Max = yl(2);
           
           this.Sub_HalfSpace.InterpolationPlot(PlotArea,true);
           
           if((nargin > 2) && (ischar(opts)))
                this.Sub_Strip.plot(V(this.mark_id(:,1)),opts); hold on;           
                this.Sub_HalfSpace.plot(V(this.mark_id(:,2)),opts);
           else
               this.Sub_Strip.plot(V(this.mark_id(:,1))); hold on;           
               this.Sub_HalfSpace.plot(V(this.mark_id(:,2)));
           end
                      
           xlim(xl);
           ylim(yl);
           pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);
          
       end       
       function IP = SubShapePtsCart(this,a_shapePts)
            pts = GetInvCartPts(this,a_shapePts.y1_kv,a_shapePts.y2_kv);
            IP  = SubShapePts(this,pts);  
            %IP = InterpolationMatrix_Pointwise(this,pts.y1_kv,pts.y2_kv);                           
       end 
       function IP = SubShapePts(this,a_shapePts)
           maskStrip = (a_shapePts.y2_kv < this.Sub_Strip.y2Max);
           maskHS    = (a_shapePts.y2_kv >= this.Sub_HalfSpace.y2Min);
           
           if( sum(~(xor(maskStrip,maskHS))) ~= 0)
               exc = MException('HalfSpace_Composed:SubShapePts',...
                         'Points outside domain or points double counted');
               throw(exc);
           end           
           
           IP = zeros(length(a_shapePts.y1_kv),this.Mstrip+this.MHS);
           
           IP(maskStrip,1:this.Mstrip)  = this.Sub_Strip.InterpolationMatrix_Pointwise(...
                                            a_shapePts.y1_kv(maskStrip),...
                                            a_shapePts.y2_kv(maskStrip));
           IP(maskHS,1+this.Mstrip:end) = this.Sub_HalfSpace.InterpolationMatrix_Pointwise(...
                                            a_shapePts.y1_kv(maskHS),...
                                            a_shapePts.y2_kv(maskHS));                                                              
             
       end
        
       %X = InterpolateAndIntegratePtsOriginGrid(this,ptsOrGrid,dataDisk,weights);
       [Int,Int1,Int2] = ComputeIntegrationVector(this);
       [AAD] = ComputeConvolutionFiniteSupport(this,area,weights,pts);
       
       function [Interp1,Interp2,mark_id] = ComputeInterpolationMatrix12(this,interp1,interp2)           
           
           Interp2 = zeros(length(interp2),this.N2);           
           mark_0  = (interp2(:,2)==1);           
           
           Interp1                              = barychebevalMatrix(this.Pts.x1,interp1); 
           Interp2(mark_0,this.mark_id_2(:,1))  = barychebevalMatrix(this.Sub_Strip.Pts.x2,interp2(mark_0,1)); 
           Interp2(~mark_0,this.mark_id_2(:,2)) = barychebevalMatrix(this.Sub_HalfSpace.Pts.x2,interp2(~mark_0,1)); 
           
           mark_id = interp2(:,2);
           
       end 
       function [x1] = CompSpace1(this,y1)           
            x1 = this.Sub_Strip.CompSpace1(y1);           
       end 
	   function [x2] = CompSpace2(this,y2)
           mark_id = GetArea_ID(this,1,y2);
          
           x2(mark_id==1,1) = this.Sub_Strip.CompSpace2(y2(mark_id==1));
           x2(mark_id==2,1) = this.Sub_HalfSpace.CompSpace2(y2(mark_id==2));
           x2(:,2)          = mark_id;          
       end    
       function [I12] = Combine12(this,I1,I2)
           %I1 - linear Operator on first dimension
           %I2 - linear Operator on second dimension           
           I12 = [kronecker(I1,I2(:,this.mark_id_2{1})),kronecker(I1,I2(:,this.mark_id_2{2}))];
       end
       function id = GetArea_ID(this,y1,y2)
           id = ones(size(y2));
           id(y2 >= this.Sub_Strip.y2Max) = 2;
       end
	   function do1DPlotNormal(this,V,sym,deltaY)
            if(nargin < 4)
                deltaY = 0;
            end
            if(nargin < 3)
                sym = 'o';
            end 
            %V: Vector of length N2
            y2Max       = 5;
            mark        = (this.Pts.y1_kv == inf);            
            y2IP        = (0:0.01:y2Max)';
            [h_1,IP]      = ComputeInterpolationMatrix12(this,1,CompSpace2(this,y2IP));
            
            PtsCart     = GetCartPts(this);
            y2IPCart    = y2IP*sin(this.alpha);            
            plot(y2IPCart+deltaY,IP*V,'k','linewidth',1.5); hold on;
            if(~isempty(sym))
                plot(PtsCart.y2_kv(mark)+deltaY,V,sym,'MarkerEdgeColor','k','MarkerFaceColor','g'); 
            end            
            xlim([min(y2IPCart) max(y2IPCart)]);
            xlabel('$y_{2,Cart}$','Interpreter','Latex','fontsize',25);            
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
       end    
       function do1DPlotParallel(this,V)  
            %V: Vector of length N2
            y1W         = 5;
            mark        = (this.Pts.y2_kv == inf);
            y1IP        = (-y1W:0.1:y1W)';
            IP          = ComputeInterpolationMatrix12(this,CompSpace1(this,y1IP),1);
            
            %IP          = ComputeInterpolationMatrix(this,(-x1W:0.01:x1W)',1,true);
            IP.InterPol = IP.InterPol(:,mark);            
            
            plot(this.Pts.y1_kv(mark),V,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); 
            hold on;
            plot(y1IP,IP.InterPol*V,'linewidth',1.5);
            xlim([min(y1IP) max(y1IP)]);
            xlabel('$y_{1}$','Interpreter','Latex','fontsize',25);        
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
       end 
       function PlotGridLines(this)
            this.Sub_HalfSpace.PlotGridLines(); hold on;
            this.Sub_Strip.PlotGridLines();
       end           
       function PlotGrid(this)
            this.Sub_HalfSpace.PlotGrid(); hold on;
            this.Sub_Strip.PlotGrid();
            plot([-1 1]*1000,min(this.Sub_HalfSpace.GetCartPts.y2_kv)*[1,1],'k','linewidth',2);
       end        
       
       function PlotIsoline(this,x,y1y2)
           this.Sub_HalfSpace.PlotIsoline(x,y1y2);
           this.Sub_Strip.PlotIsoline(x,y1y2);
       end
       function Diff   = ComputeDifferentiationMatrix(this)
           Diff_ST  = this.Sub_Strip.ComputeDifferentiationMatrix();
           Diff_HS  = this.Sub_HalfSpace.ComputeDifferentiationMatrix();
           
           Diff.Dy1 = blkdiag(Diff_ST.Dy1,Diff_HS.Dy1);
           Diff.Dy2 = blkdiag(Diff_ST.Dy2,Diff_HS.Dy2);                      
       end
       
%        function PlotGrid(this)
%             scatter(this.Sub_Strip.Pts.y1_kv,...
%                     this.Sub_Strip.Pts.y2_kv,'ob');    hold on;
%             scatter(this.Sub_HalfSpace.Pts.y1_kv,...
%                     this.Sub_HalfSpace.Pts.y2_kv,'or');
%             plot([-7 7],[1. 1.],'b');
%             plot([-7 7],[0. 0.],'k','linewidth',2);
%             plot([-7 7],[0.5 0.5],'--b','linewidth',1.2);
%             xlim([-7 7]);
%             ylim([-0.5 3]);
%             pbaspect([14 3.5 1]);
%             set(gca,'fontsize',15);
%        end
   end
end