 classdef InfCapillary_Composed < ConvolutionFiniteSupport
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Always Check: The y1-Map (PhysSpace1) has to be the same for SubStrip
% and SubHalfSpace
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   properties 
       Bottom_Strip     %Area_ID: 0
       Main_Strip      %Area_ID: 1
       Top_Strip        %Area_ID: 2
       %mark_id_1,mark_id_2,mark_id,mark_12
       Pts,N1,N2
       M,MB,MM,MT     %no of points (in both dimensions)                       
   end
   
   methods
       function this = InfCapillary_Composed(Geometry)
           %Geometry includes:
           % y2Min,h,L1,L2,N,N2bound
           
           shapeStripB     = struct('y2Min',Geometry.y2Min - Geometry.R,...
                                    'y2Max',Geometry.y2Min + Geometry.R,...
                                    'L1',Geometry.L1,...
                                    'N',[Geometry.N(1),Geometry.N2bound],...
                                    'alpha',pi/2);                                
           this.Bottom_Strip  = InfCapillarySkewed(shapeStripB); 
           
           %'L2',Geometry.L2,
           shapeMainStrip  = struct('y2Min',Geometry.y2Min + Geometry.R,...
                                    'y2Max',Geometry.y2Max - Geometry.R,...                                    
                                    'L1',Geometry.L1,...
                                    'N',Geometry.N,...
                                    'alpha',pi/2);
           this.Main_Strip = InfCapillarySkewed(shapeMainStrip);
           
           shapeStripT     = struct('y2Min',Geometry.y2Max - Geometry.R,...
                                    'y2Max',Geometry.y2Max + Geometry.R,...
                                    'L1',Geometry.L1,...
                                    'N',[Geometry.N(1),Geometry.N2bound],...
                                    'alpha',pi/2);           
                                
           this.Top_Strip  = InfCapillarySkewed(shapeStripT);                       
           
           this.Pts.y1_kv   = [this.Bottom_Strip.Pts.y1_kv;...
                               this.Main_Strip.Pts.y1_kv;...
                               this.Top_Strip.Pts.y1_kv];
                           
           this.Pts.y2_kv   = [this.Bottom_Strip.Pts.y2_kv;...
                               this.Main_Strip.Pts.y2_kv;...
                               this.Top_Strip.Pts.y2_kv];
                           
           this.Pts.y1      = this.Pts.y1_kv(this.Pts.y2_kv == min(this.Pts.y2_kv));
           this.Pts.y2      = this.Pts.y2_kv(this.Pts.y1_kv == inf);
           
           this.Pts.x1      = this.Main_Strip.Pts.x1;
           
           this.N1          = Geometry.N(1);
           this.N2          = Geometry.N(2) + 2*Geometry.N2bound;
                           
           this.MB     = length(this.Bottom_Strip.Pts.y1_kv);
           this.MM     = length(this.Main_Strip.Pts.y1_kv);
           this.MT     = length(this.Top_Strip.Pts.y1_kv);
           this.M      = this.MB + this.MM + this.MT;
           
           
           %There are four numberings:
           %(1) the numberings of the domains. The points with identifyer
           %     X are given by mark_id(:,X)
           this.mark_id      = true(this.M,3);
           this.mark_id(:,1) = [true(this.MB,1);false(this.MM,1);false(this.MT,1)];
           this.mark_id(:,2) = [false(this.MB,1);true(this.MM,1);false(this.MT,1)];
           this.mark_id(:,3) = [false(this.MB,1);false(this.MM,1);true(this.MT,1)];           
           
           %(2) the numbering of the domains in the first dimension. Points
           %in the first dimension belonging to the part with identifyer Y
           %are given by mark_id_1(:,Y)
           this.mark_id_1   = true(this.N1,1);
           
           %(3) the numbering of the domains in the second dimension. Points
           %in the second dimension belonging to the part with identifyer Z
           %are given by mark_id_2(:,Z)           
           this.mark_id_2       = true(this.N2,2);                     
           this.mark_id_2(:,1)  = [true(this.Bottom_Strip.N2,1);false(this.Main_Strip.N2,1);false(this.Top_Strip.N2,1)];
           this.mark_id_2(:,2)  = [false(this.Bottom_Strip.N2,1);true(this.Main_Strip.N2,1);false(this.Top_Strip.N2,1)];
           this.mark_id_2(:,3)  = [false(this.Bottom_Strip.N2,1);false(this.Main_Strip.N2,1);true(this.Top_Strip.N2,1)];           
           
           %mark12 makes a link between identifyers of dimensions 1,2 and 
           % the global identifyer. Points which are in the first dimension
           % located in a domain with identifyer X and in the second
           % direction in a domain with identifyer Y, belong to the global
           % domain with the identifyer Z = mark_12(X,Y).
           this.mark_12      = zeros(1,3);
           this.mark_12(1,1) = 1;
           this.mark_12(1,2) = 2;
           this.mark_12(1,3) = 3;
       end

        function ptsCart = GetCartPts(this,pts_y1,pts_y2)            
             if(nargin == 1)
                 pts_y1 = this.Pts.y1_kv;
                 pts_y2 = this.Pts.y2_kv;
             end                                              
%             [ptsCart.y1_kv,ptsCart.y2_kv] = SkewedGrid(pts_y1,pts_y2,this.alpha);
             ptsCart.y1_kv = pts_y1;
             ptsCart.y2_kv = pts_y2;             
        end    
        function pts = GetInvCartPts(this,ptsCart_y1,ptsCart_y2)                                              
             %[pts.y1_kv,pts.y2_kv] = InvSkewedGrid(ptsCart_y1,ptsCart_y2,this.alpha);
             pts.y1_kv = ptsCart_y1;
             pts.y2_kv = ptsCart_y2;             
        end              
         
        function plot(this,V,opts)                      

           xl = [-5 5];           
           yl = [this.Bottom_Strip.y2Min this.Top_Strip.y2Max]
           N  = 50;
           
           PlotArea = struct('y1Min',xl(1),'y1Max',xl(2),...
                             'y2Min',this.Bottom_Strip.y2Min,...
                             'y2Max',this.Bottom_Strip.y2Max,...
                             'N1',N,'N2',N);
           this.Bottom_Strip.InterpolationPlot(PlotArea,true);           
           
           PlotArea.y2Min = this.Main_Strip.y2Min;
           PlotArea.y2Max = this.Main_Strip.y2Max;           
           this.Main_Strip.InterpolationPlot(PlotArea,true);
                      
           PlotArea.y2Min = this.Top_Strip.y2Min;
           PlotArea.y2Max = this.Top_Strip.y2Max;           
           this.Top_Strip.InterpolationPlot(PlotArea,true);
           
           if((nargin > 2) && (ischar(opts)))
                this.Bottom_Strip.plot(V(this.mark_id(:,1)),opts); hold on;           
                this.Main_Strip.plot(V(this.mark_id(:,2)),opts);hold on;      
                this.Top_Strip.plot(V(this.mark_id(:,3)),opts);
           else
                this.Bottom_Strip.plot(V(this.mark_id(:,1))); hold on;
                this.Main_Strip.plot(V(this.mark_id(:,2)));hold on;      
                this.Top_Strip.plot(V(this.mark_id(:,3)));
           end
                      
           xlim(xl);
           ylim(yl);
           pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);
          
       end       
        function IP = SubShapePtsCart(this,a_shapePts)
             pts = GetInvCartPts(this,a_shapePts.y1_kv,a_shapePts.y2_kv);
             IP  = SubShapePts(this,pts);  
             %IP  = InterpolationMatrix_Pointwise(this,pts.y1_kv,pts.y2_kv);                           
        end 
        function IP = SubShapePts(this,a_shapePts)
            maskB = (a_shapePts.y2_kv < this.Bottom_Strip.y2Max);
            maskM = ((a_shapePts.y2_kv >= this.Main_Strip.y2Min) & ...
                     (a_shapePts.y2_kv <= this.Main_Strip.y2Max));
            maskT = (a_shapePts.y2_kv > this.Top_Strip.y2Min);
           
           if( sum(~xor(xor(maskB,maskM),maskT)) ~= 0)
               exc = MException('InfCapillary_Composed:SubShapePts',...
                         'Points outside domain or points double counted');
               throw(exc);
           end           
           
           IP = zeros(length(a_shapePts.y1_kv),this.M);
           
           
           IP(maskB,this.mark_id(:,1))  = this.Bottom_Strip.InterpolationMatrix_Pointwise(...
                                            a_shapePts.y1_kv(maskB),...
                                            a_shapePts.y2_kv(maskB));
                                        
           IP(maskM,this.mark_id(:,2))  = this.Main_Strip.InterpolationMatrix_Pointwise(...
                                            a_shapePts.y1_kv(maskM),...
                                            a_shapePts.y2_kv(maskM));
                                        
           IP(maskT,this.mark_id(:,3))  = this.Top_Strip.InterpolationMatrix_Pointwise(...
                                            a_shapePts.y1_kv(maskT),...
                                            a_shapePts.y2_kv(maskT));
                                                     
       end
%         
        %X = InterpolateAndIntegratePtsOriginGrid(this,ptsOrGrid,dataDisk,weights);
%        [Int,Int1,Int2] = ComputeIntegrationVector(this);
%        
         function [Interp1,Interp2,mark_id] = ComputeInterpolationMatrix12(this,interp1,interp2)           
             
             Interp2 = zeros(length(interp2),this.N2);           
             mark_1  = (interp2(:,2)==1);           
             mark_2  = (interp2(:,2)==2);
             mark_3  = (interp2(:,2)==3);
              
             Interp1                              = barychebevalMatrix(this.Pts.x1,interp1); 
             Interp2(mark_1,this.mark_id_2(:,1))  = barychebevalMatrix(this.Bottom_Strip.Pts.x2,interp2(mark_1,1)); 
             Interp2(mark_2,this.mark_id_2(:,2))  = barychebevalMatrix(this.Main_Strip.Pts.x2,interp2(mark_2,1)); 
             Interp2(mark_3,this.mark_id_2(:,3))  = barychebevalMatrix(this.Top_Strip.Pts.x2,interp2(mark_3,1));                           
             
             mark_id = interp2(:,2);
            
         end 
        function [x1] = CompSpace1(this,y1)
             x1 = this.Main_Strip.CompSpace1(y1);
        end 
 	    function [x2] = CompSpace2(this,y2)
            
            mark_id = GetArea_ID(this,1,y2);
           
            x2(mark_id==1,1) = this.Bottom_Strip.CompSpace2(y2(mark_id==1));
            x2(mark_id==2,1) = this.Main_Strip.CompSpace2(y2(mark_id==2));
            x2(mark_id==3,1) = this.Top_Strip.CompSpace2(y2(mark_id==3));
            x2(:,2)          = mark_id;          
        end    
%        function [I12] = Combine12(this,I1,I2)
%            %I1 - linear Operator on first dimension
%            %I2 - linear Operator on second dimension           
%            I12 = [kronecker(I1,I2(:,this.mark_id_2{1})),kronecker(I1,I2(:,this.mark_id_2{2}))];
%        end
       function id = GetArea_ID(this,y1,y2)
            id = ones(size(y2));
            
            id(y2 < this.Bottom_Strip.y2Max)    = 1;
            id((y2 >= this.Main_Strip.y2Min) & ...
               (y2 <= this.Main_Strip.y2Max))   = 2;
            id(y2 > this.Top_Strip.y2Min)       = 3;
            
        end
	   function do1DPlotNormal(this,V)
            alpha = pi/2 ; %this.alpha
            
            %V: Vector of length N2            
            mark        = (this.Pts.y1_kv == inf);            
            y2IP        = (this.Bottom_Strip.y2Min:0.1:this.Top_Strip.y2Max)';
            [h_1,IP]      = ComputeInterpolationMatrix12(this,1,CompSpace2(this,y2IP));
            
            PtsCart     = GetCartPts(this);
            y2IPCart    = y2IP*sin(alpha);            
            plot(PtsCart.y2_kv(mark),V,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); 
            hold on;
            plot(y2IPCart,IP*V,'linewidth',1.5);
            xlim([min(y2IPCart) max(y2IPCart)]);
            xlabel('$y_{2,Cart}$','Interpreter','Latex','fontsize',25);            
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
       end    
%        function do1DPlotParallel(this,V)  
%             %V: Vector of length N2
%             y1W         = 5;
%             mark        = (this.Pts.y2_kv == inf);
%             y1IP        = (-y1W:0.1:y1W)';
%             IP          = ComputeInterpolationMatrix12(this,CompSpace1(this,y1IP),1);
%             
%             %IP          = ComputeInterpolationMatrix(this,(-x1W:0.01:x1W)',1,true);
%             IP.InterPol = IP.InterPol(:,mark);            
%             
%             plot(this.Pts.y1_kv(mark),V,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); 
%             hold on;
%             plot(y1IP,IP.InterPol*V,'linewidth',1.5);
%             xlim([min(y1IP) max(y1IP)]);
%             xlabel('$y_{1}$','Interpreter','Latex','fontsize',25);        
%             set(gca,'fontsize',20);                        
%             set(gca,'linewidth',1.5);       
%        end 
%        function PlotGridLines(this)
%             this.Sub_HalfSpace.PlotGridLines(); hold on;
%             this.Sub_Strip.PlotGridLines();
%        end
   end
end