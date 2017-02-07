classdef HoleInHalfSpace < handle
   
    properties 
        OriginHole = [0;0];
        RHole
        Strip = []
        LeftSlit,RightSlit,UpperHS               
        Int
        Pts     
        sphere
    end
    
    methods 
        function this = HoleInHalfSpace(Geometry)            
            %Elements of Geometry:
            %
            %-OriginHole
            %-RHole
            %-N
            %-L  
            %-y2Wall
            this.OriginHole = Geometry.OriginHole;            
            this.RHole      = Geometry.R;
                        
            %********************************
            shapeHS.y2Min        = Geometry.OriginHole(2) + Geometry.R;
            shapeHS.L1           = Geometry.L;
            shapeHS.L2           = Geometry.L;
            shapeHS.y10          = Geometry.OriginHole(1);
            shapeHS.N            = Geometry.N;
            this.UpperHS         = HalfSpace(shapeHS);
            
            this.Pts  = struct('y1_kv',this.UpperHS.Pts.y1_kv,...
                               'y2_kv',this.UpperHS.Pts.y2_kv);            
            
            shapeSlit.N          = Geometry.N;
            shapeSlit.L1         = Geometry.L;
            shapeSlit.R          = Geometry.R;
            shapeSlit.OriginDisk = Geometry.OriginHole;
            shapeSlit.y2Min      = max(Geometry.OriginHole(2)-Geometry.R,Geometry.y2Wall);
            shapeSlit.y2Max      = Geometry.OriginHole(2) + Geometry.R;
            shapeSlit.LeftRight  = 'Left';            
            this.LeftSlit        = SlitMinusDisk(shapeSlit);
            
            this.Pts.y1_kv    = [this.Pts.y1_kv;this.LeftSlit.Pts.y1_kv];
            this.Pts.y2_kv    = [this.Pts.y2_kv;this.LeftSlit.Pts.y2_kv];                
            
            shapeSlit.LeftRight  = 'Right';
            this.RightSlit       = SlitMinusDisk(shapeSlit);
            
            this.Pts.y1_kv    = [this.Pts.y1_kv;this.RightSlit.Pts.y1_kv];
            this.Pts.y2_kv    = [this.Pts.y2_kv;this.RightSlit.Pts.y2_kv];                
            
            %********************************            
            if(Geometry.OriginHole(2) - Geometry.R > Geometry.y2Wall)
                shapeStrip.y2Min  = Geometry.y2Wall;    
                shapeStrip.y2Max  = Geometry.OriginHole(2)-Geometry.R;
                shapeStrip.N      = Geometry.N;
                shapeStrip.L1     = Geometry.L;
                shapeStrip.y10    = Geometry.OriginHole(1);
                this.Strip        = InfCapillary(shapeStrip);            
                this.Pts.y1_kv    = [this.Pts.y1_kv;this.Strip.Pts.y1_kv];
                this.Pts.y2_kv    = [this.Pts.y2_kv;this.Strip.Pts.y2_kv];                
            end
            %********************************
        end        
        function ptsCart = GetCartPts(this)            
            ptsCart = this.Pts;
        end         
        function int = ComputeIntegrationVector(this)           
            %Strip = []
            %LeftSlit,RightSlit,UpperHS               
            IntUpperHS   = this.UpperHS.ComputeIntegrationVector();
            IntLeftSlit  = this.LeftSlit.ComputeIntegrationVector();
            IntRightSlit = this.RightSlit.ComputeIntegrationVector();
            int          = [IntUpperHS,IntLeftSlit,IntRightSlit];
            if(~isempty(this.Strip))    
                IntStrip = this.Strip.ComputeIntegrationVector();                        
                int      = [int,IntStrip];
            end
        end        
        function PlotGridLines(this)
           if(~isempty(this.Strip))
               this.Strip.PlotGridLines();
           end
           if(~isempty(this.LeftSlit))
               this.LeftSlit.PlotGridLines();
           end
           if(~isempty(this.RightSlit))
               this.RightSlit.PlotGridLines();
           end
           if(~isempty(this.UpperHS))
               this.UpperHS.PlotGridLines();
           end
        end
        function PlotGrid(this,InitFig)    

            ptsCart = GetCartPts(this);                        
            h = scatter(ptsCart.y1_kv,ptsCart.y2_kv,'ob');              
            set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');

            if((nargin == 2) && InitFig)
                y1lim = [min(ptsCart.y1_kv) max(ptsCart.y1_kv)];
                y2lim = [min(ptsCart.y2_kv) max(ptsCart.y2_kv)];

                xlim(y1lim + [-0.5 0.5]);
                ylim(y2lim + [-0.5 0.5]);

                pbaspect([(y1lim(2)-y1lim(1)) (y2lim(2)-y2lim(1)) 1]);
                set(gca,'fontsize',15);
            end
        end    

    end
    
    
end    

