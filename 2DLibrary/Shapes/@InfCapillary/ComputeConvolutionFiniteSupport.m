function AD = ComputeConvolutionFiniteSupport(this,area,weights,pts,params)    
%%
%
% This function computes the convolution between a function defined on 
% an infinite capillary with a set of functions with of finite support, at a specific
% set of points 

%% Input
%
% * area - object of a class with two methods 
%                   (1) [int,A] = ComputeIntegrationVector() and
%                   (2) pts = GetCartPts(), with pts = struct(y1_kv,y2_kv)
% * weights - a cell with a list of functions which get a structure with 
%         two arrays [y1_kv,y2_kv] as input and returns one
%         array of same size. [y1_kv,y2_kv] represents a point in polar
%         coordinates, representing the radial and angular
%         component, respectively.
% * pts - structure with 'y1_kv','y2_kv','y1','y2'. 
%           (y1_kv,y2_kv) represents a point in the skewed coordinate
%           system, and is a grid defined through [y1 (X) y2]
%
%% Output
%
% Average densities to get average densities
%
% $$X_{i,:,k}\cdot \rho = \int_{A}  \rho({\bf r}_i+{\bf r}')w_k({\bf r}') d{\bf r}'$$
%
% where
% $A = area.GetCartPts() \cap this.GetCartPts()$, and ${\bf r}_i \in
% this.AD.GetCartPts()$. 
% $w_k$ is the k-th element in the weight-function
% list _weights_.
%
%           
% Note: This method is only tested if pts are in a shape of an infinite capillary.
    
    %% Initialization
    %
    % For the input points given through pts, we distinguish
    % two sets: 
    %
    % # points r in _pts_, for which r+rd with rd in area.GetCartPts() are always
    %   within the area defined through _this_ class, so the area
    %   does not have to be split. These are collected in ptsHS.
    % # points r in pts, for which r+rd with rd in area.GetCartPts() are not
    % always a subset of _this_ class. These are collected in
    % ptsStrip.
    %
    
    if(isprop(this,'alpha'))
        alpha = this.alpha; %pi/2;%
    else
        alpha = pi/2;
    end

    fprintf('Computing interpolation for matrices for averaged densities..\n');

    AD          = zeros(length(pts.y1_kv),this.M,numel(weights)+1);%always include unity weight
    areaPtsCart = area.GetCartPts();

    y2SepMin    = min(this.Pts.y2) + abs(min(areaPtsCart.y2_kv))/sin(alpha);
    y2SepMax    = max(this.Pts.y2) - abs(max(areaPtsCart.y2_kv))/sin(alpha);
                       
                       
    % Note: division through sin(this.alpha) takes into account that
    % pts.y1_kv and pts.y2_kv are the coordinates in the skewed
    % grid.
    mark_Y2{1}  = (pts.y2    < y2SepMin);
    mark_Ykv{1} = (pts.y2_kv < y2SepMin);            
    
    mark_Y2{2}  = ((pts.y2    >= y2SepMin) & (pts.y2  <= y2SepMax));
    mark_Ykv{2} = ((pts.y2_kv >= y2SepMin) & (pts.y2_kv <= y2SepMax));            

    mark_Y2{3}  = ((pts.y2    > y2SepMax) & (pts.y2    > y2SepMin));
    mark_Ykv{3} = ((pts.y2_kv > y2SepMax) & (pts.y2_kv > y2SepMin));

    for i = 1:3
        ptsStrip{i}.y1_kv = pts.y1_kv(mark_Ykv{i});
        ptsStrip{i}.y2_kv = pts.y2_kv(mark_Ykv{i});
        ptsStrip{i}.y2    = pts.y2(mark_Y2{i});
        ptsStrip{i}.y1    = pts.y1;
    end
   
    %% Computation of Intersections
    %
    % <<DraftX.png>>
    %
    % The intersection between {area+r} with r in ptsStrip and the
    % HalfSpace (the _this_ class) are computed. An
    % intersection for {area+r_i} and {area+r_j} is equivalent up to
    % translation, if the y-coordinates of r_i and r_j are the same.
    % We make use of this property, for performance purposes, as the
    % points in ptsStrip are aligned with a grid.
    
    ptsStripCart = GetCartPts(this,0,ptsStrip{1}.y2);
    ptsy2        = ptsStripCart.y2_kv;
    for iPts = 1:length(ptsy2)
        area.Origin(2) = ptsy2(iPts);
        dataAD1(iPts)  = Intersect(this,area);
    end    
    area.Origin(2) = 0;

    ptsStripCart = GetCartPts(this,0,ptsStrip{3}.y2);
    ptsy2        = ptsStripCart.y2_kv;        
    for iPts = 1:length(ptsy2)
        area.Origin(2)  = ptsy2(iPts);
        dataAD3(iPts)   = Intersect(this,area);
    end    
    area.Origin(2) = 0;
    
    %% Computation of convolution matrices
    %
    % For the point in the set ptsHS, the integration matrix (over
    % area) is the same for each point. To see this, consider that in
    %
    % $$\int_A \rho({\bf r}_i+{\bf r}') w_k({\bf r}') d{\bf r}'$$
    %
    % for r in ptsHSthe area A is independent of r_i.     
  
    if(nargin==5)
        AD(mark_Ykv{1},:,:)  = Conv_LinearGridX(this,ptsStrip{1},dataAD1,weights,params);
        AD(mark_Ykv{2},:,:)  = Conv_LinearGridXY(this,ptsStrip{2},area,weights,params);    
        AD(mark_Ykv{3},:,:)  = Conv_LinearGridX(this,ptsStrip{3},dataAD3,weights,params);
    else      
        if(~isempty(ptsStrip{1}.y2))
            AD(mark_Ykv{1},:,:)  = Conv_LinearGridX(this,ptsStrip{1},dataAD1,weights);
        end
        if(~isempty(ptsStrip{2}.y2))
            AD(mark_Ykv{2},:,:)  = Conv_LinearGridXY(this,ptsStrip{2},area,weights);
        end
        if(~isempty(ptsStrip{3}.y2))
            AD(mark_Ykv{3},:,:)  = Conv_LinearGridX(this,ptsStrip{3},dataAD3,weights);        
        end
    end
    
end      