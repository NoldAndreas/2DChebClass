function AD = ComputeConvolutionFiniteSupport(this,area,weights,pts,params)    
%%
%
% This function computes the convolution between a function defined on 
% a HalfSpace with a set of functions with of finite support, at a specific
% set of points 

%% Input
%
% * area - structure with two methods (1) [int,A] = ComputeIntegrationVector() and
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
% Note: This method is only tested if pts are in a shape of a HalfSpace.
    
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

    fprintf('Computing interpolation for matrices for averaged densities..\n');
    tic

    AD          = zeros(length(pts.y1_kv),this.M,numel(weights)+1);%always include unity weight
    areaPtsCart = area.GetCartPts();

    y2Wall   = min(this.Pts.y2);
                %this.GetInvCartPts(0,min(this.GetCartPts.y2_kv)).y2_kv; 
                %this.y2wall + this.R/sin(this.alpha)    
    y2Sep    = y2Wall + abs(min(areaPtsCart.y2_kv))/sin(this.alpha);
    % Note: division through sin(this.alpha) takes into account that
    % pts.y1_kv and pts.y2_kv are the coordinates in the skewed
    % grid.
    markY2  = (pts.y2    < y2Sep);
    markYkv = (pts.y2_kv < y2Sep);            

    ptsStrip.y1_kv = pts.y1_kv(markYkv);
    ptsStrip.y2_kv = pts.y2_kv(markYkv);
    ptsStrip.y2    = pts.y2(markY2);
    ptsStrip.y1    = pts.y1;

    ptsHS.y1_kv = pts.y1_kv(~markYkv);
    ptsHS.y2_kv = pts.y2_kv(~markYkv);
    ptsHS.y2    = pts.y2(~markY2);
    ptsHS.y1    = pts.y1;

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
    
    ptsStripCart = GetCartPts(this,0,ptsStrip.y2);
    ptsy2        = ptsStripCart.y2_kv;

    for iPts = 1:length(ptsy2)
        area.Origin(2) = ptsy2(iPts);
        dataAD(iPts)   = Intersect(this,area);
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
        if(sum(markYkv)>0)
            AD(markYkv,:,:) = Conv_LinearGridX(this,ptsStrip,dataAD,weights,params);
        end
        AD(~markYkv,:,:)  = Conv_LinearGridXY(this,ptsHS,area,weights,params);    
    else
        if(sum(markYkv)>0)
            AD(markYkv,:,:)   = Conv_LinearGridX(this,ptsStrip,dataAD,weights);        
        end
        AD(~markYkv,:,:)  = Conv_LinearGridXY(this,ptsHS,area,weights);    
    end

    t = toc;
    disp([num2str(t),'s']);  
end      