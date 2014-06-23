function plotParticles3D(r,p,opts,handles)
%plotParticles3D(r,opts,handles)
% plots spheres at given points, can take multiple inputs (given by
% multiple columns of r) and shift the data in a specified way to avoid
% overlaps.  
% Note: r must be 3D data
%
% INPUTS:
%  r        --  matrix of size (3*nParticles) x nPlots of positions
%
%  opts     --  structure of size [1 1] containing
%                 sigma         (particle diameter)
%                 colours       (colours for each plot, of the form 
%                                    {{'r',...,'b'}})
%                 shift         (3xnPlots matrix to shift the plots by)
%                 renormalize   (3x1 vector of true/false whether to 
%                                    renormalize coordinate)
%                 relabel       (3x1 vector of true/false whether to 
%                                    relabel axes)
%                 lims          (cell of three 1x2 vectors of limits,
%                                 {{[xmin,xmax],[ymin,ymax],[zmin,zmax]}} )
%                 ticks         (cell of tick mark locations given 
%                                     by row vectors, e.g.{{[0 2],[],[]}} )
%                 labs          (cell of list of tick mark labels given 
%                                  as column vectors,
%                                  e.g. {{{'HI off'; 'HI on'},{''},{''}}} )
%                 followType    ('max'/'min', whether to follow the max or
%                                 min of a coordinate when using renormalize)
%                 dpi           (integer resolution)
%
%  handles  --  structure of size [1 1] containing
%                 hf            (figure handle)
%                 ha            (axis handle)
%               or can be [] and a new figure will be used

% get figure and axis handles
if(isempty(handles))
    hf=figure;
    ha=axes;
    hold(ha,'on');
else
    hf=handles.hf;
    ha=handles.ha;
end

% default renderer has problems with background colour and axes
set(hf,'Renderer','zbuffer');
% vector graphics but no lighting:
%set(hf,'Renderer','painters');

% stop print from setting the fluid colour to white
set(hf,'InvertHardcopy','off')


% particle diameter
sigma=opts.sigma;
sigma=diag(sigma);

% colours for each set of particles
colours=opts.colours;

% shift for each set of particles, used to prevent them overlapping
shift=opts.shift;

% get dimensions
nParticles=size(r,1)/3;
nPlots=size(r,2);

% shift each plot so they don't overlap
if(~isempty(shift))
    shift=repmat(shift,nParticles,1);
    %shift=repmat(shift,nParticles,nPlots);
    r=r+shift;
end

% make a unit sphere
numSphereFaces =100;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);

% determine whic point to follow, at the moment either the fastest or
% slowest particle
followType=opts.followType;

if(strcmp(followType,'max'))
    % find maximum over all calculations for each particle
    bounds=max(r,[],2);
    % then maximum over each particle
    bounds=reshape(bounds,3,[]);
    bound=max(bounds,[],2);    
    % make the shift in each coordinate over all plots
    boundAll=repmat(bound,nParticles,nPlots);
elseif(strcmp(followType,'min'))
    bounds=min(r,[],2);
    bounds=reshape(bounds,3,[]);
    bound=min(bounds,[],2);
    boundAll=repmat(bound,nParticles,nPlots);
elseif(strcmp(followType,'com'))
    rTemp=reshape(r,[3,nParticles,nPlots]);
    rTemp=rTemp(:,:,opts.comPlot);
    com=sum(rTemp,2)/nParticles;
    boundAll=repmat(com,nParticles,nPlots);
else
    fprintf(1,'Pick a valid followType\n');
    pause;
end
    
% find which axes to renormalize
renormalize=opts.renormalize;
% copy for each particle
renormalize=repmat(renormalize,nParticles,1);

% rescale by the shift
r(renormalize,:)=r(renormalize,:)-boundAll(renormalize,:);
% note this allows us to set the axis limits to be constant which produces
% much nicer plots, stopping things like horizontal rescaling due to
% different label sizes.


axes(ha); % surface doesn't take axis handles
hold(ha,'on');

% equal axes so the spheres are actually spheres
axis equal

% set axis limits
xlim(ha,opts.lims{1});
ylim(ha,opts.lims{2});
zlim(ha,opts.lims{3});

% and viewpoint
view(ha,opts.viewPoint);

% set boxed axes which looks nicer
box(ha,'on')

for iPlot=1:nPlots
    
    % get x, y and z coordinates
    xyz=reshape(r(:,iPlot),3,[]);
    x=xyz(1,:);
    y=xyz(2,:);
    z=xyz(3,:);
    
    for iParticle=1:nParticles
        
        s=sigma(iParticle);  % assume in 3D
        
        % calculate a sphere for each particle of each plot
        sphereX = x(iParticle) + unitSphereX*s/2;
        sphereY = y(iParticle) + unitSphereY*s/2;
        sphereZ = z(iParticle) + unitSphereZ*s/2;
        
        % and plot it in the corresponding colour
        hSurf=surface(sphereX, sphereY, sphereZ);
        set(hSurf,'FaceColor',colours{iPlot}(:,iParticle),'EdgeColor','none');
        
        % set up lighting - this is somewhat trial and error
        set(hSurf, ...
        'FaceLighting','phong', ...
        'AmbientStrength',0.3, ...
        'DiffuseStrength',0.8, ...
        'SpecularStrength',0.9, ...
        'SpecularExponent',5, ... %20
        'BackFaceLighting','lit');
    
    end
    
    % get x,y,z components of velocity
    uvw=reshape(p(:,iPlot),3,[]);
    u=uvw(1,:);
    v=uvw(2,:);
    w=uvw(3,:);
    % and magnitude and remove errors when it's zero
    P=sqrt(sum(uvw.^2,1))+eps;
    
    % calculate points at which to place the arrows
    xa=x+u/P.*sigma'/2;
    ya=y+v/P.*sigma'/2;
    za=z+w/P.*sigma'/2;
    
    % plot arrows
    for iParticle=1:nParticles
        quiver3(ha,xa(iParticle),ya(iParticle),za(iParticle),u(iParticle),v(iParticle),w(iParticle) ...
            ,0,'Color',colours{iPlot}(:,iParticle),'MaxHeadSize',0.5)
    end
    
end

% set lighting position based on camera position
%camlight('headlight');
%camlight('right');
camlight('left');

% whether to relabel axes
relabel=opts.relabel;

% and tickmarks
if(relabel(1))
    set(ha,'XTick',opts.ticks{1});
    set(ha,'XTickLabel',opts.labs{1});
end
if(relabel(2))
    set(ha,'YTick',opts.ticks{2});
    set(ha,'YTickLabel',opts.labs{2});
end
if(relabel(3))
    set(ha,'ZTick',opts.ticks{3});
    set(ha,'ZTickLabel',opts.labs{3});
end

% and title
if(isfield(opts,'title'))
    title(ha,opts.title);
end

% set background colour to a light blue
set(ha,'Color',opts.bath);


% make sure hold is off before we exit
hold(ha,'off');
end
