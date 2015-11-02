% Simulating the 2-D Burgers' equation by the Finite Difference Method
%...(a time march)
% Numerical scheme used is a first order upwind in time and for the
%...first order partial space derivatives, while a central difference
%...for the second order spatial derivatives

%%
%Specifying Parameters
nx=30;                           %Number of steps in space(x)
ny=30;                           %Number of steps in space(y)       
nt=60;                           %Number of time steps 
dt=0.01;                         %Width of each time step
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
vis=0.01;                        %Diffusion coefficient/viscosity
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
v=zeros(nx,ny);                  %Preallocating v
vn=zeros(nx,ny);                 %Preallocating vn
rhsx=zeros(nx,ny);
rhsy=zeros(nx,ny);

%%
% Initial Conditions
for i=1:nx
    for j=1:ny
        if ((0.5<=y(j))&&(y(j)<=1)&&(0.5<=x(i))&&(x(i)<=1))
            u(i,j)=0;
            v(i,j)=1;
        else
            u(i,j)=1;
            v(i,j)=0;
        end
    end
end

%%
%Boundary conditions
u(1,:)=0;
u(nx,:)=0;
u(:,1)=0;
u(:,ny)=0;
v(1,:)=0;
v(nx,:)=0;
v(:,1)=0;
v(:,ny)=0;

%%
%Calculating the coefficient matrix for the implicit scheme
Ex=(vis*dt/dx^2)*sparse(2:nx,1:nx-1,1,nx,nx);
Dx=Ex+Ex'-(2*dt*vis/dx^2)*speye(nx);
Ey=(vis*dt/dy^2)*sparse(2:ny,1:ny-1,1,ny,ny);
Dy=Ey+Ey'-(2*dt*vis/dy^2)*speye(ny);
d=kron(Dy,speye(nx))+kron(speye(ny),Dx);
D=d-speye(nx*ny);

%%
%Calculating the velocity profile for each time step
i=2:nx-1;
j=2:ny-1;
m=2:nx;
n=2:ny;
for it=0:nt
    un=u;
    vn=v;
    %plotting the velocity field
    h=quiver(x,y,u',v','k');
    axis([0 2 0 2])
    axis square
    title({['2-D Burgers'' equation with {\nu} = ',num2str(vis)];'Transport property vector field {\bfu}=(u_x,u_y)';['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Spatial co-ordinate (y) \rightarrow')
    drawnow;
    refreshdata(h)
    %Uncomment as necessary
    %---------------
    %Implicit method:
    
    rhsx(m,n)=un(m,n)-(dt*un(m,n).*(un(m,n)-un(m-1,n))/dx)-(dt*vn(m,n).*(un(m,n)-un(m,n-1))/dy);
    rhsy(m,n)=vn(m,n)-(dt*un(m,n).*(vn(m,n)-vn(m-1,n))/dx)-(dt*vn(m,n).*(vn(m,n)-vn(m,n-1))/dy);
    bx=-reshape(rhsx,[],1);
    bx(1:nx)=0;bx(nx*(ny-1)+1:nx*ny)=0;         %Imposing boundary conditions
    bx(1:nx:nx*(ny-1)+1)=0;bx(nx:nx:nx*ny)=0;   ... on b itself
    solx=D\bx;
    u=reshape(solx,nx,ny);
    by=-reshape(rhsy,[],1);
    by(1:nx)=0;by(nx*(ny-1)+1:nx*ny)=0;         %Imposing boundary conditions
    by(1:nx:nx*(ny-1)+1)=0;by(nx:nx:nx*ny)=0;   ... on b itself
    soly=D\by;
    v=reshape(soly,nx,ny);
    %}
    %---------------
    %Explicit mathod
    %{
    u(i,j)=un(i,j)-(dt*(un(i,j)-un(i-1,j)).*un(i,j)/dx)-(dt*(un(i,j)-un(i,j-1)).*vn(i,j)/dy)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j-1)-2*un(i,j)+un(i,j+1))/(dy*dy));
    v(i,j)=vn(i,j)-(dt*(vn(i,j)-vn(i-1,j)).*un(i,j)/dx)-(dt*(vn(i,j)-vn(i,j-1)).*vn(i,j)/dy)+(vis*dt*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))/(dx*dx))+(vis*dt*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))/(dy*dy));
    %Boundary Conditions
    u(1,:)=0;
    u(nx,:)=0;
    u(:,1)=0;
    u(:,ny)=0;
    v(1,:)=0;
    v(nx,:)=0;
    v(:,1)=0;
    v(:,ny)=0;
    %}
end