function ddfttest(temp)

x1=temp.Interp.pts1;
x2=temp.Interp.pts2;

N1=temp.Interp.Nplot1;
N2=temp.Interp.Nplot2;

rho=temp.rho_IP(:,end);

y1 = x1.*cos(x2); 
y2 = x1.*sin(x2);

y1=reshape(y1,N1,N2);
y2=reshape(y2,N1,N2);

rho=reshape(rho,N1,N2);

surf(y1,y2,rho)

end