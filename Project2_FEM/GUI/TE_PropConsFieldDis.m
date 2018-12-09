function [wavenumber_z]=TE_PropConsFieldDis(wd,ht,Nx,Ny,epsi_r,mode,haxes2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Input necessary data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_RW=wd;
b_RW=ht;
stepX=Nx-1;
stepY=Ny-1;
m=mode;
freq = 0.1e9:0.01e9:10e9;

w = 2.*pi.*freq;

mu0=4*pi*1.0e-7;  
epsi0=8.854e-12; 
% epsi_r=1; % relative permittivity
epsi_die=epsi0*epsi_r;
c_die=1/sqrt(mu0*epsi_die);

dx = a_RW/stepX;
dy = b_RW/stepY;

[x,y] = meshgrid(0:dx:a_RW,0:dy:b_RW);     %Creating the nodes cooredinates
DT = delaunayTriangulation(x(:),y(:)); 
ConnectivityArray = DT.ConnectivityList;   %Creating the connectivity array
Nn = max(DT.ConnectivityList(:));          %Total number of nodes
Points = DT.Points;                        %Nodes coordinates
temp = size(DT.ConnectivityList);
Ne = temp(1);                              %Total number of elements
% figure(1);
% triplot(DT);
% axis equal;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Initialize matrix A & B (Nn x Nn) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = sparse(Nn,Nn);
B = sparse(Nn,Nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assemble all triangular elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ee=1:1:Ne      %looping through all the elements
    for nn=1:1:3   %looping through nodes of the element
        node_global_ID(nn) = ConnectivityArray(ee,nn);    %getting node global ID
        node_coordinates_x(nn) = Points(node_global_ID(nn),1);  %getting node x-coordinates
        node_coordinates_y(nn) = Points(node_global_ID(nn),2);  %getting node y-coordinates 
    end
    
    a(1) = node_coordinates_x(2)*node_coordinates_y(3)-node_coordinates_y(2)*node_coordinates_x(3);
    a(2) = node_coordinates_x(3)*node_coordinates_y(1)-node_coordinates_y(3)*node_coordinates_x(1);
    a(3) = node_coordinates_x(1)*node_coordinates_y(2)-node_coordinates_y(1)*node_coordinates_x(2);
    
    b(1) = node_coordinates_y(2)-node_coordinates_y(3);
    b(2) = node_coordinates_y(3)-node_coordinates_y(1);
    b(3) = node_coordinates_y(1)-node_coordinates_y(2);
    
    c(1) = node_coordinates_x(3)-node_coordinates_x(2);
    c(2) = node_coordinates_x(1)-node_coordinates_x(3);
    c(3) = node_coordinates_x(2)-node_coordinates_x(1);
    
    delta = 0.5*(b(1)*c(2)-b(2)*c(1));
    
    for rr_element=1:1:3   %looping through rows of the matrix A_element and B_element
        rr = ConnectivityArray(ee,rr_element);  %row location in the ugmented matrix A and B
        for cc_element = 1:1:3    %looping through columns of the matrix A_element and B-element
            cc = ConnectivityArray(ee,cc_element);  %column location in the ugmented matrix A and B
            A_element(rr_element,cc_element) = (0.25/delta)*(b(rr_element)*b(cc_element)+c(rr_element)*c(cc_element));
            if rr_element == cc_element
                B_element(rr_element,cc_element) = (delta/12)*(1+1);
            else
                B_element(rr_element,cc_element) = (delta/12);
            end
        A(rr,cc) = A(rr,cc) + A_element(rr_element,cc_element);
        B(rr,cc) = B(rr,cc) + B_element(rr_element,cc_element);
        end
    end 
    
end

%% Calculate propagation const
% wavenumber_cutoff = sqrt(eig(B\A));
% % sort kc in ascending order
% wavenumber_cutoff=sort(wavenumber_cutoff); 

%  [V,D]=eigs(A,B,6,'smallestreal');
 [V,D]=eigs(A,B,6,'sr');
 wavenumber_cutoff=sqrt(D);

% plot the first N_mode TE mode
N_mode=6;
str=cell(N_mode,1);
% % figure(2);
% m = 2; %size(wavenumber_cutoff)
wavenumber_z = sqrt((w/c_die).^2-wavenumber_cutoff(m,m)^2);
% plot(freq/1e9,wavenumber_z,'linewidth',3);
% hold on;
% grid on;
% str{m}=[num2str(m-1),'-th TE mode'];
% set(gca,'fontsize',16);
% xlabel('frequency ( GHz )','FontSize',16);
% ylabel('Progagation constant ( m^{-1} )','FontSize',16);
% legend(str(m),'Location','northwest');
% hold off;


%% Plot the Field distribution
Hz=reshape(V(:,m),Ny,Nx);

% choose arbitrary kz, then the frequency is determined
kc=wavenumber_cutoff(m,m);
kz=1;
k_tot=sqrt(kc^2+kz^2);
omega=k_tot*c_die;

% Hx & Hy
Hx=zeros(Ny,Nx);
Hy=zeros(Ny,Nx);
Hx(:,2:Nx-1)=-kz*1.0i/kc/kc*(Hz(:,3:Nx)-Hz(:,1:Nx-2))/2/dx;
Hy(2:Ny-1,:)=-kz*1.0i/kc/kc*(Hz(3:Ny,:)-Hz(1:Ny-2,:))/2/dy;
Htot=sqrt(Hx.*conj(Hx)+Hy.*conj(Hy)+Hz.*conj(Hz));

% Ex & Ey
Ex=zeros(Ny,Nx);
Ey=zeros(Ny,Nx);
Ex(2:Ny-1,:)=-1.0i*omega*mu0/kc/kc*(Hz(3:Ny,:)-Hz(1:Ny-2,:))/2/dy;
Ey(:,2:Nx-1)=1.0i*omega*mu0/kc/kc*(Hz(:,3:Nx)-Hz(:,1:Nx-2))/2/dx;
Etot=sqrt(Ex.*conj(Ex)+Ey.*conj(Ey));

% normalize field value based on max(Etot)
Emax=max(max(Etot));
Ex=Ex/1.0i/Emax;
Ey=Ey/1.0i/Emax;
Hx=Hx/1.0i/Emax;
Hy=Hy/1.0i/Emax;

% set(haxes2);
axes(haxes2);

y_s=round(linspace(1,Nx,15)); % only plot several points
x_s=round(linspace(1,Ny,10));
quiver(x(x_s,y_s),y(x_s,y_s),Ex(x_s,y_s),Ey(x_s,y_s),'Color','r','AutoScaleFactor',1.2);
hold on;
quiver(x(x_s,y_s),y(x_s,y_s),Hx(x_s,y_s),Hy(x_s,y_s),'b:','AutoScaleFactor',1.2);
axis equal;
set(gca,'fontsize',16);
xlabel('Width ( m )','FontSize',16);
ylabel('Height ( m )','FontSize',16);
legend('E-field lines', 'H-field lines','Location','northwest');
title('Modal Field Distribution');
hold off;

% Etot=Etot/Emax;
% Htot=Htot/Emax;
% figure;
% contour(x,y,Htot,'ShowText','off','LineStyle',':','LineColor','b','LineWidth',1);
% hold on;
% contour(x,y,Etot,'ShowText','on','LineColor','r','LineWidth',1);
% axis equal;
% xlabel('x ( m )','FontSize',16);
% ylabel('y ( m )','FontSize',16);
% hold off;
end