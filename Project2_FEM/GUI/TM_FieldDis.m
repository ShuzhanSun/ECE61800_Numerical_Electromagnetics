function [kc_TM]=TM_FieldDis(wd,ht,Nx,Ny,epsi_r,haxes2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Input necessary data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_RW=wd;
b_RW=ht;
stepX=Nx-1;
stepY=Ny-1;
freq = 0.1e9:0.01e9:7e9;

w = 2.*pi.*freq;
freq = 0.1e9:0.01e9:10e9;
w = 2.*pi.*freq;
eps0 = 1e-9/(36*pi);
u0 = 4*pi*1e-7;
eps0=eps0*epsi_r;
c0=1/sqrt(eps0*u0);

w0 = 2*pi*8e9;
% Nx = 51;
% Ny = 51;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Initialize matrix A & B (Nn x Nn) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = zeros(Nn,Nn);
B = zeros(Nn,Nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assemble all triangular elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_coordinates_x = zeros(Ne,3);
node_coordinates_y = zeros(Ne,3);
delta = zeros(Ne,1);
a = zeros(Ne,3);
b = zeros(Ne,3);
c = zeros(Ne,3);
A_element = zeros(3);
B_element = zeros(3);

for ee=1:1:Ne      %looping through all the elements
    for nn=1:1:3   %looping through nodes of each element
        node_coordinates_x(ee,nn) = Points(ConnectivityArray(ee,nn),1);  %getting node x-coordinates
        node_coordinates_y(ee,nn) = Points(ConnectivityArray(ee,nn),2);  %getting node y-coordinates 
    end
    
    a(ee,1) = node_coordinates_x(ee,2)*node_coordinates_y(ee,3)-node_coordinates_y(ee,2)*node_coordinates_x(ee,3);
    a(ee,2) = node_coordinates_x(ee,3)*node_coordinates_y(ee,1)-node_coordinates_y(ee,3)*node_coordinates_x(ee,1);
    a(ee,3) = node_coordinates_x(ee,1)*node_coordinates_y(ee,2)-node_coordinates_y(ee,1)*node_coordinates_x(ee,2);
    
    b(ee,1) = node_coordinates_y(ee,2)-node_coordinates_y(ee,3);
    b(ee,2) = node_coordinates_y(ee,3)-node_coordinates_y(ee,1);
    b(ee,3) = node_coordinates_y(ee,1)-node_coordinates_y(ee,2);
    
    c(ee,1) = node_coordinates_x(ee,3)-node_coordinates_x(ee,2);
    c(ee,2) = node_coordinates_x(ee,1)-node_coordinates_x(ee,3);
    c(ee,3) = node_coordinates_x(ee,2)-node_coordinates_x(ee,1);
    
    delta(ee) = 0.5*(b(ee,1)*c(ee,2)-b(ee,2)*c(ee,1));
    
    for rr_element=1:1:3   %looping through rows of the matrix A_element and B_element
        rr = ConnectivityArray(ee,rr_element);  %row location in the ugmented matrix A and B
        for cc_element = 1:1:3    %looping through columns of the matrix A_element and B-element
            cc = ConnectivityArray(ee,cc_element);  %column location in the ugmented matrix A and B
            A_element(rr_element,cc_element) = (0.25/delta(ee))*(b(ee,rr_element)*b(ee,cc_element)+c(ee,rr_element)*c(ee,cc_element));
            if rr_element == cc_element
                B_element(rr_element,cc_element) = (delta(ee)/12)*(1+1);
            else
                B_element(rr_element,cc_element) = (delta(ee)/12);
            end
        A(rr,cc) = A(rr,cc) + A_element(rr_element,cc_element);
        B(rr,cc) = B(rr,cc) + B_element(rr_element,cc_element);
        end
    end 
    
end

clear node_coordinates_x;
clear node_coordinates_y;
clear A_element;
clear B_element;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Imposing Dirichlet BCs in matrices A & B %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for node_counter=1:1:Nn
    if  Points(node_counter,1)==0 || Points(node_counter,1)==a_RW || Points(node_counter,2)==0 || Points(node_counter,2)==b_RW
        A(node_counter,:)=0;
        A(:,node_counter)=0;
        A(node_counter,node_counter)=0;
        B(node_counter,:)=0;
        B(:,node_counter)=0;
        B(node_counter,node_counter)=1;
    else
    end
end
clear node_counter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Solving the Eigen-value Problem %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k0 = 2*pi.*freq./c0;
TEMP = A*inv(B);
[Ez,D] = eig(TEMP);

kc_TM = sqrt(diag(D));


%%%%%Claculating and Plotting kz_TM vs frequency %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(2);
% for eigen_value_counter = 1:1:size(kc_TM)
%     if kc_TM(eigen_value_counter)==0
%         continue;
%     else
%         kz_TM = sqrt((w.^2).*(u0*eps0)-(kc_TM(eigen_value_counter))^2);
%         plot(freq/1e9,kz_TM);
%         hold on;
%     end
% end
% 
% xlabel('Freq (GHz)');
% ylabel('kz (TM modes)');
% grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Claculating and Plotting total E-Field %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TM11 mode %%%
counter2 = 1;
for counter=1:1:size(kc_TM)
    if kc_TM(counter)==0
        continue;
    else
        kc_TM_shrink(counter2)=kc_TM(counter);
        counter2 = counter2+1;
    end
end

for counter=1:1:size(kc_TM)
    if kc_TM(counter)==min(kc_TM_shrink)
        kc_TM_min_index = counter;
        break;
    else
        continue;
    end
end

kc_TM11 = sqrt(min(kc_TM_shrink));
kz_TM11 = sqrt(w0^2*u0*eps0-kc_TM11^2);


%%%%First it is useful to arrange the Ez eigen vector into an array similar
%%%%to the connectivity array, in which each row represents an element and each
%%%%column represents a local node in the corresponding element, the reason
%%%%for this is that the interpolation happens within each element

Ez_TM11 = zeros(Ne,3);
Ex_TM11 = zeros(Ne,3);
Ey_TM11 = zeros(Ne,3);
Hx_TM11 = zeros(Ne,3);
Hy_TM11 = zeros(Ne,3);
for e_counter = 1:1:Ne
    for n_counter = 1:1:3
        Ez_TM11(e_counter,n_counter)=Ez(ConnectivityArray(e_counter,n_counter),kc_TM_min_index);     
    end
end

%%%Now we are ready to find the rest of field components Ex,Ey,Hx, and Hy

for e_counter = 1:1:Ne
    for n_counter = 1:1:3
        Ex_TM11(e_counter,n_counter) = ((-1i.*kz_TM11)/(2*delta(e_counter)*kc_TM11^2))*(b(e_counter,1)*Ez_TM11(e_counter,1)+b(e_counter,2)*Ez_TM11(e_counter,2)+b(e_counter,3)*Ez_TM11(e_counter,3));
        Ey_TM11(e_counter,n_counter) = ((-1i.*kz_TM11)/(2*delta(e_counter)*kc_TM11^2))*(c(e_counter,1)*Ez_TM11(e_counter,1)+c(e_counter,2)*Ez_TM11(e_counter,2)+c(e_counter,3)*Ez_TM11(e_counter,3));
        Hx_TM11(e_counter,n_counter) = ((1i.*w0*eps0)/(2*delta(e_counter)*kc_TM11^2))*(c(e_counter,1)*Ez_TM11(e_counter,1)+c(e_counter,2)*Ez_TM11(e_counter,2)+c(e_counter,3)*Ez_TM11(e_counter,3));
        Hy_TM11(e_counter,n_counter) = ((-1i.*w0*eps0)/(2*delta(e_counter)*kc_TM11^2))*(b(e_counter,1)*Ez_TM11(e_counter,1)+b(e_counter,2)*Ez_TM11(e_counter,2)+b(e_counter,3)*Ez_TM11(e_counter,3));
    end
end

%%%%Now we are ready to clculate the total field at each node


Etot_TM11 = sqrt((abs(Ex_TM11)).^2+(abs(Ey_TM11)).^2+(abs(Ez_TM11)).^2);
Htot_TM11 = sqrt((abs(Hx_TM11)).^2+(abs(Hy_TM11)).^2);

%%%%Now we need to arrange the fields arrays in such a way to be combatible
%%%%with the cross section of the waveguide

Etot_TM11_arranged = zeros(Ny,Nx);
Htot_TM11_arranged = zeros(Ny,Nx);

Ex_TM11_arranged = zeros(Ny,Nx);
Hx_TM11_arranged = zeros(Ny,Nx);
Ey_TM11_arranged = zeros(Ny,Nx);
Hy_TM11_arranged = zeros(Ny,Nx);

for x_counter=1:1:Nx
    for y_counter=1:1:Ny
        x_coordinate = (x_counter-1)*dx;
        y_coordinate = (y_counter-1)*dy;
        for global_node_counter=1:1:Nn
            if abs(Points(global_node_counter,1)-x_coordinate)<=0.001 &&  abs(Points(global_node_counter,2)-y_coordinate)<=0.001
                break;
            end
        end
        
        for element_counter=1:1:Ne            
            if ConnectivityArray(element_counter,1)==global_node_counter
                local_node_counter=1;
                break;
            elseif ConnectivityArray(element_counter,2)==global_node_counter
                local_node_counter=2;
                break;
            elseif ConnectivityArray(element_counter,3)==global_node_counter
                local_node_counter=3;
                break;
            end
        end       
        Etot_TM11_arranged(y_counter,x_counter) = Etot_TM11(element_counter,local_node_counter);
        Htot_TM11_arranged(y_counter,x_counter) = Htot_TM11(element_counter,local_node_counter);
        
        Ex_TM11_arranged(y_counter,x_counter) = Ex_TM11(element_counter,local_node_counter);
        Hx_TM11_arranged(y_counter,x_counter) = Hx_TM11(element_counter,local_node_counter);
        Ey_TM11_arranged(y_counter,x_counter) = Ey_TM11(element_counter,local_node_counter);
        Hy_TM11_arranged(y_counter,x_counter) = Hy_TM11(element_counter,local_node_counter);
    end
end


%% Plot the Field distribution
% normalize field value based on max(Etot)
Emax=max(max(Etot_TM11_arranged));
Ex=Ex_TM11_arranged/1.0i/Emax;
Ey=Ey_TM11_arranged/1.0i/Emax;
Hx=Hx_TM11_arranged/1.0i/Emax;
Hy=Hy_TM11_arranged/1.0i/Emax;

% figure(3);
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
end