function [wavenumber_z,any_K_z]=TM_PropCons(wd,ht,Nx,Ny)
freq = 0.1e9:0.01e9:10e9;
w = 2.*pi.*freq;
a_RW = wd;
b_RW = ht;
stepX=Nx-1;
stepY=Ny-1;
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

A = zeros(Nn);
B = zeros(Nn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Assemble all triangular elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ee=1:1:Ne      %looping through all the elements
    for nn=1:1:3   %looping through nodes of each element
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
    
    delta(ee) = 0.5*(b(1)*c(2)-b(2)*c(1));
    
    
    
    
    for rr_element=1:1:3   %looping through rows of the matrix A_element and B_element
        rr = ConnectivityArray(ee,rr_element);  %row location in the ugmented matrix A and B
        for cc_element = 1:1:3    %looping through columns of the matrix A_element and B-element
            cc = ConnectivityArray(ee,cc_element);  %column location in the ugmented matrix A and B
            A_element(rr_element,cc_element) = (0.25/delta(ee))*(b(rr_element)*b(cc_element)+c(rr_element)*c(cc_element));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Imposing Dirichlet BCs in matrices A & B %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for counter=1:1:Nn
    if  Points(counter,1)==0 || Points(counter,1)==a_RW || Points(counter,2)==0 || Points(counter,2)==b_RW
        A(counter,:)=0;
        A(:,counter)=0;
        A(counter,counter)=0;
        B(counter,:)=0;
        B(:,counter)=0;
        B(counter,counter)=0;
    else
    end
end

%% remove the rows & columns corresponding to Dirichlet BC
A_shrink=A;
B_shrink=B;
% Remove zero rows in K
A_shrink( all(~A_shrink,2), : ) = [];
B_shrink( all(~B_shrink,2), : ) = [];
% Remove zero columns
A_shrink( :, all(~A_shrink,1) ) = [];
B_shrink( :, all(~B_shrink,1) ) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Solving the Eigen-value Problem %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavenumber_cutoff = 1./(sqrt(eig(inv(A)*B)));
wavenumber_cutoff = sqrt(eig(A_shrink*inv(B_shrink)));
%[Ez,D] = eig(A*inv(B));
% sort kc in ascending order
wavenumber_cutoff=sort(wavenumber_cutoff); 


m=1;
% figure(2);
% for m = 1:1:size(wavenumber_cutoff)
    wavenumber_z = sqrt((w.^2).*(4*pi*1e-7).*(1e-9/(36*pi))-wavenumber_cutoff(m).^2);
    any_K_z= sqrt((w.^2).*(4*pi*1e-7).*(1e-9/(36*pi))-(((1*pi)./a_RW).^2)-(((1*pi)./b_RW).^2));
%     plot(freq/1e9,wavenumber_z);
%     hold on;
%     grid on;
% end
