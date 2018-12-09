function[DT]=MeshGen(wd,ht,Nx,Ny)
a_RW=wd;
b_RW=ht;
stepX=Nx;
stepY=Ny;
freq = 0.1e9:0.01e9:20e9;
w = 2.*pi.*freq;
% a_RW = 5e-2;
% b_RW = a_RW/2;
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

