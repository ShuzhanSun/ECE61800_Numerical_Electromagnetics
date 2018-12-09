function FDTD_3D_Waveguide_For_Plot_TE(mu_r, epsi_r, f_max, W, L, d, m, hfdtd, hanaly, ax3, initWithAnaly)
% This function "uses 3D FDTD algorithm to solve a wave propagation problem
% in a parallel plate waveguide". 

% This code adopts such BCs: top & bottom are PEC, front and back are PMC,
% left and right are known analytic field values.

mu0=4*pi*1.0e-7;  
epsilon0=8.854e-12; 

Nt=1e3; % number of time marching steps

mu=mu0*mu_r;
epsi=epsilon0*epsi_r;
c_wave=1/sqrt(mu*epsi);
omega = 2*pi*f_max;

B = 1e-5;

%% Grid size determined by the signal min wavelength 
dy=c_wave/f_max/10; dy=dy/2;
dx=dy;
dz=dy;

Nx=round(W/dx);
Ny=round(L/dy);
Nz=round(d/dz);
if Nx<10 || Ny<10 || Nz<10
    msg='Error. Grid Size is too Large to Resolve the Structure!';
    error(msg);
end

% Time step satisfying stability
dt=1/c_wave/sqrt(dx^(-2)+dy^(-2)+dz^(-2));
dt=dt/2;

%% Core 3D FDTD code
% Relevant equations and discretizations can be found in Jian-Ming Jin's
% book "Theory and Computation of Electromagnetic Fields", pp 309-314.

% Initialize all fields as 0
Ex=zeros(Nx,Ny+1,Nz+1);
Ey=zeros(Nx-1,Ny,Nz+1);
Ez=zeros(Nx-1,Ny+1,Nz);
Hx=zeros(Nx-1,Ny,Nz);
Hy=zeros(Nx,Ny+1,Nz);
Hz=zeros(Nx,Ny,Nz+1);

% The x-y-z grids of E&H 
[X_Ex,Y_Ex,Z_Ex]=ndgrid(dx/2:dx:(Nx-.5)*dx,0:dy:Ny*dy,0:dz:Nz*dz);
[X_Ey,Y_Ey,Z_Ey]=ndgrid(dx:dx:(Nx-1)*dx,dy/2:dy:(Ny-.5)*dy,0:dz:Nz*dz);
[X_Ez,Y_Ez,Z_Ez]=ndgrid(dx:dx:(Nx-1)*dx,0:dy:Ny*dy,dz/2:dz:(Nz-.5)*dz);
[X_Hx,Y_Hx,Z_Hx]=ndgrid(dx:dx:(Nx-1)*dx,dy/2:dy:(Ny-.5)*dy,dz/2:dz:(Nz-.5)*dz);
[X_Hy,Y_Hy,Z_Hy]=ndgrid(dx/2:dx:(Nx-.5)*dx,0:dy:Ny*dy,dz/2:dz:(Nz-.5)*dz);
[X_Hz,Y_Hz,Z_Hz]=ndgrid(dx/2:dx:(Nx-.5)*dx,dy/2:dy:(Ny-.5)*dy,0:dz:Nz*dz);

% Initialize all fields as the analytic value

% TE mode
if(initWithAnaly)
    [Ex,~,~,~,Hy,Hz]=analy_TE(B,m,d,omega,mu,epsi,0,-dt/2,Y_Ex,Y_Hy,Y_Hz,Z_Ex,Z_Hy,Z_Hz);
end

tic

for n=1:Nt % time marching
    
    % Time instant of E & H fields, offset by dt/2
    t_E=(n-1)*dt;
    t_H=t_E-dt/2;
    
    % Enforce BCs at every time step. PEC and PMC are naturally enforced.
    % So, only the left and right BCs need to be given here.
    
    % TE mode BC
    [Ex(:,1,:),Ey(:,1,:),Ez(:,1,:),Hx(:,1,:),Hy(:,1,:),Hz(:,1,:)]=...
        analy_TE(B,m,d,omega,mu,epsi,t_E,t_H,0,0,dy/2,Z_Ex(:,1,:),Z_Hy(:,1,:),Z_Hz(:,1,:));
    [Ex(:,Ny+1,:),Ey(:,Ny,:),Ez(:,Ny+1,:),Hx(:,Ny,:),Hy(:,Ny+1,:),Hz(:,Ny,:)]=...
        analy_TE(B,m,d,omega,mu,epsi,t_E,t_H,Ny*dy,Ny*dy,(Ny-.5)*dy,Z_Ex(:,Ny+1,:),Z_Hy(:,Ny+1,:),Z_Hz(:,Ny,:));

    %% Renew the fields at next time step
    % Renew H first
    Hx(:,2:Ny-1,:)=Hx(:,2:Ny-1,:)-dt/mu/dy*(Ez(:,3:Ny,:)-Ez(:,2:Ny-1,:))...
        +dt/mu/dz*(Ey(:,2:Ny-1,2:Nz+1)-Ey(:,2:Ny-1,1:Nz));
    Hy(2:Nx-1,2:Ny,:)=Hy(2:Nx-1,2:Ny,:)-dt/mu/dz*(Ex(2:Nx-1,2:Ny,2:Nz+1)...
        -Ex(2:Nx-1,2:Ny,1:Nz))+dt/mu/dx*(Ez(2:Nx-1,2:Ny,:)-Ez(1:Nx-2,2:Ny,:));
    Hz(2:Nx-1,2:Ny-1,2:Nz)=Hz(2:Nx-1,2:Ny-1,2:Nz)-dt/mu/dx*(Ey(2:Nx-1,2:Ny-1,2:Nz)-Ey(1:Nx-2,2:Ny-1,2:Nz))...
        +dt/mu/dy*(Ex(2:Nx-1,3:Ny,2:Nz)-Ex(2:Nx-1,2:Ny-1,2:Nz));
    % Renew E next
    Ex(2:Nx-1,2:Ny,2:Nz)=Ex(2:Nx-1,2:Ny,2:Nz)+dt/epsi/dy*(Hz(2:Nx-1,2:Ny,2:Nz)...
        -Hz(2:Nx-1,1:Ny-1,2:Nz))-dt/epsi/dz*(Hy(2:Nx-1,2:Ny,2:Nz)-Hy(2:Nx-1,2:Ny,1:Nz-1));
    Ey(:,2:Ny-1,2:Nz)=Ey(:,2:Ny-1,2:Nz)+dt/epsi/dz*(Hx(:,2:Ny-1,2:Nz)-Hx(:,2:Ny-1,1:Nz-1))...
        -dt/epsi/dx*(Hz(2:Nx,2:Ny-1,2:Nz)-Hz(1:Nx-1,2:Ny-1,2:Nz));
    Ez(:,2:Ny,:)=Ez(:,2:Ny,:)+dt/epsi/dx*(Hy(2:Nx,2:Ny,:)-Hy(1:Nx-1,2:Ny,:))...
        -dt/epsi/dy*(Hx(:,2:Ny,:)-Hx(:,1:Ny-1,:));
    
    %% Plot at every time instant
    addpoints(hfdtd,t_H,Hy(10,3,10));
    drawnow;
    
    %% compute the relative error
    [~,~,~,~,HyA,~] = ParallelPlateAnalyOnePointTE(B, m, d, omega, t_H, mu, epsi, 2*dy, 19/2*dz);
    % show analytical result:
    addpoints(hanaly,t_H,HyA);
    drawnow;
    
    imagesc(ax3, squeeze(Hz(round(Nx/2),:,:)))
    
end
toc

end