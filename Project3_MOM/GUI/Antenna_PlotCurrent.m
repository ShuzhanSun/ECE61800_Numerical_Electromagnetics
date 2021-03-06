function Antenna_PlotCurrent(handles)

% Constitutive parameters
mu0 = 4*pi*1.0e-7;  
epsi0 = 8.854e-12; 
epsi_r = 1; % relative permittivity
epsi = epsi0*epsi_r;
Z0 = sqrt(mu0/epsi);

% Segment size
diameter = str2double(get(handles.edit_2a,'String'));%0.01; % thickness (m)
if diameter <= 0
    msg = 'Error. 2a should be positve!';
    errordlg(msg);
    error(msg);
end
contents = (get(handles.popup_LBy2a,'String'));
LByDiam = str2double(contents{get(handles.popup_LBy2a,'Value')});

L = LByDiam*diameter; % length (m)

NSeg = str2double(get(handles.edit_NSeg,'String'));%49;
if NSeg-floor(NSeg)~=0.0  || NSeg<=0
    msg = 'Error. Number of segments should be positve integer!';
    errordlg(msg);
    error(msg);
end

dz = L/NSeg;
z = -L/2:dz:L/2;

% Define the delta-gap source @ mFed
V0 = 1; % voltage (V)
if get(handles.popup_Pattern,'Value') == 1
    mFed = round(NSeg/2); % center-fed
else
    mFed = 1;  % end-fed
end
V = zeros(NSeg-1,1);
V(mFed)= V0;
Zmatrix = zeros(NSeg-1,NSeg-1);

%%
LByLamb = str2double(get(handles.edit_LByLamb,'String'));
if LByLamb <= 0
    msg = 'Error. L/\lambda should be positve!';
    errordlg(msg);
    error(msg);
end

lambda = L/LByLamb;
k0 = 2*pi/lambda;
if dz > lambda/5
    msg = 'Error. Segment Length is too Large to Resolve Current Wavelength!';
    errordlg(msg);
    error(msg);
end

% Define 3-D Green's function
rho = @(z1,z2) sqrt((z1-z2).^2+diameter^2/4);
Green0 = @(z1,z2) exp(-1.0i*k0*rho(z1,z2))/4/pi./rho(z1,z2);

for m = 1:NSeg-1
    for n = m:NSeg-1
        % Define basis function
        Basism = @(z1) 1-abs((z1-z(m+1))./(z(m+2)-z(m+1)));
        Basisn = @(z2) 1-abs((z2-z(n+1))./(z(n+2)-z(n+1)));
        dBasism = @(z1) 1/dz*(heaviside(z(m+1)-z1)-heaviside(z1-z(m+1)));
        dBasisn = @(z2) 1/dz*(heaviside(z(n+1)-z2)-heaviside(z2-z(n+1)));
        inteFun1 = @(z1,z2) Basism(z1).*Basisn(z2).*Green0(z1,z2);
        inteFun2 = @(z1,z2) dBasism(z1).*dBasisn(z2).*Green0(z1,z2);

        integrand = @(z1,z2) 1.0i*k0*Z0*inteFun1(z1,z2)-1.0i*Z0/k0*inteFun2(z1,z2);
        Zmatrix(m,n) = integral2(integrand,z(m),z(m+2),z(n),z(n+2));
        Zmatrix(n,m) = Zmatrix(m,n);
    end
end

% current and admittance
I=Zmatrix\V;
I=[0,I',0];

axes(handles.axes_Current);
set(gca,'fontsize',16);
xlabel('Position z ( m ) ');%,'FontSize',16);
hold on;

yyaxis left;
plot(z,abs(I),'Color',[0,0.45,0.74],'LineWidth',3);
ylabel('Current magnitutude ( A )');%,'FontSize',16);

yyaxis right;
plot(z,angle(I),'Color',[0.85,0.33,0.1],'LineWidth',3);
ylabel('Current phase angle ( rad )');%,'FontSize',16);

legend('abs( I )','angle( I )');
drawnow;

end