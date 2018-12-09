function [Ex,Ey,Ez,Hx,Hy,Hz]=analy_TE(B,m,a,omega,mu,epsilon,t_E,t_H,Y_Ex,Y_Hy,Y_Hz,Z_Ex,Z_Hy,Z_Hz)
% Function analy_TE provides the analytic values of TE mode
% m is mode number
% a is distance between 2 plates
% omega is signal frequency

v = 1 / sqrt(epsilon * mu);
kc = m * pi / a;
omegac = m * pi * v / a;
beta = omega * sqrt(epsilon * mu) * sqrt(1 - (omegac/omega)^2);
if omega<=omegac
    msg='Error. Signal Frequency < Cutoff Frequency, Unable to Transfer';
    error(msg);
end

Hy = B * cos(kc * Z_Hy) .* cos(omega * t_H - beta * Y_Hy);
Hz = - beta / kc * B * sin(kc * Z_Hz) .* sin(omega * t_H - beta * Y_Hz);
Ex = omega * mu / kc * B * sin(kc * Z_Ex) .* sin(omega * t_E - beta * Y_Ex);

Ey=0;
Ez=0;
Hx=0;
end