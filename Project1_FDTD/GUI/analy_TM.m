function [Ex,Ey,Ez,Hx,Hy,Hz]=analy_TM(B,m,a,omega,mu,epsilon,t_E,t_H,Y_Ey,Y_Ez,Y_Hx,Z_Ey,Z_Ez,Z_Hx)
% Function analy_TM provides the analytic values of TM mode
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

Ey = B * sin(kc * Z_Ey) .* cos(omega * t_E - beta * Y_Ey);
Ez =  beta / kc * B * cos(kc * Z_Ez) .* sin(omega * t_E - beta * Y_Ez);
Hx = omega * epsilon / kc * B * cos(kc * Z_Hx) .* sin(omega * t_H - beta * Y_Hx);

Hy=0;
Hz=0;
Ex=0;
end