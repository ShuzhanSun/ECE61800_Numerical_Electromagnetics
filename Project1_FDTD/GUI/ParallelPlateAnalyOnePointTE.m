function [Ex,Ey,Ez,Hx,Hy,Hz] = ParallelPlateAnalyOnePointTE(B, m, a, omega, t, mu, epsilon, y, z)
%m is mode
v = 1 / sqrt(epsilon * mu);
kc = m * pi / a;
omegac = m * pi * v / a;
beta = omega * sqrt(epsilon * mu) * sqrt(1 - (omegac/omega)^2);

Hy = B * cos(kc * z) * cos(omega * t - beta * y);
Hz = - beta / kc * B * sin(kc * z) * sin(omega * t - beta * y);
Ex = omega * mu / kc * B * sin(kc * z) * sin(omega * t - beta * y);
Ey = 0;
Ez = 0;
Hx = 0;
end