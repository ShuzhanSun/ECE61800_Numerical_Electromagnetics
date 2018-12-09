function [Ex,Ey,Ez,Hx,Hy,Hz] = ParallelPlateAnalyOnePointTM(B, m, a, omega, t, mu, epsilon, y, z)
%m is mode
v = 1 / sqrt(epsilon * mu);
kc = m * pi / a;
omegac = m * pi * v / a;
beta = omega * sqrt(epsilon * mu) * sqrt(1 - (omegac/omega)^2);

Ey = B * sin(kc * z) .* cos(omega * t - beta * y);
Ez =  beta / kc * B * cos(kc * z) .* sin(omega * t - beta * y);
Hx = omega * epsilon / kc * B * cos(kc * z) .* sin(omega * t - beta * y);

Hy=0;
Hz=0;
Ex=0;
end