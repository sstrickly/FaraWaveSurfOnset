function [ YVEC ] = LMBuildYKelvin( xVec , bVec , cVec )
%LMBUILDYKelvin is an assistance function for running a Levenburg
% Marquardt algorithm to best fit wave number for the finite-depth Kelvin
% dispersion relaiton.

% yVec = 1;
% xVec = [ 1 ];
% bVec = [ k0g ];
% cVec = [ g , omega0 , sigma0 , rho , H0 ];

k0G = bVec(1);

g = cVec(1);
omega0 = cVec(2);
sigma0 = cVec(3);
rho = cVec(4);
H0 = cVec(5);

YVEC = xVec.*tanh(k0G.*H0).*k0G.*( (g/omega0^2) + (sigma0./(rho.*omega0.^2)).*k0G.^2 );

YVEC = YVEC';


end

