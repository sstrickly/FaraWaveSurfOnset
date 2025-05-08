function [ f ] = KelvinDispersionRelationTester( g , omega0 , sigma0 , rho , H0 , k0G )
%KELVINDISPERSIONRELATIONTESTER is for plugging in the guessed wavenumber
% to see if the RHS nears 1.
% 
% 1 == tanh( k0*H0 ) * ( g*k0/omega0^2 + sigma0*k0^3/(rho*omega0^2)

f = tanh( k0G.*H0 ) .* ( g.*k0G./omega0^2 + (sigma0./(rho.*omega0^2)).*k0G.^3 );


end

