% This Master Example Script file gives an example of how to use the onset acceleration functions.
% At the end, the script generates figures for the onset acceleration vs Marangoni & Boussinesq numbers
% as well as a 3-d plot of the wave speed vs depth and frequency (i.e. the dispersion relation).


clear all;


%% Test the Aon functions to ensure that they all work correctly

omega0 = pi*30; % rad/s
g = 9.8; % m/s
rho = 1e3; % kg/m^3
mu = 1e-3; % kg/m/s
Diffusion = 1e-9; % m/s^2    % Negligible but finite diffusion
sigma0 = 0.07; % N/m
H0 = 7.9e-3; % m

% Solve for l0 and k0
k0 = KelvinDispersionRelationSolver( g , omega0 , sigma0 , rho , H0 );
l0 = 1/k0;

G0 = 1e-6; % kg/m^2
dsigdGamma = -4.8500e+04 ; % N m / kg
E0 = - G0 * dsigdGamma ; % N / m 
Mu = 1e-5 ; % kg / s
% Typical values of the surface elasticity and viscosity:
% E0 = [ 4.68e-4 , 0.007 , 0.0937 ]; % N / m
% M0 = [ 0 , 5.0948e-6 , 5.145e-5 ]; % kg / s

G = g / ( l0 * omega0^2 );
Sigma = sigma0 / ( rho * omega0^2 * l0^3 );
M = E0 * sqrt( 1 / ( mu * rho * omega0^3 * l0^4 ) );
B = Mu * sqrt( 1 / ( mu * rho * omega0 * l0^4 ) );
Di = Diffusion / ( omega0 * l0^2 ); % Recirpocal peclet number
epsilon = sqrt( mu / ( rho * omega0 * l0^2 ) );
H = H0*k0; % dimensionless fluid depth





% Infinite depth no surfactant (i.e. Chen & Vinals limit)
Aon0SurfHInfOrd1( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G ) 
Aon0SurfHInfOrd2( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )
Aon0SurfHInfOrd3( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )
Aon0SurfHInfOrd4( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )
Aon0SurfHInfOrd5( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G ) 

% Finite depth no surfactant
Aon0SurfHFinOrd1( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G ) 
Aon0SurfHFinOrd2( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )

% Infinite depth with surfactant
AonGSurfHInfOrd1( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )
AonGSurfHInfOrd2( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )

% Finite depth with surfactant
AonGSurfHFinOrd1( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )
AonGSurfHFinOrd2( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )

% Finite depth with surfactant No Diffusion
AonGSurfHFinPInfOrd1( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )
AonGSurfHFinPInfOrd2( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G )





%% Clear parameters

clearvars;

%% 3-d Plot where the Marangoni and Boussinesq numbers are varried for a driving frequency of 120 Hz and fluid depth of 1.5e-3 m

% Constants
g = 9.8; % m/s
rho = 1e3; % kg/m^3
mu = 1e-3; % kg/m/s
sigma0 = 0.07; % N/m % a mean estimate but does not vary in this figure
H0 = 1.5e-3; % m;
omega = 2*pi*120; % rad/s
omega0 = omega/2; % rad/s

% Solve for the characteristic wave numbers
k0 = KelvinDispersionRelationSolver( g , omega0 , sigma0 , rho , H0 );
l0 = 1/k0;

MMin = 1e-3;
MMax = 1e+3;
BMin = 1e-4;
BMax = 1e+2;
AMin = 0;
AMax = 2.5e0;

Pe = 7.991e-2;
% Pe = [ 7.991e-2 , 7.991e-1 , 7.991e0 , 7.991e4 ]; % Typical values for the Peclet number.
X = logspace(log10(MMin),log10(MMax),50); % E0
Y = logspace(log10(BMin),log10(BMax),50); % M0

Z = zeros(length(X),length(Y));
for i = 1:length(X)
%     E0 = X(i); % Use if iterating over surface elasticity instead of Marangoni number.
    M = X(i);

    for j = 1:length(Y)
%         M0 = Y(j); % Use if iterating over surface viscosity instead of Boussinesq number.
        B = Y(j);
    
        % Calculate non-dim constants
        G = g / ( l0 * omega0^2 );
        Sigma = sigma0 / ( rho * omega0^2 * l0^3 );
        epsilon = sqrt( mu / ( rho * omega0 * l0^2 ) );
        H = H0*k0;
        Di = 1 / Pe ;
%         M = E0/sqrt(mu*rho*omega0^3*l0^4); % Use if iterating over surface elasticity instead of Marangoni number.
%         B = M0/sqrt(mu*rho*omega0*l0^4); % Use if iterating over surface viscosity instead of Boussinesq number.

        Z(j,i) = AonGSurfHFinOrd2( omega0 , k0 , epsilon , M , B , Di , H , Sigma , G ) / g; % (in unitys of g)
    end
    i,
end

C = Z
hFig = figure('units','normalized','outerposition',[0 0 0.60 0.60]); hold on;
surf( X , Y , Z , C );

ax = get(hFig,'CurrentAxes');
set(ax,'XScale','log');
set(ax,'YScale','log');
set(ax,'ZScale');
axis([MMin,MMax,BMin,BMax,AMin,AMax]);
box on;
xlh = xlabel('$\mathcal{M}^{\dagger}$','Interpreter','latex');
ylh = ylabel('$\mathcal{B}^{\dagger}$','Interpreter','latex');
zlabel('$a_{onset}$ (g)','Interpreter','latex')
az = -38.4; el = 28.0; view(az,el);

caxis([AMin AMax]);
colorbar()

close(hFig);





%% Clear parameters

clearvars;

%% 3-d Plot to find Gravity & Capillary regimes via wave speed when the driving frequency and fluid depth are varried
% Note the four regimes in the final diagram:
%   Hi H & Lo w ; gravity waves
%   Hi H & Hi w ; capillary waves
% 
%   Lo H & Lo w ; Depth-restricted gravity waves
%   Lo H & Hi w ; Depth-restricted capillary waves

% Constants
g = 9.8; % m/s
rho = 1e3; % kg/m^3
mu = 1e-3; % kg/m/s

E0 = 0;
M0 = 0;
Pe = 7.991*1e4;
sigma0 = 7e-2;
Diffusion = ( (pi*120)*(8e-4)^2 ) ./ Pe;



HMin = 1e-6 ;
HMax = 1e1 ;
fMin = 1e-2 ;
fMax = 1e4 ;

    omega = 2*pi*logspace(log10(fMin),log10(fMax),40); % rad/s
    X = logspace(log10(HMin),log10(HMax),40); % m
    Y = omega; % rad/s
    
    Z = zeros(length(X),length(Y));
    for j = 1:length(Y)
        omega0 = Y(j)/2;
        k0G = 1;
        
%         for i = 1:length(X)
        for i = length(X):-1:1
            % Starting at the deepest depths and moving to the shallower will help the k0tmp solver. 
            H0 = X(i);
            
            k0 = KelvinDispersionRelationSolver( g , omega0 , sigma0 , rho , H0 , k0G );
            k0G = k0;
            l0 = 1/k0;
            
            Z(i,j) = omega0/k0; %phase velocity in units of m/s
        end
        j,
    end
    
    ZMin = min(Z(:));
    ZMax = max(Z(:));
    
    Msk1 = isinf(Z) + isnan(Z);
    Z(logical(Msk1)) = NaN;
    
    C = log10(Z);
    hFig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    surf( X , Y/(2*pi) , Z' , C' );
    colormap(hot)
    
    ax = get(hFig,'CurrentAxes');
    set(ax,'XScale','log');
    set(ax,'YScale','log');
    set(ax,'ZScale','log');
    axis([HMin,HMax,fMin,fMax,ZMin,ZMax]);
    box on;
    xlabel('$H$ (m)','Interpreter','latex')
    ylabel('$\omega / 2 \pi$ (Hz)','Interpreter','latex')
    zlabel('$c$ (m/s)','Interpreter','latex')
    az = 70.8; el = 24.4; view(az,el);
    
    caxis(log10([ZMin ZMax]));
    htmp = colorbar('FontSize',15);
    ylabel( htmp , '$c$ (m/s)','Interpreter','latex')

%     az = -105.0; el = 24.4; view(az,el);

%     az = 90; el = 90; view(az,el);

close(hFig);

    






