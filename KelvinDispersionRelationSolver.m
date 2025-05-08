function [ k0 , k0Err ] = KelvinDispersionRelationSolver( varargin )
%KELVINDISPERSIONRELATIONSOLVER to numerically solve the finite depth Kelvin dispersion relation so that the length scale l0 = 1/k0 is well defined.
%
% The finite depth Kelvin dispersion relation is:
% 1 == tanh( k0*H0 ) * ( g*k0/omega0^2 + sigma0*k0^3/(rho*omega0^2)
% This function uses the Levenberg-Marquardt minimization algorithm to numerically find the optimal  k0  .
% We minimize the error function (least squares fitting):  (YVEC - yVec)^2  where  YVEC = xVec * tanh( k0*H0 ) * ( g*k0/omega0^2 + sigma0*k0^3/(rho*omega0^2)  AND  yVec = 1  .
% Effectively, we are using the "data"  (x,y) = (1,1)  to fit the function  y(x) = x * tanh( k0*H0 ) * ( g*k0/omega0^2 + sigma0*k0^3/(rho*omega0^2)  where the only fitting parameter is  k0  .
% The only correct fit to this is the numerical solution to the Kelvin dispersion relation (above).
% 
% An example of what to run:  
%
% g = 9.8; % m/s
% omega0 = pi*30; % rad/s
% sigma0 = 0.07; % N/m
% rho = 1e3; % kg/m^3
% H0 = 7.9e-3; % m
% k0G = 1; % 1/m
% MeasurementSigma = 1e-6; % 1/m
% [ k0 , k0Err ] = KelvinDispersionRelationSolver( g , omega0 , sigma0 , rho , H0 , k0G , MeasurementSigma )

switch nargin
    case 5
        g = varargin{1};
        omega0 = varargin{2};
        sigma0 = varargin{3};
        rho = varargin{4};
        H0 = varargin{5};
        k0G = 1; % 1/m
        MeasurementPercentSigma = 1e-6;
    case 6
        g = varargin{1};
        omega0 = varargin{2};
        sigma0 = varargin{3};
        rho = varargin{4};
        H0 = varargin{5};
        k0G = varargin{6};
        MeasurementPercentSigma = 1e-6;
    case 7
        g = varargin{1};
        omega0 = varargin{2};
        sigma0 = varargin{3};
        rho = varargin{4};
        H0 = varargin{5};
        k0G = varargin{6};
        MeasurementPercentSigma = varargin{7};
    otherwise
        k0 = Nan;
        k0Err = Nan;
        return
end

%% Levenberg-Marquardt algorithm
    %% Setup
    % In this fitting algorithm, the constants (cVec) are:    g , omega0 , sigma0 , rho , H0
    %
    % In this fitting algorithm, the fitting parameters (bVec) are:    k0
    %
    % In this fitting algorithm, the variables (xVec) are:
    %
    % In this fitting algorithm, the fit data (yVec) is:
    %
    % In the code,  YVEC  is the projected fit data for the present state of the parameters.
    % This method will perturb the fit parameters (vVec), determine how the goodness of fit is affected, and then determine the next best guess for the parameters.
    
    % Setup bVec and cVec
    bVec = [ k0G ];
    bVec = bVec(:);
    cVec = [ g , omega0 , sigma0 , rho , H0 ];
    
    % Measurement error
    % derived from the step size of the increment
    W = 1/(MeasurementPercentSigma^2); % m^2
    
    % Setup xVec and yVec 
    xVec = [ 1 ];
    yVec = [ 1 ];
    yVec = yVec(:);




    %% Execute LM method
    Lcnt = 1;
    LcntMax = 31;
    Llogscale = logspace(-10,20,LcntMax);
    cnt = 1;
    cond = true;
    lambda = Llogscale(Lcnt);
    DirVary = 100 * eps/MeasurementPercentSigma;
    
    YVEC2 = LMBuildYKelvin(xVec,bVec,cVec);
    while(cond)
        % Obtain Y Values for current bVec
        YVEC = YVEC2;
        
        % Test to see how far away we are.  These numbers should gradually decrease.
%         sum((YVEC - yVec).^2)/length(yVec)
        
        % Build Jacobian
        J = zeros([length(YVEC),length(bVec)]);
        for j = 1:length(bVec)
            dbVec = zeros(size(bVec));
            dbVec(j) = DirVary;
            YVEC2 = LMBuildYKelvin(xVec,bVec+dbVec,cVec);
            J(:,j) = (YVEC2 - YVEC)/DirVary;
        end
        
        % Obtain delta
        D = diag(J'*J);
        Diag = zeros(length(D),length(D));
        for j = 1:length(D)
            Diag(j,j) = D(j);
        end
        Lcnt = 1;
        lambda = Llogscale(Lcnt);
        delta = ( J'*J + lambda*Diag )^(-1)*((J')*(yVec - YVEC));
        while( ( Lcnt <= LcntMax ) && ( sum((YVEC - yVec).^2)/length(yVec) < sum((LMBuildYKelvin(xVec,bVec+delta,cVec) - yVec).^2)/length(yVec) ) )
            Lcnt = Lcnt + 1;
            lambda = Llogscale(Lcnt);
            delta = ( J'*J + lambda*Diag )^(-1)*((J')*(yVec - YVEC));
        end
%         Lcnt,
%         delta,
%         if(isnan(delta))
%             delta,
%         end
        
        % Obtain new bVec
        bVec = bVec + delta;
        YVEC2 = LMBuildYKelvin(xVec,bVec,cVec);

        if( (sum((YVEC2 - yVec).^2) < MeasurementPercentSigma^2) || isnan(delta) )
            cond = false;
        else
            cnt = cnt+1;
        end
    end
    
    
    
    %% Evaluate the goodness of fit
    DoF = length(yVec) - length(bVec) + 1;
    ChiSquared = (YVEC2 - yVec)' * W * (YVEC2 - yVec);
    ReducedChiSquared = ChiSquared/DoF;

    %% Find parametres, calculate errors, and calculate percent errors
    ERRORs = sqrt( diag( inv(J'*W*J) ) );  % Error estimates are based on the paper at http://people.duke.edu/~hpgavin/ce281/lm.pdf

    k0 = bVec(1);

    k0Err = ERRORs(1);

end

