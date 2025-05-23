function [ aon ] = AonGSurfHFinPInfOrd1( omega0 , k0 , epsilon , M , B , ~ , H , ~ , ~ )
% AonGSurfHFinPInfOrd1 returns the dimensional onset acceleration for the first-order surfactant-covered but no diffusion finite-depth case.

VS = 1 ...
    + sqrt(2) * B ...
    + B^2 ...
    - M * sqrt(2) ...
    + M^2 ...
    ;



aon = ( omega0^2 / k0 ) * ( ...
    epsilon * ( ...
        sqrt(2) * ( ...
            VS * csch(H)^2 ...
            + coth(H)^2 * ( ...
                VS ...
                - 1 ...
                + sqrt(2) * M ...
            ) ...
        )/( ...
            VS ...
        ) ...
    ) ...
    );


end


