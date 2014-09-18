function alpha = mgso4_alpha(S,T,P)

% MGSO4_ALPHA    Coefficient of thermal expansion of MgSO4 brines
%=========================================================================
%
% USAGE:  alpha = mgso4_alpha(S,T,P)
%
% DESCRIPTION:
%    Thermal expansion coefficient of magnesium sulfate brine using 
%    Vance and Brown data
%
% INPUT: 
%   S = salinity    [molal]
%   T = temperature [Kelvin]
%   P = pressure    [bars]
%
% OUTPUT:
%   alpha = Coefficient of thermal expansion  [1/K] 
%
% Definition:
%   alpha = - (1/rho) d(rho)/dT at constant P,S.
% 
% AUTHOR:  Jason Goodman (goodman_jason@wheatonma.edu)
%
% REFERENCES:  
% Marion et al (2005), Geochimica et Cosmochimica Acta, Vol. 69, No. 2, pp. 259?274
% Vance and Brown (2005), doi:10.1016/j.icarus.2005.06.005
% Vance and Brown (2011), pers. comm.
%=========================================================================

% CALLER: general purpose
% CALLEE:  mgso4_loader polyvaln

global alphainterpolant

if isempty(alphainterpolant)
    disp('Loading MgSO4 data')
    mgso4_loader
end

alpha = alphainterpolant(P,S,T);