% ======================================================
% Run a advection-diffusion equation
% ======================================================
function [t,z,P] = AdvectionDiffusion(dz) % dz is the grid spacing
%
% Parameters:
%
param.D = 3; % Diffusivity (m2/day)
param.u = 1; % Settling velocity (m/day)

param.depth = 100; % meter
param.dz = dz; % Grid spacing (m)
param.z = param.dz/2:param.dz:(param.depth-param.dz/2);
param.nGrid = length(param.z);  % No. of grid cells
%
% Initialization:
%
P0 = exp(-(param.z-param.depth/2).^2/5);
%
% Run model
%
[t, P] = ode45(@PmodelDeriv, 0:200, P0, [], param);
z = param.z;

    % ---------------------------------------------------------
    % Derivative function
    % ---------------------------------------------------------
    function dPdt = PmodelDeriv(t,P,param)
        %
        % Advective fluxes
        %
        ix = 2:(param.nGrid);
        Jadv(ix) = param.u*P(ix-1);
        Jadv(1) = 0;  % No input from the surface
        Jadv(param.nGrid+1) = 0; % Closed bottom
        %
        % Diffusive fluxes:
        %
        Jdiff(ix) = -param.D*(P(ix)-P(ix-1))/param.dz;
        Jdiff(1) = 0; % No flux at the surface...
        Jdiff(param.nGrid+1) = 0;  % ...or the bottom
        %
        % Rate-of-change due to advection and diffusion:
        %
        J = Jadv + Jdiff;
        dPdt = -(J(2:(param.nGrid+1))-J(1:param.nGrid))/param.dz;
        % Make dPdt a column vector:
        dPdt = dPdt';
    end

end
