function [YYt ui vi Pi rhoi] = DOorthnorm(app, dx, dy, YYt, ui, uid, vi, vid, Pi, pid, rhoi)
% function [YYt u v P rho a] = DOorthnorm(YYt, ui, uid, vi, vid, Pi, pid, rhoi)
% This function orthonormalizes the DO modes while preserving the original
% realizations (upto a correction for the stochastic energy) and the
% stochastic energy.
%
% INPUTS:
%   dx:         Size of x discretization
%   dy:         Size of y discretization
%   YYt(MC,S):  Stochastic coefficients
%   ui(:,S):    Mode shapes for u-velocity
%   uid:        Vector of ids for interior u-velocities
%   vi(:,S):    Mode shapes for v-velocity
%   vid:        Vector of ids for interior v-velocities
%   Pi(:,S):    Mode shapes for Pressure
%   Pid:        Vector of ids for interior Pressures
%   rhoi(Optional): Mode shapes for density
%
% OUTPUTS:
%   YYt(MC,S):  Stochastic coefficients
%   ui(:,S):    Mode shapes for u-velocity
%   vi(:,S):    Mode shapes for v-velocity
%   Pi(:,S):    Mode shapes for Pressure
%   rhoi(Optional): Mode shapes for density
%
% Written by: Matt Ueckermann

if ~isfield(app,'IP'),app.IP=[1,1,1];end
if ~isfield(app,'Sb'), app.Sb = 0;end

%% Old orthonormalization procedure
if isfield(app,'GSorth')
%     display('Using old orthonormalization');
    for i=1:app.S(1)
        for j=1:app.S(1)
            if i~=j
                if nargin == 11
                    coef =app.IP(1) * InProd(ui(:,i), ui(:,j), dx, dy, uid)...
                        + app.IP(2) * InProd(vi(:,i), vi(:,j), dx ,dy, vid)...
                        + app.IP(3) * InProd(rhoi(:,i), rhoi(:,j), dx, dy, pid);
                    rhoi(:,i) = rhoi(:,i) - coef * rhoi(:,j);
                else
                    coef = InProd(ui(:,i), ui(:,j), dx, dy, uid)...
                        + InProd(vi(:,i), vi(:,j), dx, dy, vid);
                end
                ui(:,i) = ui(:,i) - coef * ui(:,j);
                vi(:,i) = vi(:,i) - coef * vi(:,j);
            end
        end
        if nargin == 11
            nrm = sqrt(app.IP(1)*InProd(ui(:,i), ui(:,i), dx, dy, uid)...
                + app.IP(2)*InProd(vi(:,i), vi(:,i), dx ,dy, vid)...
                + app.IP(3)*InProd(rhoi(:,i), rhoi(:,i), dx, dy, pid));
            rhoi(:,i) = rhoi(:,i) / nrm;
        else
            nrm = sqrt(InProd(ui(:,i), ui(:,i), dx, dy, uid)...
                + InProd(vi(:,i), vi(:,i), dx, dy, vid));
        end        
        ui(:,i) = ui(:,i) / nrm;
        vi(:,i) = vi(:,i) / nrm;
    end
else
    %% First rotate the stochastic coefficients so that we have
    %% uncorrelated samples    
    Sintid = app.Sb+1:app.S;
    CYY=cov(YYt(:, Sintid));
    [VC, tmp] = eig(CYY) ; %Diagonal Covariance, and Diagonal Vectors
    ui(:, Sintid) = fliplr(ui(:, Sintid) * VC);
    vi(:, Sintid) = fliplr(vi(:, Sintid) * VC);
    Pi(:, Sintid) = fliplr(Pi(:, Sintid) * VC);
    if nargin == 11
        rhoi(:,Sintid) = fliplr(rhoi(:, Sintid) * VC);
    else
        rhoi = 0;
    end
    YYt(:,Sintid) = fliplr(YYt(:,Sintid) * VC);
    
    %% Next orthonormalize the modes
    if nargin == 11
        M = (app.IP(1)*ui(uid, :)' * ui(uid, :) + app.IP(2)*vi(vid, :)' * vi(vid, :) ...
            + app.IP(3) * rhoi(pid, :)' * rhoi(pid, :)) * dx * dy;
    else
        M = (ui(uid, :)' * ui(uid, :) + vi(vid, :)' * vi(vid, :)) * dx * dy;
    end
    [VC DC] = eig(M);
    YYt = fliplr(YYt * VC * sqrt(DC));
%     YYt = fliplr(YYt * VC);
    DC = diag (1 ./ diag(sqrt(DC)));
    ui = fliplr(ui * VC * DC);
    vi = fliplr(vi * VC * DC);
    if nargin == 11
        rhoi = fliplr(rhoi * VC * DC);
    end
    Pi = fliplr(Pi * VC * DC);
    Nbcs = length(ui) - length(uid);
    %Now rotate the modes such that the boundary is orthogonal once again
    %Now do the boundary-inner product to separate edge modes from non-boundary
    %modes
    if nargin == 11
        [VC,DC] = eig(ui(1:Nbcs,:)'*ui(1:Nbcs,:) + vi(1:Nbcs,:)'*vi(1:Nbcs,:) + ...
            rhoi(1:Nbcs,:)'*rhoi(1:Nbcs,:));
        rhoi = fliplr(rhoi*(VC));
    else
         [VC,DC] = eig(ui(1:Nbcs,:)'*ui(1:Nbcs,:) + vi(1:Nbcs,:)'*vi(1:Nbcs,:));
    end
    ui = fliplr(ui*(VC));
    vi = fliplr(vi*(VC));
    Pi = fliplr(Pi*(VC));
    YYt = fliplr(YYt*(VC));

    %% Rotate the stochastic coefficients back to the uncorrelated case
    CYY=cov(YYt(:, Sintid));
    [VC,DC] = eig(CYY) ; %Diagonal Covariance, and Diagonal Vectors
    ui(:, Sintid) = fliplr(ui(:, Sintid) * VC);
    vi(:, Sintid) = fliplr(vi(:, Sintid) * VC);
    Pi(:, Sintid) = fliplr(Pi(:, Sintid) * VC);
    if nargin == 11
        rhoi(:, Sintid) = fliplr(rhoi(:, Sintid) * VC);
    end
    %Also correct for, or make sure that, stochastic energy is preseved.
    YYt(:, Sintid) = fliplr(YYt(:, Sintid) * VC * ...
        sqrt(sum(tmp(:)) / sum(DC(:))));
%     YYt=fliplr(YYt * VC);
end
