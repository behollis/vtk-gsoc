clear all, clc, clf ,close all

%This is the root directory where the saved files are stored
path = '../../Save/UCSC_vis/S15';

%load param file
load(sprintf('%s/param.mat', path));

%Load the desire output file
sol = load(sprintf('%s/00050.mat', path));

%Orthonormalize the solution (to plot solution in a consistent rotation of
%the coordinate system). 
[sol.YYt sol.ui sol.vi sol.Pi, sol.rhoi] = ...
        DOorthnorm(app, dx, dy, sol.YYt, sol.ui, uid, sol.vi, vid,...
        sol.Pi, pid, sol.rhoi);
    
    
%Plot the mean density, mean u-velocity, and mean v-velocity
figure(1)
pcolor(XP, YP, sol.rho(NodeP(2:app.Ny+1, 2:app.Nx+1))), shading interp
figure(2)
pcolor(XU, YU, sol.u(Nodeu(2:app.Ny+1, 2:app.Nx))), shading interp
figure(3)
pcolor(XV, YV, sol.v(Nodev(2:app.Ny, 2:app.Nx+1))), shading interp

%Plot one of the modes for density, u-velocity, and v-velocity
%Which mode to plot
modeno = 3;
%Plots
figure(4)
rhoi = sol.rhoi(:, modeno);
ui = sol.ui(:, modeno);
vi = sol.vi(:, modeno);
pcolor(XP, YP, rhoi(NodeP(2:app.Ny+1, 2:app.Nx+1))), shading interp
figure(5)
pcolor(XU, YU, ui(Nodeu(2:app.Ny+1, 2:app.Nx))), shading interp
figure(6)
pcolor(XV, YV, vi(Nodev(2:app.Ny, 2:app.Nx+1))), shading interp
%plot the 1D marginal of the same mode
figure(7)
ksdensity(sol.YYt(:, modeno))

%Now, create a realization, and plot the density, u-velocity and v-velocity
%of that realization
realno = 500;
rho_r = sol.rho + sol.rhoi * YYt(realno,:)';
u_r = sol.u + sol.ui * YYt(realno,:)';
v_r = sol.v + sol.vi * YYt(realno,:)';
%Plots
figure(8)
pcolor(XP, YP, rho_r(NodeP(2:app.Ny+1, 2:app.Nx+1))), shading interp
figure(9)
pcolor(XU, YU, u_r(Nodeu(2:app.Ny+1, 2:app.Nx))), shading interp
figure(10)
pcolor(XV, YV, v_r(Nodev(2:app.Ny, 2:app.Nx+1))), shading interp