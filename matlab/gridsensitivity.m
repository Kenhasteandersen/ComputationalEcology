dz = [0.2 0.4 0.8, 1.6, 3.2, 6.4];

clf
for i = 1:length(dz)
    [t,z,P] = AdvectionDiffusion(dz(i));
    plot(P(end,:), -z, 'o-')
    drawnow
    hold on
end
%%
legend()
xlabel('Concentration (X/m^3)')
ylabel('Depth (m)')
