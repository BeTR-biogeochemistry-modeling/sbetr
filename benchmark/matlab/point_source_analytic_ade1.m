function y=point_source_analytic_ade1(z,L,u,D,time)
%point source analytical solution to equation
%\p_C/\p_t=\p_z(D\p_C\p_z)-u\p_C/p_z
%Tang et al., 2013, GMD., equations (22-23)
%z: depth, [m]
%L: total depth, [m]
%u: advection velocity, [m/s]
%D: diffusivity, [m2/s]
%time: time, [s]
y=0.5.*erfc((z-u.*time)/(2.*sqrt(D.*time)))+0.5.*exp(u.*z./D)...
    .*erfc((z+u.*time)/(2.*sqrt(D.*time)))+(1+u./(2.*D).*(2.*L-z+u.*time))...
    .*exp(u.*L./D).*erfc((2.*L-z+u.*time)./(2.*sqrt(D.*time)))-...
    sqrt(u*u*time/(pi.*D)).*exp(u.*L./D-(2.*L-z+u.*time).^2./(4.*D.*time));




end