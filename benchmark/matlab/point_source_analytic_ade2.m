function y=point_source_analytic_ade2(z,L,u,D,time)
%point source analytical solution to equation
%\p_C/\p_t=\p_z(D\p_C\p_z)-u\p_C/p_z
%Tang et al., 2013, GMD., equations (22-23)
%z: depth, [m]
%L: total depth, [m]
%u: advection velocity, [m/s]
%D: diffusivity, [m2/s]
%time: time, [s]
y=0.5.*eerfcs((z-u.*time)/(2.*sqrt(D.*time)),zeros(size(z)))+0.5...
    .*eerfcs((z+u.*time)/(2.*sqrt(D.*time)),u.*z./D)+...
    (1+u./(2.*D).*(2.*L-z+u.*time))...
    .*eerfcs((2.*L-z+u.*time)./(2.*sqrt(D.*time)),u.*L./D.*ones(size(z)))-...
    sqrt(u*u*time/(pi.*D)).*exp(u.*L./D-(2.*L-z+u.*time).^2./(4.*D.*time));




end