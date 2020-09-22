function y=trpoint_source_analytic_ade(z,u,D,time)
%point source analytical solution to equation
%\p_C/\p_t=\p_z(D\p_C\p_z)-u\p_C/p_z
%Tang et al., 2013, GMD., equations (22-24)
%z: depth, [m]
%u: advection velocity, [m/s]
%D: diffusivity, [m2/s]
%time: time, [s]
C0=12./23.;
A1=9./23.;
A2=2./23.;
w1=2.*pi/(365.*86400);
w2=2.*pi/86400;
y=C0+A1.*funaj(u,D,w1,z,time)+A2.*funaj(u,D,w2,z,time);




end

function aj=funaj(u,D,w,z,time)
tmp=sqrt(u.*u+sqrt(u.^4+16.*D*D*w*w));
aj=exp(-u./(2.*D)-sqrt(2)./(4.*D).*tmp).*sin(w.*time-sqrt(2.).*w.*z./tmp);
end