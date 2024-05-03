function phi_m = geomagLatEProjOfIonoIntPoint(phi_i, lambda_i)
phi_m = phi_i + 0.064*cos(lambda_i-1);
end