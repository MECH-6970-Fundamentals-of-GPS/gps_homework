function lambda_i = geodetLonEProjOfIonoIntPoint(lambda_u, phi_i, psi, A)
lambda_i = lambda_u + (psi*sin(A))/cos(phi_i);
end