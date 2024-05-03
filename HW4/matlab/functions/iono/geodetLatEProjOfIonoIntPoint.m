function phi_i = geodetLatEProjOfIonoIntPoint(phi_u, psi, A)
phi_i = phi_u + psi*cos(A);
if phi_i > 0.416
    phi_i = 0.416;
elseif phi_i < -0.416
    phi_i = -0.416;
end
end