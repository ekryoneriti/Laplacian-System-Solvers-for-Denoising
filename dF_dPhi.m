function der_phi = dF_dPhi(v_p, u, phi, Sigma_p, Sigma_u)
    der_phi = ((v_p - phi)/Sigma_p) + ((u - g(phi))/Sigma_u)*g_prime(phi);
end