    ode_temp.Q_V = Q_V;
    ode_temp.Q_VL = Q_VL;
    ode_temp.T_sat = T_sat;
    ode_temp.deltaH = deltaH;
    ode_temp.Cv = Cv;
    ode_temp.Cp_V = Cp_V;
    ode_temp.k_V = k_V;
    ode_temp.rho_L = rho_L;
    ode_temp.V_vap = V_vap;
    ode_temp.qL_int_flux = qL_int_flux;
    ode_temp.n_evaporated = n_evaporated;
    ode_temp.rho_gas = rho_gas;
    ode_temp.T_vap = T_vap;
    ode_temp.p_new = p_new;
    ode_temp.r_q = r_q;
    ode_temp.h_VL = h_VL;
    ode_temp.A_T = A_T;
    // Conjugate heat transfer
    ode_temp.U_dry = U_dry;
    ode_temp.eta_evap = eta_evap;
    ode_temp.d_o = d_o;
    ode_temp.d_i = d_i;
    ode_temp.dmdt = dmdt;
    ode_temp.T_air = T_air.value();
    ode_temp.nz_v = nz_v;
    ode_temp.f_c = f_c;
    ode_temp.debug = false;

    // Initialise vapor temperature at the saturation temperature
    // of the stored cryogen at the initial conditions

    // Initialise dimensionless mesh
    for(label i=0; i < nz_v; i++)
    {
        vapour_mesh[i] = 1.0/(nz_v-1) * i;
    }

    for(label i=0; i <= temp_v.size(); i++)
    {
        // temp_v[i] = T_vap.value();
        temp_v[i] = temp_v_a0.value() + temp_v_a1.value() * vapour_mesh[i] +
        temp_v_a2.value() * pow(vapour_mesh[i], 2) +
        temp_v_a3.value() * pow(vapour_mesh[i], 3);
    }
    ode_temp.dz_v = dz_v;