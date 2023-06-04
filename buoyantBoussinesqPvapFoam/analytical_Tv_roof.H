// Calculate roof temperature
Tv_roof = c_1 * Foam::exp(V_vap.value()/A_T.value() * chi_minus) +
    c_2 * Foam::exp(V_vap.value()/A_T.value() * chi_plus) + T_air.value();
