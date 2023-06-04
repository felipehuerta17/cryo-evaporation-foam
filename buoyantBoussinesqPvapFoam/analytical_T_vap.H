// Calculate average vapour temperature
// "T_vap" in this code
T_vap = T_air.value() + 1 / (V_vap.value()/A_T.value()) * 
(c_1 / chi_minus * ( Foam::exp(V_vap.value()/A_T.value() * chi_minus) - 1) +
 c_2 / chi_plus * (Foam::exp(V_vap.value()/A_T.value() * chi_plus) - 1) );

// Update cryogen vapour thermal conductivity
k_V = kV_a0 + kV_a1 * T_vap + kV_a2 * pow(T_vap, 2) + 
kV_a3 * pow(T_vap, 3);
k_V_roof = k_V;

Info << "T_v_avg_analytical = " << T_vap << " K " << endl;
// Debug prints
// if(debug == true)
// {
//    Info << " k_V =" << k_V << " W/mK " << endl; 
// }
