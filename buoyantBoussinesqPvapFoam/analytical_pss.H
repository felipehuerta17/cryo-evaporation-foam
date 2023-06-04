// Analytical solutions for pseudo-steady state heat transfer
// In the vapour phase. Calculates Q_VI and the average
// vapour temperature. Based in Supp. Material from Huerta and Vesovic (2019)
// IJHMT article

// Advective velocity based on evaporation rate

// No advective velocity!
H = 0;
// do/d_i^2 ~ 1/d_i
S = 4 * U_dry.value() * d_o.value() / pow(d_i.value(),2);
// T_air is assumed constant and read from the transportProperties
// dictionary, see readTransportProperties.H
E = S * T_air.value();
// Constant for the calculation of the analytical solutions
chi_plus = (H + Foam::sqrt( pow(H,2) + 4 * k_V.value() * S)) / (2 * k_V.value()) ;
chi_minus = (H - Foam::sqrt( pow(H,2) + 4 * k_V.value() * S)) / (2 * k_V.value());

b_plus = Foam::exp(V_vap.value()/A_T.value() * chi_plus);
b_minus = Foam::exp(V_vap.value()/A_T.value() * chi_minus);

a_plus = chi_plus * b_plus;
a_minus = chi_minus * b_minus;

// Boundary conditions
if(name_BC == "Neumann")
{
    c_1 = (a_plus * (T_sat.value() - T_air.value()) - gradT_roof.value())  / (a_plus - a_minus);
    c_2 = (a_minus * (T_air.value() - T_sat.value()) + gradT_roof.value())  / (a_plus - a_minus);
}
else
{
    gamma_roof = U_roof.value() / k_V_roof.value();
    c_1 = (T_sat.value() - T_air.value()) * (a_plus + gamma_roof * b_plus) /
        (a_plus + gamma_roof*b_plus - (a_minus + gamma_roof*b_minus));
    c_2 = (T_air.value() - T_sat.value()) * (a_minus + gamma_roof * b_minus) /
        (a_plus + gamma_roof*b_plus - (a_minus + gamma_roof * b_minus) );
}
