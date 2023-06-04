// Updates vapour volume considering liquid thermal expansion
auto dTL_dt = (T.weightedAverage(mesh.V()) - TL_av)/runTime.deltaTValue();
// Update TL_av
TL_av = T.weightedAverage(mesh.V());
// Use density polinomials in kg m^-3 and transform in mol m ^-3
auto mol_over_kg = 1e3 /  MW.value();     
// drhoLdt: read density polynomials and calculate the derivative
scalar drho_L_dTL = rhoL_a1.value() + 2 * rhoL_a2.value() + 3* rhoL_a3.value();
auto drho_L_dt = mol_over_kg * (dTL_dt * drho_L_dTL) ;
dimensionedScalar rho_L_exp = mol_over_kg * ( scalar(rhoL_a0.value()) + rhoL_a1 * TL_av.value()  );
// Note that the evaporation rate and all units are in dmdt
// Calculate vapour density derivative
auto drho_v_dt = (rho_gas - rho_old)/runTime.deltaTValue();
auto V_V_exp = (dmdt.value() * (1 - rho_old.value()/rho_L_exp) -
                rho_old.value()/rho_L_exp * (drho_L_dt.value()) * V_tank ) /
            (drho_v_dt.value() - rho_old.value()/rho_L.value() * drho_L_dt.value() );
rho_L = rho_L_exp;

if (debug == true)
{
    // Liquid thermal expansion debug prints
    Info << "TL_av: " << TL_av << endl;
    Info << "dTL_dt: " << dTL_dt << endl;
    Info << "rho_L_exp" << rho_L_exp.value() << endl;
    Info << "drho_L_dt = " << drho_L_dt.value() << endl;
    Info << "drho_v_dt = " << drho_v_dt.value() << " m^3 s^-1 " << endl;
    Info << "V_V_exp = " << V_V_exp.value() << endl;
    Info << "V_tank = " << V_tank.value() << endl;
    Info << "V_vap = " << V_vap.value() << endl;
}
