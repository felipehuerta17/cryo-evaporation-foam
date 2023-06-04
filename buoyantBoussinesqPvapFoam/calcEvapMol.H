// Wall evaporation
scalar nvap_wall = 0;
double relax = 1;

// Define scalar to track saturation temperature as a function of height
// considering the hydrostatic head
dimensionedScalar T_sat_z = 0;
// prgh = p + rho * g * h
double prgh = 0;

forAll(mesh.C(),i)
{
   // Convert liquid density to kg m ^-3 and then ad hydrostatic head rho * g * (z_int - z)
   prgh = p_RK[0] + (rho_L.value() * MW.value() / 1e3) * 9.81 * (tank_height.value() * LF.value() - mesh.C()[i].z());
   // Calculate height dependent saturation temperature
   T_sat_z = B_ant / (A_ant - Foam::log10(prgh/1e5)) - C_ant;
   // If the temperature of the cell center exceeds the local saturation temperature, evaporate the liquid
   // that will make the cell be exactly at the saturation temperature, and progress to the next time-step

   if(T[i] > T_sat_z.value())
   {
      // Convert sensible heat in latent heat and evaporation rate
       nvap_wall = nvap_wall + relax * (T[i] - T_sat_z.value()) * mesh.V()[i] * (rho_L.value()) * Cp_l.value()/deltaH.value();
      // Set cell temperature to T_sat
      T[i] = T[i] + (T_sat_z.value()-T[i])*relax;
      // Debug prints
      // Info << " T[i] = " << T[i] << endl;
      // Info << "Mesh cell center " << mesh.C()[i].z() << endl;
      // Info << "P_rgh = " << prgh << " Pa " << endl;
      // Info << "T_sat = " << T_sat_z.value() << "K" << endl;
   }
}
// As calculations are made in the whole cylindrical tank, the cylinder
// must be reconstructed from the wedge

nvap_wall *= 360/5;
Info << " Surface evaporation rate  = " << nvap_wall << " mol / s " << endl;
