singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// string vapourHeatTransferModel("vapourHeatTransferModel");
dimensionedScalar vapourHeatTransferModel("vapourHeatTransferModel", dimless, laminarTransport);

// Vapour temperature polynomial
dimensionedScalar temp_v_a0("temp_v_a0", dimless, laminarTransport);
dimensionedScalar temp_v_a1("temp_v_a1", dimless, laminarTransport);
dimensionedScalar temp_v_a2("temp_v_a2", dimless, laminarTransport);
dimensionedScalar temp_v_a3("temp_v_a3", dimless, laminarTransport);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

dimensionedScalar T_air("T_air", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);
    
// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Specific heat capacity at constant volume
// Note that the derivation was performed
// dimensionedScalar Cv("Cv", dimensionSet(1, 2, -2, -1, -1, 0, 0), laminarTransport);
dimensionedScalar Cp_l(64.614);

dimensionedScalar Cv("Cv", dimless, laminarTransport);

// Vapour heat capacity at constant pressure
dimensionedScalar Cp_V("Cp_V", dimless, laminarTransport);
// dimensionedTensor Cp_V_coeffs("Cp_V_coeffs", dimless, laminarTransport);
// Latent heat of vaporization
dimensionedScalar deltaH("deltaH", dimless, laminarTransport);

// Antonie coefficients in log_10 basis
// T_sat = B / (A - log10(P)) - C
dimensionedScalar A_ant("A", dimless, laminarTransport);
dimensionedScalar B_ant("B", dimless, laminarTransport);
dimensionedScalar C_ant("C", dimless, laminarTransport);

// Defining saturation temperature
dimensionedScalar T_sat = B_ant / (A_ant - Foam::log10(pvap[1]/1e5)) - C_ant;
// Molar weight
dimensionedScalar MW("MW", dimensionSet(1, 0, 0, 0, -1, 0, 0), laminarTransport);

// Liquid thermophysical properties

// Molar liquid density
dimensionedScalar rho_L("rho_L", dimless, laminarTransport);
// Transform to molar basis
// Liquid viscosity
rho_L = rho_L / MW.value() * 1000;
// const volScalarField& nu = laminarTransport.nu();

// Liquid thermal conductivity W m^-2 K^-1
dimensionedScalar k_L("k_L", dimless, laminarTransport);
// Vapour thermal conductivity W m^-2 K^-1
dimensionedScalar k_V("k_V", dimless, laminarTransport);
// Vapour to liquid to total vapour heat ingress ratio
dimensionedScalar r_q("r_q", dimless, laminarTransport);
// Convection coefficient
dimensionedScalar h_VL("h_VL", dimless, laminarTransport);

dimensionedScalar T_vap(T_sat);

if (r_q.value() < 1){
    T_vap = Tvap[1];
} 

// Calculate average vapour temperature
if (vapourHeatTransferModel.value() == 5)
{
    T_vap = temp_v_a0 + temp_v_a1/2 + temp_v_a2/3 + temp_v_a3/4;
}

dimensionedScalar Tv_roof(T_vap);

// Vapour temperature

// Error
dimensionedScalar pbuild_err = 1;

// Required for iterations
dimensionedScalar p_old, p_new;
dimensionedScalar T_old, T_new;

dimensionedScalar Q_VL = 0; // Estimation of the vapour to liquid heat transfer
// Initialice evaporation rate, >0 evap, <0 cond.
dimensionedScalar dmdt = 0;

// Tank transport properties

dimensionedScalar U_dry("U_dry", dimless, laminarTransport);
// Overall heat transfer coefficient for the liquid phase
dimensionedScalar U_L("U_L", dimless, laminarTransport);
//dimensionedScalar U_dry = 0.039;


dimensionedScalar d_i("d_i", dimless, laminarTransport); // Internal diameter [m]
dimensionedScalar d_o("d_o", dimless, laminarTransport); // External diameter [m]

dimensionedScalar tank_height("tank_height", dimless, laminarTransport); // Tank height
dimensionedScalar LF("LF", dimless, laminarTransport); // Initial liquid filling
dimensionedScalar A_T = constant::mathematical::pi *pow(d_i,2)/4;
dimensionedScalar V_tank = A_T*tank_height;
// Vapour volume
dimensionedScalar V_vap = V_tank*(1-LF); //0.00174163;


// Initialize density with ideal gas equation
dimensionedScalar rho_gas = pvap[1]/(8.314*T_vap);
// dimensionedScalar V_vap = V_tank * (1-LF);

volScalarField V_vapfield
(
    IOobject
    (
        "V_vap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    (V_tank * (1-LF))
);

volScalarField V_tankfield
(
    IOobject
    (
        "V_tank",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    V_tank
);

// Defining number of moles evaporated
dimensionedScalar n_evaporated = pvap[1] * V_vap / (8.314 * T_sat);
// dimensionedScalar n_evaporated = n_evap[1];
// Number of moles evaporated at each time-step

dimensionedScalar nvap_dt = 0; 
dimensionedScalar A_dry = V_vap/A_T * constant::mathematical::pi * d_i;

dimensionedScalar A_wet = (V_tank-V_vap)/A_T * constant::mathematical::pi * d_i;

dimensionedScalar qL_int_flux = 0;

dimensionedScalar Q_V = U_dry * A_dry*(T_air.value()-T_vap);// Vapour heat ingress
// Heat transfer from the wall that produces evaporation
// through wall boiling / W
dimensionedScalar Q_L_wb(0);

// eta evap: heat transfer partitioning q = q_e/q_tot (add quenching in future models)
dimensionedScalar eta_evap("eta_evap", dimless, laminarTransport); // Internal diameter [m]

// Non-isobaric evaporation model
// Initialise constants as 0
scalar H(0);
scalar S(0);
scalar E(0);
scalar chi_plus(0);
scalar chi_minus(0);
scalar a_plus(0);
scalar a_minus(0);
scalar b_plus(0);
scalar b_minus(0);
scalar c_1(0);
scalar c_2(0);
scalar gamma_roof(0); // Convection to conduction heat transfer ratio
// Setup
word name_BC("Robin");

// Temperature gradient at the tank roof if Neumann
// boundary conditions are selected
dimensionedScalar gradT_roof("gradT_roof", dimless, laminarTransport);
// Roof overall heat transfer coefficient
dimensionedScalar U_roof("U_roof", dimless, laminarTransport);
// Cryogen vapour thermal conductivity at the tank roof
// Suggestion: use k_V = k_V(Tv_av) or less accurately k_V = k_V(Tsat)
// If the profile is not known in advance
dimensionedScalar k_V_roof("k_V_roof", dimless, laminarTransport);
// Heat transfer from the roof
dimensionedScalar Q_roof =  U_roof * A_T * (T_air.value()-T_vap);
// q_roof for Neumann BC
dimensionedScalar q_roof("q_roof", dimless, laminarTransport);

// Density polynomial, the input should be in kg / m^3 but
// in internal calculation mol / m^3 are used.
dimensionedScalar rhoL_a0("rhoL_a0", dimless, laminarTransport);
dimensionedScalar rhoL_a1("rhoL_a1", dimless, laminarTransport);
dimensionedScalar rhoL_a2("rhoL_a2", dimless, laminarTransport);
dimensionedScalar rhoL_a3("rhoL_a3", dimless, laminarTransport);

// Vapour thermal conductivity polynomial (W/mK)
dimensionedScalar kV_a0("kV_a0", dimless, laminarTransport);
dimensionedScalar kV_a1("kV_a1", dimless, laminarTransport);
dimensionedScalar kV_a2("kV_a2", dimless, laminarTransport);
dimensionedScalar kV_a3("kV_a3", dimless, laminarTransport);




// Number of nodes in the vertical subdomain
dimensionedScalar nz_V("nz_V", dimless, laminarTransport);

IOdictionary transportProperties
(
        IOobject
        (
            "transportProperties",  //name of what to read
            runTime.constant(), //where to read it from 
            mesh,               //what objectRegistry it is a part of (mesh)
            IOobject::MUST_READ, //Does it read? (yes)
            IOobject::NO_WRITE   //Does it write? (no)
        )
);

// Non-isobaric evaporation
word roof_bc = word(transportProperties.lookup("roof_BC"));

// Coefficients for bulk phase change based in non-eq. thermodynamics
// Condensation coefficient
scalar f_c = readScalar(transportProperties.lookup("f_c"));
// Evaporation coefficient
scalar f_e = readScalar(transportProperties.lookup("f_e"));