/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantBoussinesqPimpleFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "trapezoidal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class pvODE
:
    public ODESystem
{

public:

    // Empty constructor
    pvODE()
    {}

    label nEqns() const
    {
        return 1;
    }

    // Define public variables that can be accessed
    // In the ODE system
    dimensionedScalar r_q;
    dimensionedScalar Q_V;
    dimensionedScalar Q_VL;
    dimensionedScalar T_sat;
    dimensionedScalar deltaH;
    dimensionedScalar Cv;
    dimensionedScalar Cp_V;
    dimensionedScalar rho_L;
    dimensionedScalar V_vap;
    dimensionedScalar T_vap;
    dimensionedScalar rho_gas;
    dimensionedScalar n_evaporated;
    dimensionedScalar qL_int_flux;
    dimensionedScalar p_new;
    dimensionedScalar vapourHeatTransferModel;
    dimensionedScalar A_T;
    dimensionedScalar h_VL;
    dimensionedScalar dmdt_wall;
    dimensionedScalar Q_VI;
    dimensionedScalar Q_L_wb;
    dimensionedScalar eta_evap;
    //
    dimensionedScalar dmdt;
    bool debug;


    void derivatives
    (
        const scalar t,
        const scalarField& p,
        scalarField& dpdt
    ) const
    {
        if (r_q.value() < 1) { 
            // The vapour is allowed to be superheated, and the
            // vapour to liquid heat ingress is controlled by
            // r_q = Q_VL/Q_V_in. or a convection coefficient
            dimensionedScalar Q_VI(0);
            // Define vapour to liquid heat ingress
            if (vapourHeatTransferModel.value() == 2)
                {
                    // Fixed r_q
                    Q_VI = r_q*Q_V;
                }
            else if (vapourHeatTransferModel.value() == 3)
            {
                // convectionCoefficient model
                Q_VI = h_VL * A_T * (T_vap - T_sat);
                Info << "Q_VI convection = " << Q_VI.value() << endl;
            }
  
            dimensionedScalar LHS =  (Q_V - Q_VI) + (Q_VI + qL_int_flux + dmdt_wall * deltaH + Q_L_wb)*Cp_V*T_sat/(deltaH);

            // Info << " LHS = " << LHS.value() << endl;
            // Info << " Q_VI = " << Q_VI.value() << endl;

            dimensionedScalar dpdtdim = 8.314 / (Cp_V*V_vap) * LHS - p[0] * (Q_VI + qL_int_flux + dmdt_wall*deltaH + Q_L_wb) / (deltaH * rho_L);

            dpdt[0] = dpdtdim.value();
        } else {
            // The vapour is at equilibrium with the liquid
            dimensionedScalar F = deltaH/ V_vap * pow((Cv*T_sat + (deltaH/(8.314*T_sat)-1) * rho_L / (rho_L - rho_gas) *
            (deltaH - p[0] * (1/rho_gas - 1/rho_L) )), -1);
            dimensionedScalar dpdtdim = (F * (Q_V + qL_int_flux + Q_L_wb + dmdt_wall*deltaH));

            dpdt[0] = dpdtdim.value();
            
            // Optional debug prints
            if(debug == true)
            {
                Info << "Q_V_ode = " << Q_V.value() << "W" << endl;
                Info << "qL_int_flux = " << qL_int_flux.value() << "W" << endl;
                Info << "Q_L_wb = " << Q_L_wb.value() << endl;
            }
        }
    }

    // No jacobian at the beginning
    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {

    }
};

class tempv_ODE
:
    public ODESystem
{

public:

    // Empty constructor
    tempv_ODE()
    {}

    label nEqns() const
    {
        return nz_v;
    }

    // Define public variables that can be accessed
    // In the ODE system
    // Number of nodes
    label nz_v;
    // Heat flux ratio
    dimensionedScalar r_q;
    dimensionedScalar Q_V;
    dimensionedScalar Q_VL;
    dimensionedScalar T_sat;
    dimensionedScalar deltaH;
    dimensionedScalar Cv;
    dimensionedScalar Cp_V;
    dimensionedScalar rho_L;
    dimensionedScalar V_vap;
    dimensionedScalar T_vap;
    dimensionedScalar rho_gas;
    dimensionedScalar n_evaporated;
    dimensionedScalar qL_int_flux;
    dimensionedScalar p_new;
    dimensionedScalar vapourHeatTransferModel;
    dimensionedScalar A_T;
    dimensionedScalar h_VL;
    dimensionedScalar dmdt_wall;
    dimensionedScalar Q_VI;
    dimensionedScalar Q_L_wb;
    dimensionedScalar eta_evap;
    // Grid spacing in the vertical direction, uniform mesh
    dimensionedScalar dz_v;
    // Vapour phase thermal conductivity
    dimensionedScalar k_V;
    // Vapour phase heat transfer coefficient
    dimensionedScalar U_dry;
    // Tank internal and external diameters
    dimensionedScalar d_o;
    dimensionedScalar d_i;
    // Temperature of the surroundings
    dimensionedScalar T_air;
    dimensionedScalar dmdt;
    // Bulk evaporation and condesation coefficients
    scalar f_e;
    scalar f_c;
    bool debug;

    void derivatives
    (
        const scalar t,
        const scalarField& temp_v,
        scalarField& dTv_dt
    ) const
    {
    // Huerta and Vesovic (2019) vapour phase 1-D model extended with evaporation
    // and condensation switches
    // Vapour thermal conductivity
    // Initial vertical spacing
    // Info << "deltaz_v = " << dz_v << "m" << endl; 
    // Finite difference coefficients
    double alpha_v = k_V.value() / (rho_gas.value() * Cp_V.value());
    // Advective velocity
    double v_z = dmdt.value() / (A_T.value() * rho_gas.value());
    // Info << "Vapour velocity "  << v_z << "m/s" << endl;
    // Info << "alpha_V  "  << alpha_v << "m^2/s" << endl;
    // Just for completeness, BC declared outside
    dTv_dt[0] = 0;
    
    // Wall heating, dependent on vapour temperature 
    dimensionedScalar S_wall;
    dimensionedScalar dtv_dt;
    // Bulk condensation
    dimensionedScalar S_c(0);

    // Info << "Vapour temperature array " <<  temp_v << endl;

    for(label i=1; i < temp_v.size() - 1; i++) 
    {
        // Bulk condensation source term
        if(temp_v[i] < temp_v[0]){
            S_c = f_c * (temp_v[0]- temp_v[1]);
        } else {
            S_c = 0;
        }

        S_wall = 4 * U_dry.value() * (1-eta_evap.value()) * d_o.value() / pow(d_i.value(),2) *
        (T_air.value()-temp_v[i]);
        // Heat diffusion
        dtv_dt = alpha_v * (temp_v[i+1] - 2*temp_v[i] + temp_v[i-1]) / pow(dz_v,2)
        // Advection
        // - v_z * (temp_v[i+1] - temp_v[i-1]) / (2*dz_v) +
        + 
        // Source term
        alpha_v/k_V.value() * S_wall - S_c;
                // Info << "dtv_dt" << dtv_dt << "K/s" << endl;
        dTv_dt[i] = dtv_dt.value();
        
        if(debug)
            {
            Info << "Diffusion term = " << alpha_v * (temp_v[i+1] - 2*temp_v[i] + temp_v[i-1]) / pow(dz_v,2) 
            << "K/s" << endl; 
            Info << "Advection term = " << v_z * (temp_v[i+1] - temp_v[i-1]) / (2*dz_v) <<
            " K/s " << endl;
            Info << "S_wall term = " << alpha_v/k_V.value() * S_wall << "K/s" << endl;
            } 
    }

    // Extrapolate temperature derivative, independent of Robin or Neumann BCs
    dTv_dt[temp_v.size()-1] = 1/3.0 * ( 4*dTv_dt[temp_v.size()-2] - dTv_dt[temp_v.size()-3]);

    }

    // No jacobian at the beginning
    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {

    }
};

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    // Liquid thermal conductivity
    Info << "k_L input = " << k_L << endl;

    Info << "\nReading non-standard transport properties" << endl;
    Info << "Cv = " << Cv.value() << endl;
    Info << "DH_LV = " << deltaH.value() << endl;
    Info << "MW = " << MW.value() << endl;
    Info<< "\n Reading tank heat transfer properties" << endl;
    Info << "U_dry = " << U_dry.value() << endl;
    Info << "A_T =" << A_T.value() << endl;
    Info << "A_dry = " << A_dry.value() << endl;
    Info << "Creating ODE system for pressure evolution integration" << endl;
    Info << "V_tank = " << V_tank.value() << endl;
    // Vapour to liquid to total vapour heat ingress ratio
    Info << "r_q = " << r_q.value() << endl;
    Info << "rho_L = " << rho_L.value() << endl;
    
    Info << " Vapour heat transfer model = " << vapourHeatTransferModel << endl;
    Info << "Roof boundary condition: " << roof_bc << endl;
    // Non-equilibrium model prints
    Info << "f_e = " << f_e << endl;
    Info << "f_c = " << f_c << endl;


    // Pressure ODE
    pvODE ode;
    // Temperature ODE system
    tempv_ODE ode_temp;

    // Create ODEsystem
    ode.Q_V = Q_V;
    ode.Q_VL = Q_VL;
    ode.T_sat = T_sat;
    ode.deltaH = deltaH;
    ode.Cv = Cv;
    ode.Cp_V = Cp_V;
    ode.rho_L = rho_L;
    ode.V_vap = V_vap;
    ode.qL_int_flux = qL_int_flux;
    ode.n_evaporated = n_evaporated;
    ode.rho_gas = rho_gas;
    ode.T_vap = T_vap;
    ode.p_new = p_new;
    ode.r_q = r_q;
    ode.h_VL = h_VL;
    ode.A_T = A_T;
    ode.vapourHeatTransferModel = vapourHeatTransferModel;

    dictionary dict;
    dict.add("solver", "RKF45");
    // Create ODE solver for pressure
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode,dict);
    Info << "Initialisiting ODE solver for temperature" << endl;
    // Read vapour properties
    const label nz_v = nz_V.value();
    // Initialise mesh
    scalar dz_v = (V_vap.value()/A_T.value()) / (nz_v);
    // Vapour temperature 1-D array initialisation
    scalarField temp_v(nz_v);
    // Create dimensionless vapour mesh with the initial vapour volume
    scalarField vapour_mesh(nz_v);
    #include "odetemp_init.H"
    Info << "Create temperature object" << endl;
    // Create ODE solver for temperature
    autoPtr<ODESolver> odeSolver_temp = ODESolver::New(ode_temp,dict);
    // Vapour temperature derivative
    scalarField dTv_dt(ode_temp.nEqns());

    // Pressure and pressure derivative
    scalarField p_RK(ode.nEqns()); // Initialize Runge Kutta integration
    scalarField dpdt(ode.nEqns());
    // Initial pressure

    p_RK[0] = pvap[1];
    // Initialize temperature
    if (r_q.value() == 1)
    {
        T_vap = p_RK[0] / (8.314 * rho_gas);
    }

    // If the analytical solutions are being used, initialise the temperature
    // to the pseudo-steady state value for vapour temperature
    if(vapourHeatTransferModel.value() == 4)
    {
        T_old = T_vap;
        scalar err = 1;
        while (err > 1e-4) {
            #include "analytical_pss.H"
            #include "analytical_T_vap.H"
            err = T_vap.value() - T_old.value();
            T_old = T_vap;
        }
    }
    rho_gas = pvap[1] / (8.314 * T_vap);
    // Update n_evaporated to provide consistent
    // initial Pressure
    n_evaporated = V_vap * rho_gas;
    Info << "p0 = " << p_RK[0] << "Pa" << endl;
    Info << "T_vap_0 " << T_vap.value() << endl;
    // Initialize liquid temperature
    auto TL_av = T.weightedAverage(mesh.V());

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        // Access interphase temperature values

        // Find Patch ID for the interphase
        label interface = mesh.boundaryMesh().findPatchID("interphase");

        // Compute T gradients at the interface
        surfaceScalarField SnGrad
        (
            fvc::snGrad(T)
        );
        scalarField gradTL_int = SnGrad.boundaryField()[interface];

        // Get the areas of each face in the interphase boundary
        const scalarField magInterfaceAreas (mesh.magSf().boundaryField()[interface]);
        
        // Add the heat flux along all the sides of the interface

        dimensionedScalar qL_int_flux = 0;
        forAll(magInterfaceAreas, i)
        {
            qL_int_flux += -magInterfaceAreas[i] * gradTL_int[i];
        }
        
        // Recover the whole disk from the 5 degree wedge
        qL_int_flux *= 360/5 * k_L.value() ;
        // Liquid to interface heat transfer rate
        Info << "Q_LI = " << qL_int_flux.value() << " W "<< endl;

        T_old = T_sat;
        //T_sat = 264.651 / (3.7362 - Foam::log10(p_RK[0]/1e5) ) + 6.788;
        T_sat = B_ant / (A_ant - Foam::log10(p_RK[0]/1e5)) - C_ant;
        // Construct vapour temperature polynomials
        #include "analytical_pss.H"

        // Calculate roof temperature
        if (vapourHeatTransferModel.value() == 4)
        {
            // Calculate roof temperature
            #include "analytical_Tv_roof.H"
        }
        else if (vapourHeatTransferModel.value() == 5)
        {
            // Set roof temperature as the vapour temperature of the last node
            Tv_roof = temp_v[nz_v-1];
        }
        
        #include "calcEvapMol.H" // Calculate evaporated moles
        T_new = T_sat;

        // P-V-T isochoric loop
        p_old = p_RK[0];
        // Info << "p_old = " << p_old.value() << " Pa " << endl;
        // dimensionedScalar rho_gas = n_evaporated / V_vap;

        A_dry = V_vap/A_T * constant::mathematical::pi * d_o;
        A_wet = (V_tank-V_vap)/A_T * constant::mathematical::pi * d_o;
        Info << "V_vap = " << V_vap.value() << " m3 " << endl;
        Info << "A_dry = " << A_dry.value() << " m2 " << endl;
        Info << "A_wet = " << A_wet.value() << " m2 " << endl;

        // Vapor heat ingress through the roof and bottom
        Q_V = U_dry * A_dry * (1-eta_evap.value()) * (T_air.value()-T_vap.value()) + 
            U_roof * A_T * (1-eta_evap.value()) * (T_air.value() - Tv_roof.value());

        Info << "Tv_roof = " << Tv_roof.value() << " K " << endl;
        Info << "Q_roof = " << U_roof * A_T * (T_air.value() - Tv_roof.value()) << "W" << endl;

        Q_VL = - qL_int_flux;

        // Calculate wall boiling assuming that a eta_evap fraction
        // of the vapour heat ingress is transferred directly
        // to the interface as latent heat.
        Q_L_wb = eta_evap * U_dry.value() * A_dry.value() * (T_air.value() - T_vap);

        Info << "Q_L_wb = " << Q_L_wb.value() << endl;

        // Store old density and volume
        dimensionedScalar rho_old = rho_gas;
        dimensionedScalar V_old = V_vap;
        // Update vapour volume
        bool debug = false;

        // --- Solve ODE system ---
        // Update ODE coefficients
        ode.Q_V = Q_V;
        ode.T_sat = T_sat;
        ode.V_vap = V_vap;
        ode.qL_int_flux = qL_int_flux;
        ode.Q_L_wb = Q_L_wb;
        ode.n_evaporated = n_evaporated;
        ode.rho_gas = rho_gas;
        ode.T_vap = T_vap;
        ode.dmdt_wall = nvap_wall/runTime.deltaTValue();
        ode.debug = debug;

        // Estimate a delta t estimate
        scalar dtEst = 0.1 * runTime.deltaTValue();
        // Calculate saturation pressure form average density and vapour temperature
        if (vapourHeatTransferModel.value() == 4 ||vapourHeatTransferModel.value() == 5)
            {
                p_RK[0] = rho_gas.value() * 8.314 * T_vap.value();
            }
        else
        // Solve ODE; p_RK is automatically updated
        {
        odeSolver->solve(runTime.timeOutputValue(), runTime.deltaTValue()+runTime.timeOutputValue(), p_RK, dtEst);
        }

        // Vapour temperature profile 
        // BC at the vapour liquid interface
        temp_v[0] = T_sat.value();

        // Roof boundary condition
        // Input: Choose: Neumann or Robin
        // If Neumann, dT(roof) is calculated from
        // an input parameter q_roof

        if(roof_bc == "Neumann")
        {
           temp_v[nz_v-1] = 1.0/3.0 * (4 * dz_v * q_roof.value()/k_V.value() + 4 * temp_v[nz_v-2] - temp_v[nz_v-3]);
        }
        else if (roof_bc == "Robin")
        {
           temp_v[nz_v-1] = (U_roof.value() * (1-eta_evap.value()) /k_V.value() * 2.0 * dz_v * T_air.value() +
            2.0 * temp_v[nz_v-2] - 1.0/2.0 * temp_v[nz_v-3]) /
           (3.0/2.0 + U_roof.value() * (1-eta_evap.value()) /k_V.value() * 2.0 * dz_v);
        }
        

        odeSolver_temp->solve(runTime.timeOutputValue(), runTime.deltaTValue()+runTime.timeOutputValue(), temp_v, dtEst);
        
        // rho_gas calculation
        if (r_q.value() < 1){
            // Calculate evaporation rate for non-equilibrium evaporation
            dimensionedScalar Q_VI(0);
            if (vapourHeatTransferModel.value() == 2)
                {
                    // Fixed r_q
                    Q_VI = r_q*Q_V;
                    // Update average vapour temperature
                    T_vap = p_RK[0] / (8.314 * rho_gas);
                }
            else if (vapourHeatTransferModel.value() == 3)
            {
                // convectionCoefficient model
                Q_VI = h_VL * A_T * (T_vap - T_sat);
                Info << "Q_VI convection = " << Q_VI.value() << endl;
                // Update average vapour temperature
                T_vap = p_RK[0] / (8.314 * rho_gas);

            }
            else if (vapourHeatTransferModel.value() == 4)
            {
                // Calculate Q_VI using Huerta and Vesovic (2019) analytical
                // solutions setting v_z = 0
                #include "analytical_QVI.H"
                // Calculate T_vap using Huerta and Vesovic (2019) analytical
                // solutions setting v_z = 0
                #include "analytical_T_vap.H"
            }
            // Vapour 1-D model
            else if (vapourHeatTransferModel.value() == 5)
            {
                // Calculate average vapour temperature
                T_vap = trapezoidalRule(temp_v, dz_v);
                // Vapour to liquid heat transfer rate
                Q_VI = k_V.value() * A_T.value() * (-3/2.0 * temp_v[0] + 2*temp_v[1] - 1/2.0*temp_v[2]) / dz_v;
                // Update vapour thermal conductivity and density
                k_V = kV_a0 + kV_a1 * T_vap + kV_a2 * pow(T_vap, 2) + kV_a3 * pow(T_vap, 3);
                k_V_roof = k_V;
                ode_temp.k_V = k_V;
                Info << "ode_temp_rho" << ode_temp.rho_gas.value() << endl;
                // Info << "Q_VL " << Q_VL_neq << "W" << endl;
                // Info << "Average vapour temperature using trapz = " << temp_average <<  " K " << endl;
            }


            // Add phase change through wall boiling and local superheating
            dmdt = (Q_VI + qL_int_flux + Q_L_wb)/deltaH + nvap_wall/runTime.deltaTValue();
            
            nvap_dt = dmdt * runTime.deltaTValue();
            // Update numbers of moles evaporated
            n_evaporated += nvap_dt;
            // Update vapour volume
            #include "liquidThermalExpansion.H"
            V_vap += nvap_dt/rho_L + (V_tank - V_vap) / rho_L * drho_L_dt.value() * runTime.deltaTValue(); 
            Info << "thermal exp term: " << (V_tank - V_vap) / rho_L * drho_L_dt.value() * runTime.deltaTValue() << " m^3 s^-1" << endl;

            // Update vapour density
            rho_gas = n_evaporated / V_vap;

            // Alternative mass balance calculation
            dimensionedScalar rho_gas_alt = rho_old + (dmdt / V_old) * (rho_L - rho_old) / rho_L * runTime.deltaTValue();


            // Alternative rho_Gas calculation
            // rho_gas = p_RK[0] / (8.314 * T_vap);
            Info << "Non-equilibrium evaporation rate = " << dmdt.value() << "mol s^-1" << endl;
            Info << "rho_gas_mb = " << rho_gas_alt.value() << "mol m^-3" << endl;
        }
        else {
            // Equilibrium routine
            // Update vapour temperature after pressure build-up
            T_sat = 264.651 / (3.7362 - Foam::log10(p_RK[0]/1e5) ) + 6.788;
            // Update vapour density after pressure build-up
            rho_gas = p_RK[0] / (8.314 * T_sat);
            // Calculate evaporation rate using first order forward differences
            dmdt = (rho_gas - rho_old)/runTime.deltaTValue() * V_old * (rho_L / (rho_L - rho_old) );
            // Add wall and surface evaporation rate
            // dmdt += nvap_wall/runTime.deltaTValue();
            // Calculate moles evaporated in a time-step
            nvap_dt = dmdt * runTime.deltaTValue();  
            // Acumulate vapour
            n_evaporated += dmdt;
            // Calculate rho_L and drho_L/dt
            #include "liquidThermalExpansion.H"
            // Update vapour volume considering liquid thermal expansion
            V_vap += nvap_dt/rho_L + (V_tank - V_vap) / rho_L * drho_L_dt.value() * runTime.deltaTValue(); 
            // Set vapour temperature equal to saturation owing to thermal equilibrium assumption
            T_vap = T_sat;
        }

        if(debug == true)
            {
                Info << "Equilibrium evaporation rate = " << dmdt.value() << "mol s^-1 " << endl;
                Info << "Mass balance debug" << endl;
                Info << "rho_old = " << rho_old.value() << "mol m^-3 " << endl;
                Info << "rho_new = " << rho_gas.value() << "mol m^-3 " << endl;
                Info << "V_old = " << V_old.value() << "m^-3 " << endl;
                Info << "V_new = " << V_vap.value() << "m^-3 " << endl;
            }


        Info << "number of moles evaporated = " << nvap_dt.value() << endl;
        Info << "Q_V =" << Q_V.value() << endl;

        // Output important variables of the solver
        Info << "pvap = " << p_RK << endl;
        Info << "T_sat = " << T_sat.value() << endl;
        Info << "T_vap = " << T_vap.value() << endl;
        Info << "rho_gas = " << rho_gas.value() << endl;
        Info << "T_new - T_old = " << T_new.value() - T_old.value() << endl;
        Info << "n_evaporated = " << n_evaporated.value() << endl;

        // Trial, delete this asap: update vapour volume dodgily
        // V_vap = V_V_exp;
 
        // Set the vapour pressure at the new pressure in all nodes
        volScalarField& pv = pvap;
        // This is the average vapour temperature
        volScalarField& Tv = Tvap;
        volScalarField& nvap = n_evap;

        forAll(pv, i)
        {
            pv[i] = p_RK[0];
            Tv[i] = T_vap.value();
            nvap[i] = n_evaporated.value();
            V_vapfield[i] = V_vap.value();
        }

        runTime.write();

        // Write vapour temperature
        // If data must be read in this time-step, based on 
        // writeinterval...
        
        if(runTime.outputTime())
        {
            // Define a fileName for the vapour tempreature object
            fileName outDir(runTime.path()/runTime.timeName());
            // Create OFstream
            OFstream os(outDir/"vap_temp.csv");
            // Write the array to os
            os << "Vapour Temperature" << temp_v << endl;
            // Mesh
            OFstream os2(outDir/"vap_mesh.csv");
            os2 << "Vapour mesh" << (vapour_mesh* V_vap.value()/A_T.value()) << endl;
        }        
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
