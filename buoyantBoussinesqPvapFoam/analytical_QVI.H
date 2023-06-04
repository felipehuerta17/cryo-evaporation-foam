// Calculate Q_VI = Q_VL using the analytical solution
Q_VI = constant::mathematical::pi * pow(d_i.value(), 2) / 4 *
        k_V.value() * (c_1 * chi_minus + c_2 * chi_plus);

Info << "Q_VI_analytical = " << Q_VI << " W " << endl;
