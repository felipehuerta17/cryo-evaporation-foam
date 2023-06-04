#include "Field.H"
#include "trapezoidal.H"
namespace Foam {
scalar trapezoidalRule (scalarField f, double dx) {
    scalar int_f (0);
    // Integrate over the mesh
    for (int i = 0; i < f.size()-1; i++){
        int_f += (f[i]+f[i+1])/2 * dx;
    }
// Divide by domain length
// Info << "f = " << int_f / (dx * (f.size()-1))<< endl;
return int_f / (dx * (f.size()-1)) ;
};
} // End namespace Foam