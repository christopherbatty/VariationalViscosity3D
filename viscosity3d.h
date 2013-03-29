#ifndef VISCOSITY3D_H
#define VISCOSITY3D_H

#include "array3.h"

void advance_viscosity_implicit_weighted(Array3f& u, Array3f& v, Array3f& w, 
                                         const Array3f& vol_u, const Array3f& vol_v, const Array3f& vol_w, 
                                         const Array3f& vol_c, const Array3f& vol_ex, const Array3f& vol_ey, const Array3f& vol_ez,
                                         const Array3f& solid_phi,
                                         const Array3f& viscosity, float dt, float dx);
#endif

