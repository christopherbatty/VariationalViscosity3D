#include "viscosity3d.h"
#include "array3.h"
//#include "sparse/sparsematrix.h"
//#include "sparse/cgsolver.h"
#include "pcgsolver/pcg_solver.h"
#include <fstream>
#include <cmath>

int u_ind(int i, int j, int k, int nx, int ny) {
   return i + j*(nx+1) + k*(nx+1)*ny;
}

int v_ind(int i, int j, int k, int nx, int ny, int nz){
   return i + j*nx + k*nx*(ny+1) + (nx+1)*ny*nz ;   
}

int w_ind(int i, int j, int k, int nx, int ny, int nz){
   return i + j*nx + k*nx*ny + (nx+1)*ny*nz + nx*(ny+1)*nz;   
}

Array3c u_state;//(nx+1,ny,nz,0);
Array3c v_state;//(nx,ny+1,nz,0);
Array3c w_state;//(nx,ny,nz+1,0);


SparseMatrixd matrix;//(dim,dim,15);
SparseMatrixd matrix2;//(dim,dim
std::vector<double> rhs;//(dim);
std::vector<double> soln;//(dim);

void advance_viscosity_implicit_weighted(Array3f& u, Array3f& v, Array3f& w, 
                                         const Array3f& vol_u, const Array3f& vol_v, const Array3f& vol_w, 
                                         const Array3f& vol_c, const Array3f& vol_ex, const Array3f& vol_ey, const Array3f& vol_ez,
                                         const Array3f& solid_phi,
                                         const Array3f& viscosity, float dt, float dx) {
   float over_dx = 1.0f/dx;
   int nx = solid_phi.ni;
   int ny = solid_phi.nj;
   int nz = solid_phi.nk;
   printf("Creating state arrays.\n");
   std::cout << "Phi-size:" << nx << " " << ny << " " << nz << std::endl;
   int dim = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1);
   if(u_state.ni != u.ni) {
      printf("Creating matrices and vectors.\n");
      u_state.resize(nx+1,ny,nz);
      v_state.resize(nx,ny+1,nz);
      w_state.resize(nx,ny,nz+1);
      matrix.resize(dim);
      matrix2.resize(dim);
      rhs.resize(dim);
      soln.resize(dim);
      printf("Done that for good.");
   }
   
   u_state.assign(0);
   v_state.assign(0);
   w_state.assign(0);
   matrix.zero();
   //matrix2.zero();
   rhs.assign(dim, 0);
   soln.assign(dim, 0);

   const int SOLID = 3;
   const int FLUID = 2;
   //const int AIR = 1;

   //check if interpolated velocity positions are inside solid
   for(int k = 0; k < nz; ++k) for(int j = 0; j < ny; ++j) for(int i = 0; i < nx+1; ++i) {
      if(i - 1 < 0 || i >= nx || solid_phi(i-1,j,k) + solid_phi(i,j,k) <= 0)
         u_state(i,j,k) = SOLID;
      else 
         u_state(i,j,k) = FLUID;
   }

   for(int k = 0; k < nz; ++k) for(int j = 0; j < ny+1; ++j) for(int i = 0; i < nx; ++i) {
      if(j - 1 < 0 || j >= ny || solid_phi(i,j-1,k) + solid_phi(i,j,k) <= 0)
         v_state(i,j,k) = SOLID;
      else 
         v_state(i,j,k) = FLUID;
   }

   for(int k = 0; k < nz+1; ++k) for(int j = 0; j < ny; ++j) for(int i = 0; i < nx; ++i) {
      if(k - 1 < 0 || k >= nz || solid_phi(i,j,k-1) + solid_phi(i,j,k) <= 0)
         w_state(i,j,k) = SOLID;
      else 
         w_state(i,j,k) = FLUID;
   }
   
   float factor = dt*sqr(over_dx);
   //u-terms
   //2u_xx+ v_xy +uyy + u_zz + w_xz
   printf("Building u-components.\n");
   for(int k = 1; k < nz; ++k) for(int j = 1; j < ny; ++j) for(int i = 1; i < nx; ++i) {
      
      if(u_state(i,j,k) == FLUID) {
         int index = u_ind(i,j,k,nx,ny);
         
         rhs[index] = vol_u(i,j,k)*u(i,j,k);
         matrix.set_element(index,index,vol_u(i,j,k));

         float visc_right = viscosity(i,j,k);
         float visc_left = viscosity(i-1,j,k);
         float vol_right = vol_c(i,j,k);
         float vol_left = vol_c(i-1,j,k);

         float visc_top = 0.25f*(viscosity(i-1,j+1,k) + viscosity(i-1,j,k) + viscosity(i,j+1,k) + viscosity(i,j,k));
         float visc_bottom = 0.25f*(viscosity(i-1,j,k) + viscosity(i-1,j-1,k) + viscosity(i,j,k) + viscosity(i,j-1,k));
         float vol_top = vol_ez(i,j+1,k);
         float vol_bottom = vol_ez(i,j,k);

         float visc_front = 0.25f*(viscosity(i-1,j,k+1) + viscosity(i-1,j,k) + viscosity(i,j,k+1) + viscosity(i,j,k));
         float visc_back = 0.25f*(viscosity(i-1,j,k) + viscosity(i-1,j,k-1) + viscosity(i,j,k) + viscosity(i,j,k-1));
         float vol_front = vol_ey(i,j,k+1);
         float vol_back = vol_ey(i,j,k);

         //u_x_right
         matrix.add_to_element(index,index, 2*factor*visc_right*vol_right);
         if(u_state(i+1,j,k) == FLUID)
            matrix.add_to_element(index,u_ind(i+1,j,k,nx,ny), -2*factor*visc_right*vol_right);
         else if(u_state(i+1,j,k) == SOLID)
            rhs[index] -= -2*factor*visc_right*vol_right*u(i+1,j,k);

         //u_x_left
         matrix.add_to_element(index,index, 2*factor*visc_left*vol_left);
         if(u_state(i-1,j,k) == FLUID)
            matrix.add_to_element(index,u_ind(i-1,j,k,nx,ny), -2*factor*visc_left*vol_left);
         else if(u_state(i-1,j,k) == SOLID)
            rhs[index] -= -2*factor*visc_left*vol_left*u(i-1,j,k);

         //u_y_top
         matrix.add_to_element(index,index, +factor*visc_top*vol_top);
         if(u_state(i,j+1,k) == FLUID)
            matrix.add_to_element(index,u_ind(i,j+1,k,nx,ny), -factor*visc_top*vol_top);
         else if(u_state(i,j+1,k) == SOLID)
            rhs[index] -= -u(i,j+1,k)*factor*visc_top*vol_top;

         //u_y_bottom
         matrix.add_to_element(index,index, +factor*visc_bottom*vol_bottom);
         if(u_state(i,j-1,k) == FLUID)
            matrix.add_to_element(index,u_ind(i,j-1,k,nx,ny), -factor*visc_bottom*vol_bottom);
         else if(u_state(i,j-1,k) == SOLID)
            rhs[index] -= -u(i,j-1,k)*factor*visc_bottom*vol_bottom;

         //u_z_front
         matrix.add_to_element(index,index, +factor*visc_front*vol_front);
         if(u_state(i,j,k+1) == FLUID)
            matrix.add_to_element(index,u_ind(i,j,k+1,nx,ny), -factor*visc_front*vol_front);
         else if(u_state(i,j,k+1) == SOLID)
            rhs[index] -= -u(i,j,k+1)*factor*visc_front*vol_front;

         //u_z_back
         matrix.add_to_element(index,index, +factor*visc_back*vol_back);
         if(u_state(i,j,k-1) == FLUID)
            matrix.add_to_element(index,u_ind(i,j,k-1,nx,ny), -factor*visc_back*vol_back);
         else if(u_state(i,j,k-1) == SOLID)
            rhs[index] -= -u(i,j,k-1)*factor*visc_back*vol_back;

         //v_x_top
         if(v_state(i,j+1,k) == FLUID)
            matrix.add_to_element(index,v_ind(i,j+1,k,nx,ny,nz), -factor*visc_top*vol_top);
         else if(v_state(i,j+1,k) == SOLID)
            rhs[index] -= -v(i,j+1,k)*factor*visc_top*vol_top;
         
         if(v_state(i-1,j+1,k) == FLUID)
            matrix.add_to_element(index,v_ind(i-1,j+1,k,nx,ny,nz), factor*visc_top*vol_top);
         else if(v_state(i-1,j+1,k) == SOLID)
            rhs[index] -= v(i-1,j+1,k)*factor*visc_top*vol_top;
         
         //v_x_bottom
         if(v_state(i,j,k) == FLUID)
            matrix.add_to_element(index,v_ind(i,j,k,nx,ny,nz), +factor*visc_bottom*vol_bottom);
         else if(v_state(i,j,k) == SOLID)
            rhs[index] -= v(i,j,k)*factor*visc_bottom*vol_bottom;
         
         if(v_state(i-1,j,k) == FLUID)
            matrix.add_to_element(index,v_ind(i-1,j,k,nx,ny,nz), -factor*visc_bottom*vol_bottom);
         else if(v_state(i-1,j,k) == SOLID)
            rhs[index] -= -v(i-1,j,k)*factor*visc_bottom*vol_bottom;
         
         //w_x_front
         if(w_state(i,j,k+1) == FLUID)
            matrix.add_to_element(index,w_ind(i,j,k+1,nx,ny,nz), -factor*visc_front*vol_front);
         else if(w_state(i,j,k+1) == SOLID)
            rhs[index] -= -w(i,j,k+1)*factor*visc_front*vol_front;
         
         if(w_state(i-1,j,k+1) == FLUID)
            matrix.add_to_element(index,w_ind(i-1,j,k+1,nx,ny,nz), factor*visc_front*vol_front);
         else if(w_state(i-1,j,k+1) == SOLID)
            rhs[index] -= w(i-1,j,k+1)*factor*visc_front*vol_front;
         
         //w_x_back
         if(w_state(i,j,k) == FLUID)
            matrix.add_to_element(index,w_ind(i,j,k,nx,ny,nz), +factor*visc_back*vol_back);
         else if(w_state(i,j,k) == SOLID)
            rhs[index] -= w(i,j,k)*factor*visc_back*vol_back;
         
         if(w_state(i-1,j,k) == FLUID)
            matrix.add_to_element(index,w_ind(i-1,j,k,nx,ny,nz), -factor*visc_back*vol_back);
         else if(w_state(i-1,j,k) == SOLID)
            rhs[index] -= -w(i-1,j,k)*factor*visc_back*vol_back;
      }
   }

   //v-terms
   //vxx + 2vyy + vzz + u_yx + w_yz
   printf("Building v-components.\n");
   for(int k = 1; k < nz; ++k) for(int j = 1; j < ny; ++j) for(int i = 1; i < nx; ++i) {
      if(v_state(i,j,k) == FLUID) {
         int index = v_ind(i,j,k,nx,ny,nz);      
         
         rhs[index] = vol_v(i,j,k)*v(i,j,k);
         matrix.set_element(index, index, vol_v(i,j,k));

         float visc_right = 0.25f*(viscosity(i,j-1,k) + viscosity(i+1,j-1,k) + viscosity(i,j,k) + viscosity(i+1,j,k));
         float visc_left = 0.25f*(viscosity(i,j-1,k) + viscosity(i-1,j-1,k) + viscosity(i,j,k) + viscosity(i-1,j,k));
         float vol_right = vol_ez(i+1,j,k);
         float vol_left = vol_ez(i,j,k);
         
         float visc_top = viscosity(i,j,k);
         float visc_bottom = viscosity(i,j-1,k);
         float vol_top = vol_c(i,j,k);
         float vol_bottom = vol_c(i,j-1,k);
         
         float visc_front = 0.25f*(viscosity(i,j-1,k) + viscosity(i,j-1,k+1) + viscosity(i,j,k) + viscosity(i,j,k+1));
         float visc_back = 0.25f*(viscosity(i,j-1,k) + viscosity(i,j-1,k-1) + viscosity(i,j,k) + viscosity(i,j,k-1));
         float vol_front = vol_ex(i,j,k+1);
         float vol_back = vol_ex(i,j,k);
         
         //v_x_right
         matrix.add_to_element(index,index, +factor*visc_right*vol_right);
         if(v_state(i+1,j,k) == FLUID)
            matrix.add_to_element(index,v_ind(i+1,j,k,nx,ny,nz), -factor*visc_right*vol_right);
         else if(v_state(i+1,j,k) == SOLID)
            rhs[index] -= -v(i+1,j,k)*factor*visc_right*vol_right;
      
         //v_x_left
         matrix.add_to_element(index,index, +factor*visc_left*vol_left);
         if(v_state(i-1,j,k) == FLUID)
            matrix.add_to_element(index,v_ind(i-1,j,k,nx,ny,nz), -factor*visc_left*vol_left);
         else if(v_state(i-1,j,k) == SOLID)
            rhs[index] -= -v(i-1,j,k)*factor*visc_left*vol_left;

         //vy_top
         matrix.add_to_element(index,index, +2*factor*visc_top*vol_top);
         if(v_state(i,j+1,k) == FLUID)
            matrix.add_to_element(index,v_ind(i,j+1,k,nx,ny,nz), -2*factor*visc_top*vol_top);
         else if (v_state(i,j+1,k) == SOLID)
            rhs[index] -= -2*factor*visc_top*vol_top*v(i,j+1,k);
         
         //vy_bottom
         matrix.add_to_element(index,index, +2*factor*visc_bottom*vol_bottom);
         if(v_state(i,j-1,k) == FLUID)
            matrix.add_to_element(index,v_ind(i,j-1,k,nx,ny,nz), -2*factor*visc_bottom*vol_bottom);
         else if(v_state(i,j-1,k) == SOLID)
            rhs[index] -= -2*factor*visc_bottom*vol_bottom*v(i,j-1,k);

         //v_z_front
         matrix.add_to_element(index,index, +factor*visc_front*vol_front);
         if(v_state(i,j,k+1) == FLUID)
            matrix.add_to_element(index,v_ind(i,j,k+1,nx,ny,nz), -factor*visc_front*vol_front);
         else if(v_state(i+1,j,k) == SOLID)
            rhs[index] -= -v(i,j,k+1)*factor*visc_front*vol_front;
         
         //v_z_back
         matrix.add_to_element(index,index, +factor*visc_back*vol_back);
         if(v_state(i,j,k-1) == FLUID)
            matrix.add_to_element(index,v_ind(i,j,k-1,nx,ny,nz), -factor*visc_back*vol_back);
         else if(v_state(i,j,k-1) == SOLID)
            rhs[index] -= -v(i,j,k-1)*factor*visc_back*vol_back;
         
         //u_y_right
         if(u_state(i+1,j,k) == FLUID)
            matrix.add_to_element(index,u_ind(i+1,j,k,nx,ny), -factor*visc_right*vol_right);
         else if(u_state(i+1,j,k) == SOLID)
            rhs[index] -= -u(i+1,j,k)*factor*visc_right*vol_right;
         
         if(u_state(i+1,j-1,k) == FLUID)
            matrix.add_to_element(index,u_ind(i+1,j-1,k,nx,ny), factor*visc_right*vol_right);
         else if(u_state(i+1,j-1,k) == SOLID)
            rhs[index] -= u(i+1,j-1,k)*factor*visc_right*vol_right;
      
         //u_y_left
         if(u_state(i,j,k) == FLUID)
            matrix.add_to_element(index,u_ind(i,j,k,nx,ny), factor*visc_left*vol_left);
         else if(u_state(i,j,k) == SOLID)
            rhs[index] -= u(i,j,k)*factor*visc_left*vol_left;
          
         if(u_state(i,j-1,k) == FLUID)
            matrix.add_to_element(index,u_ind(i,j-1,k,nx,ny), -factor*visc_left*vol_left);
         else if(u_state(i,j-1,k) == SOLID)
            rhs[index] -= -u(i,j-1,k)*factor*visc_left*vol_left;
      
         //w_y_front
         if(w_state(i,j,k+1) == FLUID)
            matrix.add_to_element(index,w_ind(i,j,k+1,nx,ny,nz), -factor*visc_front*vol_front);
         else if(w_state(i,j,k+1) == SOLID)
            rhs[index] -= -w(i,j,k+1)*factor*visc_front*vol_front;
         
         if(w_state(i,j-1,k+1) == FLUID)
            matrix.add_to_element(index,w_ind(i,j-1,k+1,nx,ny,nz), factor*visc_front*vol_front);
         else if(w_state(i,j-1,k+1) == SOLID)
            rhs[index] -= w(i,j-1,k+1)*factor*visc_front*vol_front;
         
         //w_y_back
         if(w_state(i,j,k) == FLUID)
            matrix.add_to_element(index,w_ind(i,j,k,nx,ny,nz), factor*visc_back*vol_back);
         else if(w_state(i,j,k) == SOLID)
            rhs[index] -= w(i,j,k)*factor*visc_back*vol_back;
          
         if(w_state(i,j-1,k) == FLUID)
            matrix.add_to_element(index,w_ind(i,j-1,k,nx,ny,nz), -factor*visc_back*vol_back);
         else if(w_state(i,j-1,k) == SOLID)
            rhs[index] -= -w(i,j-1,k)*factor*visc_back*vol_back;
      }
   }

   //w-terms
   //wxx+ wyy+ 2wzz + u_zx + v_zy
   printf("Building w-components.\n");
   for(int k = 1; k < nz; ++k) for(int j = 1; j < ny; ++j) for(int i = 1; i < nx; ++i) {
      if(w_state(i,j,k) == FLUID) {
         int index = w_ind(i,j,k,nx,ny,nz);
         rhs[index] = vol_w(i,j,k)*w(i,j,k);
         matrix.set_element(index,index, vol_w(i,j,k));

         float visc_right = 0.25f*(viscosity(i,j,k) + viscosity(i,j,k-1) + viscosity(i+1,j,k) + viscosity(i+1,j,k-1));
         float visc_left = 0.25f*(viscosity(i,j,k) + viscosity(i,j,k-1) + viscosity(i-1,j,k) + viscosity(i-1,j,k-1));
         float vol_right = vol_ey(i+1,j,k);
         float vol_left = vol_ey(i,j,k);

         float visc_top = 0.25f*(viscosity(i,j,k) + viscosity(i,j,k-1) + viscosity(i,j+1,k) + viscosity(i,j+1,k-1));;
         float visc_bottom = 0.25f*(viscosity(i,j,k) + viscosity(i,j,k-1) + viscosity(i,j-1,k) + viscosity(i,j-1,k-1));;
         float vol_top = vol_ex(i,j+1,k);
         float vol_bottom = vol_ex(i,j,k);

         float visc_front = viscosity(i,j,k);   
         float visc_back = viscosity(i,j,k-1); 
         float vol_front = vol_c(i,j,k);
         float vol_back = vol_c(i,j,k-1);

         //w_x_right
         matrix.add_to_element(index,index, +factor*visc_right*vol_right);
         if(w_state(i+1,j,k) == FLUID)
            matrix.add_to_element(index,w_ind(i+1,j,k,nx,ny,nz), -factor*visc_right*vol_right);
         else if (w_state(i+1,j,k) == SOLID)
            rhs[index] -= -factor*visc_right*vol_right*w(i+1,j,k);
         
         //w_x_left
         matrix.add_to_element(index,index, factor*visc_left*vol_left);
         if(w_state(i-1,j,k) == FLUID)
            matrix.add_to_element(index,w_ind(i-1,j,k,nx,ny,nz), -factor*visc_left*vol_left);
         else if (w_state(i-1,j,k) == SOLID)
            rhs[index] -= -factor*visc_left*vol_left*w(i-1,j,k);

         //w_y_top
         matrix.add_to_element(index,index, +factor*visc_top*vol_top);
         if(w_state(i,j+1,k) == FLUID)
            matrix.add_to_element(index,w_ind(i,j+1,k,nx,ny,nz), -factor*visc_top*vol_top);
         else if (w_state(i,j+1,k) == SOLID)
            rhs[index] -= -factor*visc_top*vol_top*w(i,j+1,k);
         
         //w_y_bottom
         matrix.add_to_element(index,index, factor*visc_bottom*vol_bottom);
         if(w_state(i,j-1,k) == FLUID)
            matrix.add_to_element(index,w_ind(i,j-1,k,nx,ny,nz), -factor*visc_bottom*vol_bottom);
         else if (w_state(i,j-1,k) == SOLID)
            rhs[index] -= -factor*visc_bottom*vol_bottom*w(i,j-1,k);
         
         //w_z_front
         matrix.add_to_element(index,index, +2*factor*visc_front*vol_front);
         if(w_state(i,j,k+1) == FLUID)
            matrix.add_to_element(index,w_ind(i,j,k+1,nx,ny,nz), -2*factor*visc_front*vol_front);
         else if (w_state(i,j,k+1) == SOLID)
            rhs[index] -= -2*factor*visc_front*vol_front*w(i,j,k+1);
         
         //w_z_back
         matrix.add_to_element(index,index, +2*factor*visc_back*vol_back);
         if(w_state(i,j,k-1) == FLUID)
            matrix.add_to_element(index,w_ind(i,j,k-1,nx,ny,nz), -2*factor*visc_back*vol_back);
         else if (w_state(i,j,k-1) == SOLID)
            rhs[index] -= -2*factor*visc_back*vol_back*w(i,j,k-1);

         //u_z_right
         if(u_state(i+1,j,k) == FLUID)
            matrix.add_to_element(index,u_ind(i+1,j,k,nx,ny), -factor*visc_right*vol_right);
         else if(u_state(i+1,j,k) == SOLID)
            rhs[index] -= -u(i+1,j,k)*factor*visc_right*vol_right;
         
         if(u_state(i+1,j,k-1) == FLUID)
            matrix.add_to_element(index,u_ind(i+1,j,k-1,nx,ny), factor*visc_right*vol_right);
         else if(u_state(i+1,j,k-1) == SOLID)
            rhs[index] -= u(i+1,j,k-1)*factor*visc_right*vol_right;
      
         //u_z_left
         if(u_state(i,j,k) == FLUID)
            matrix.add_to_element(index,u_ind(i,j,k,nx,ny), factor*visc_left*vol_left);
         else if(u_state(i,j,k) == SOLID)
            rhs[index] -= u(i,j,k)*factor*visc_left*vol_left;
          
         if(u_state(i,j,k-1) == FLUID)
            matrix.add_to_element(index,u_ind(i,j,k-1,nx,ny), -factor*visc_left*vol_left);
         else if(u_state(i,j,k-1) == SOLID)
            rhs[index] -= -u(i,j,k-1)*factor*visc_left*vol_left;
         
         //v_z_top
         if(v_state(i,j+1,k) == FLUID)
            matrix.add_to_element(index,v_ind(i,j+1,k,nx,ny,nz), -factor*visc_top*vol_top);
         else if(v_state(i,j+1,k) == SOLID)
            rhs[index] -= -v(i,j+1,k)*factor*visc_top*vol_top;
         
         if(v_state(i,j+1,k-1) == FLUID)
            matrix.add_to_element(index,v_ind(i,j+1,k-1,nx,ny,nz), factor*visc_top*vol_top);
         else if(v_state(i,j+1,k-1) == SOLID)
            rhs[index] -= v(i,j+1,k-1)*factor*visc_top*vol_top;
     
         //v_z_bottom
         if(v_state(i,j,k) == FLUID)
            matrix.add_to_element(index,v_ind(i,j,k,nx,ny,nz), +factor*visc_bottom*vol_bottom);
         else if(v_state(i,j,k) == SOLID)
            rhs[index] -= v(i,j,k)*factor*visc_bottom*vol_bottom;
         
         if(v_state(i,j,k-1) == FLUID)
            matrix.add_to_element(index,v_ind(i,j,k-1,nx,ny,nz), -factor*visc_bottom*vol_bottom);
         else if(v_state(i,j,k-1) == SOLID)
            rhs[index] -= -v(i,j,k-1)*factor*visc_bottom*vol_bottom;

      }
   }

   ////strip out near zero entries to speed this thing up!
   //printf("Stripping out near-zeros.\n");
   ////SparseMatrixd matrix2(matrix.m,matrix.n,15);
   //for(unsigned int row = 0; row < matrix.m; ++row){
   //   for(unsigned int col = 0; col < matrix.index[row].size(); ++col) {
   //      int index = matrix.index[row][col];
   //      double val = matrix(row,index);
   //      if(std::abs(val) > 1e-10)
   //         matrix2.set_element(row,index, val);
   //   }
   //}
   printf("Solving sparse system.\n");
   PCGSolver<double> solver;
   double res_out;
   int iter_out;
   solver.set_solver_parameters(1e-7, 500);

   printf("Launching CG\n");
   solver.solve(matrix, rhs, soln, res_out, iter_out);

   if(iter_out >= 500)
      printf("\n\n\n**********FAILED**************\n\n\n");
   
   printf("Copying back.\n");
   for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
         for(int i = 0; i < nx+1; ++i) {
            if(u_state(i,j,k) == FLUID) {
               u(i,j,k) = (float)soln[u_ind(i,j,k,nx,ny)];
            }
            else {
               u(i,j,k) = 0;
            }
         }

   for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny+1; ++j)
         for(int i = 0; i < nx; ++i) {
            if(v_state(i,j,k) == FLUID) {
               v(i,j,k) = (float)soln[v_ind(i,j,k,nx,ny,nz)];
            }
            else
               v(i,j,k) = 0;
         }

   for(int k = 0; k < nz+1; ++k)
      for(int j = 0; j < ny; ++j)
         for(int i = 0; i < nx; ++i) {
            if(w_state(i,j,k)== FLUID) {
               w(i,j,k) = (float)soln[w_ind(i,j,k,nx,ny,nz)];
            }
            else
               w(i,j,k) = 0;
         }
   printf("Done copying back.\n");
   
}

