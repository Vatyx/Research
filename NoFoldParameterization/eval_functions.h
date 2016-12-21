//#include "BarrierMethod.h"

#include "iso_pre.h"
#include <vector>
#include <Eigen/Sparse>

//Jason Variant
// s1^2 + s2^2 + 1 / s1^2 + 1 / s2^2
//Scott Construction
double jason_function_tri(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void jason_function_grad(double *newgradient, double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);


double paper_jason_function_tri(iso_pre *precomps, double u0, double v0, double u1, double v1, double u2, double v2);
void paper_jason_function_grad(double *newgradient, iso_pre *precomps, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);

void paper_jason_function_grad(double *newgradient, iso_pre *precomps, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);
void paper_jason_function_grad_quad_nonlin(double *gradQuad, double *gradNonlin, iso_pre *precomps, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);
void paper_jason_function_hessian_quad(std::vector<Eigen::Triplet<double> > &hessianTrips, iso_pre *precomps, int i, int ip1, int j);


//Yarons
// s1^2 + 1 / s2^2
//Scott construction
//doesnt move
double yaron_function_tri(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void yaron_function_grad(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);
//Yaron construction
double yaron_function_tri2(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void yaron_function_grad2(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);


double yaron_function_tri_sqrt_int(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void yaron_function_grad_sqrt_int(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);



//Scott's Simplified Yaron's
//Function going negative....
//s1 + 1/s2
//Yaron Construction
double syaron_function_tri(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void syaron_function_grad(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);

//Conformal Error
//Function going negative
//s1 /s2
//Yaron Construction
double confromal_function_tri(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void conformal_function_grad(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);

//Mips error
//Function going negative
//s1 /s2 + s2 /s1
//Yaron Construction
double mips_function_tri(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void mips_function_grad(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);


double lscm_function_tri(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void lscm_function_grad(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);


double olga_function_tri(double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2);
void olga_function_grad(double *newgradient, double L, double x, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);



//double jason_function_tri(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2);
//void jason_function_grad(double *newgradient, double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j);