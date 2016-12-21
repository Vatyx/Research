#pragma once
#include "BarrierMethod.h"

using namespace std;
#include <string>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Core>
#include <Eigen/PardisoSupport>

#define MAX_ITERS 5
#define MPARAM 5

#define CHECK_BOUNDARY_SA
#define CULL_EDGE_CHECKS
#define FASTER_INTERVALS

#define USE_LBFGS_LINESEARCH

//#define ZERO_MAX_TRIS
#define GRAD_MAG 1
#define DIS_MAG .0000001

//#define CONFORMAL
//#define MIPS
//#define OLGAS
//#define YARONS
//#define YARONSSQRT

#define JASONS_VARIANT


#define JASON_PRE_COMP

//#define SYARONS



//#define LSCM


extern int  max_iters_t;
#include "eval_functions.h"

//typedef double (*eval_tri_func)(
//	double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2
//	);
//
//typedef void (*eval_gradient)(double *pos, double *newgradient);

class JasonFull : public BarrierMethod
{
public:

	bool firstTime;// = true;
	double **threadGrad;
	double *m_x;
	Eigen::SparseMatrix<double> hessian;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > hessianSolver;
	Eigen::SparseMatrix<double> kvf;
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > kvfSolver;
	Eigen::PardisoLLT<Eigen::SparseMatrix<double> > kvfSolver;


	double *temp_pos;
	//Array<vect3i> ind_tex;
	double time_sA_int,time_line_search,time_golden,time_func_int,time_grad_int;
	//JasonFull(string n, bool isboundary):ParamMethod(n, isboundary);
	JasonFull(string n, int num_uv, vector<vector<int>> bo,  vector<vector<vect2d>> it);
	void set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it);

	int run(int iters, double *pos, double *ans);
	double calc_error(double* pos, double step, double* grad);
	void calc_full_gradient( double *pos, double *g );
	void calc_interior_gradient( double *pos, double *g );
	void display_error( double *e, double *pos );
	void print_timings();

	void init_lbfgs(int n);
	
double JasonFull::max_param_interior(
	double *f,
    double *x,
	double *gg,
    double *q,
	int n
	);
	//lbfgs interface
	double evaluate(double *x, double *gg, int n,double step );
	int progress(double *x, double *g, double fx,
				 double xnorm, double gnorm, double step,
				 int n, int k, int ls
                );
	double max_param( double *f, double *x,
							double *gg,  double *q, int n );

	double golden_ratio(double x1, double y1, double x2, double y2, double x3, double y3, double *x, double *q);

	
	double max_parameter(double *f,
						double *x,
						double *gg,
						double *q,
						int n);

	void is_movement(double dis, double grad);

	double barrier_jason_function_tri(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2);
	double barrier_jason_function(double *pos);
	void barrier_jason_gradient(double *pos, double *newgradient);
	void jason_gradient_quad_nonlin(double *pos, double *quad, double *nonlin);

	void calc_quad_nonlin_gradient(double *pos, double *quad, double *nonlin);

	double (*e_tri_func) (double , double, double, double, double, double, double, double, double);
	void (*e_gradient) (double*, double , double, double, double, double, double, double, double, double, int, int, int);

	void jason_precompute_quad_hessian(void);
	void runHessian(int iters, double *pos, double *ans);
	void runKVF(int iters, double *pos, double *ans);

	void computeKVF_First_Time(double *x);
	void computeKVF(double *x);
	void localKVF ( vector<Eigen::Triplet<double> > &trips, double *x, int i0, int i1, int i2 );
	void localKVFAdd ( double *x, int i0, int i1, int i2 );

	//<index, distance from start>
	vector<pair<unsigned int, unsigned int>> diagIndexes;
	Eigen::SparseMatrix<double, 0, int>::Scalar** addresses;
	unsigned long addressCounter;

	void doStoreAddresses(int i1, int i2, int i3 );
	void storeAddresses();
	void doQuadrantAddresses(int i, int j);
	void storeDiagIndexes();

	void innerFlapOptimization(double* x);
};


