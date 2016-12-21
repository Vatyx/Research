#include "BarrierMethod.h"

using namespace std;
#include <string>

#define MAX_ITERS 1
#define MPARAM 3

#define CHECK_BOUNDARY_SA
#define CULL_EDGE_CHECKS
#define FASTER_INTERVALS

#define USE_LBFGS_LINESEARCH

//#define ZERO_MAX_TRIS
#define GRAD_MAG 1
#define DIS_MAG .000001

#define JASONS_VARIANT
//#define YARONS
//#define SYARONS
//#define CONFORMAL
//#define MIPS


//#define LSCM

#include "eval_functions.h"

//typedef double (*eval_tri_func)(
//	double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2
//	);
//
//typedef void (*eval_gradient)(double *pos, double *newgradient);

class JasonSingle : public BarrierMethod
{
public:

	bool firstTime;// = true;
	double **threadGrad;

	double *temp_pos;
	//Array<vect3i> ind_tex;
	double time_sA_int,time_line_search,time_golden,time_func_int,time_grad_int;
	//JasonFull(string n, bool isboundary):ParamMethod(n, isboundary);
	JasonSingle(string n, int num_uv, vector<vector<int>> bo,  vector<vector<vect2d>> it);
	void set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it);

	int run(int iters, double *pos, double *ans);
	double calc_error(double* pos, double step, double* grad);
	void calc_full_gradient( double *pos, double *g );
	void calc_interior_gradient( double *pos, double *g );
	void display_error( double *e, double *pos );
	void print_timings();
	
double max_param_interior(
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

	double (*e_tri_func) (double , double, double, double, double, double, double, double, double);
	void (*e_gradient) (double*, double , double, double, double, double, double, double, double, double, int, int, int);
};


