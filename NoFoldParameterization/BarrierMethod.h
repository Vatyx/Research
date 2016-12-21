#include "ParamMethod.h"
#include "PointGrid.h"
#include "ResizeArray.h"
#include <vector>
#include <float.h>
#include "timer.h"

#include "iso_pre.h"

//#define SCOTTS_CODE
#define DO_FULL_TIMINGS

#define EPSION_FUNC 1 //Slope
//#define CONST_CLAMP .04  //Distance away Use this for paper
#define CONST_CLAMP .4  //Distance away Use this for paper

//Min Iso Boundary Edge
//camel
//#define CONST_CLAMP .37  //Distance away Use this for paper
//cow
//#define CONST_CLAMP .1825  //Distance away Use this for paper
//tricera
//#define CONST_CLAMP .173
//hose
//#define CONST_CLAMP .077



//#define CONST_CLAMP 10  //Distance away
#define BBOX_INTER_THRES CONST_CLAMP*1.01

//
//#define EPSION_FUNC 4 //Slope
//#define CONST_CLAMP .0004   //Distance away
#define STARTING_EXP 2 //Starting sharpening exponent for func and gradient
#define POINT_GRID_DISC 200 //Grid discritization for boundary tests

#define GOLDEN_EPS .0000000001
#define SA_COEFF .999

#define USE_INTERVAL_TREE
#include "IntervalTree3.h"
#define INTERVAL_TREE_SIZE 200

//#define PERFORM_SEAMLESS
//#define SEAMLESS_WITHOUT_BOUND

//#define TURN_OFF_BOUNDARY
#define ALPHA_INC .01

//#define GLOBAL_BOUNDARY_SUPPORT
//#define ALPHA_TIMES_X
#ifdef ALPHA_TIMES_X
#define ALPHA .05
#else
#define ALPHA 1
#endif

#define NUM_PROCS 4
using namespace std;
	

class BarrierMethod : public ParamMethod
{
public:
	double func_exp;
	double currentsA;
	double seamless_exp;

	double alpha_val;
	PointGrid<POINT_GRID_DISC> pg;

#ifdef USE_INTERVAL_TREE
	IntervalTree<INTERVAL_TREE_SIZE> tree;
#endif

	double time_sA_intervalsetup,time_sA_sort,time_sA_bound;
	double time_func_bound, time_grad_bound, time_funcgrad_intervalsetup;
	
	int total_tex_size;
	int tex_size;//# of tex coors * dimension
	vector<int> *indchange;
	vector<vector<int>> *boundaryorder;
	vector<vector<vect2d>> *iso_tris;

	double calc_seamless_error(double a, double b);
	double calc_seamless_error(double u1, double v1, double u2, double v2);
	vector<iso_pre> precomp;
	
	ResizeArray<int> int_pos[NUM_PROCS];
	ResizeArray<int> in_pos;

	

	//Used for max param
	vector<std::pair<int, double>> mins_interval;
	vector<std::pair<int, double>> pmins_interval;
	vector<double> maxs_interval;
	vector<double> pmaxs_interval;
	void change_sA(double &sA, double *x, double *q, int i, int ip1, int j);
	//max param

	BarrierMethod(string n, int num_uv, vector<vector<int>> bo, vector<vector<vect2d>> it);
	void set_seamless_param(double t);

//Only callable functions needed
	void set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it);
	void calc_boundary_gradient( double* pos, double *g );
	int run(int iters, double *pos, double *ans);
	void runHessian(int iters, double *pos, double *ans) {}

	double calc_error(double* pos, double step, double* grad);
	void calc_full_gradient( double* pos, double *g );
	void display_error(double *e, double *pos);
	virtual void calc_interior_gradient( double *pos, double *g );
	void print_timings();
	//Increase exponent on function by power of 2
	void sharpen();
	void sharpen_alpha();
	void reduce_alpha();
	void init_lbfgs(int n);
	void is_movement(double dis, double grad);

	double max_parameter(double *f,
    double *x,
	double *gg,
    double *q,
	int n);


//Mayeb make these private!!!!!!
	//Global Function Evaluation
	double global_barrier_eval(vect2d a, vect2d b, vect2d x);
	void global_barrier_eval(double &func, const double *pos, int i, int ip1,  int j );
	void global_barrier_eval(double &func, const double *pos, const double *grad, double step, int i, int ip1,  int j );

	void global_barrier_grad(double *newgradient, const double *pos, int i, int ip1,  int j );

	void fill_point_grid(const double *pos);
	double perform_full_global_barrier_func(double ff, const double *pos);
	
	//Call These Two
	//double perform_global_fold_barrier_func(double ff, const double *pos);
	//double perform_global_fold_barrier_func(double ff, const double *pos, double step, const double *grad);

	virtual double max_param( double *f, double *x, double *gg,  double *q, int n );
};
