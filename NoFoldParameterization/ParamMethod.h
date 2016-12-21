#include <string>
#include "vect.h"
#include "array.h"
#include "timer.h"
#include <vector>

using namespace std;

class MatchingBoundary
{
public:
	int u1,u2,v1,v2;
	double p1,p2,p3,p4;
};

class MatchingV
{
public:
	int u1,v1;
};


class ParamMethod
{
public:
	string name;
	vector<vect3i> *ind_tex;

	double stime, etime;
	int iteration_count;

	ParamMethod(string n);

	vector<MatchingBoundary> match_boundary;
	//vector<MatchingV> match_v;

	double *full_mesh;
	bool is_full_created;

	virtual void set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it);
	virtual int run(int iters, double *pos, double *ans);
	void runHessian(int iters, double *pos, double *ans) {}
	virtual double calc_error(double* pos, double step, double* grad);
	virtual void calc_full_gradient(double* pos,  double *g );
	virtual void display_error( double *e, double *pos );
	virtual void print_timings();
	virtual void set_seamless_param(double t);
	virtual double max_parameter(double *f,
    double *x,
	double *gg,
    double *q,
	int n);
	//For barriers
	virtual void sharpen();
	virtual void sharpen_alpha();
	virtual void reduce_alpha();
	virtual void init_lbfgs(int n);
	virtual void is_movement(double dis, double grad);
};
