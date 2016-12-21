#include "ParamMethod.h"

#include <iostream>

/*************************************************/
//Constructors
/*************************************************/
ParamMethod::ParamMethod(string n)
{ 
	name = n; 
	is_full_created = false;
}

/*************************************************/
//Virtual Functions
/*************************************************/
void ParamMethod::init_lbfgs(int n)
{
}

void ParamMethod::print_timings()
{
	//printf("Timings\n");
}
int ParamMethod::run(int iters, double *pos, double *ans)
{
	iteration_count = 0;
	//cout << "Run method: "<<   name << '\n';

	return 0;
}

double ParamMethod::calc_error(double* pos, double step, double* grad)
{
	//cout << name << "Error = ";
	return 0;
}

void ParamMethod::calc_full_gradient(double* pos,  double *g )
{
}

void ParamMethod::set_seamless_param(double t)
{
}
double ParamMethod::max_parameter(double *f,
    double *x,
	double *gg,
    double *q,
	int n)
{
	return 0;
}

void ParamMethod::display_error( double *e, double *pos )
{
}
void ParamMethod::is_movement(double dis, double grad)
{
}
void ParamMethod::set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it)
{
	this->ind_tex = face_ind;
//	printf("\nnum faces = %d\n", this->ind_tex->size());

	/*this->ind_tex.clear();

	for (int i = 0; i < face_ind.size(); i++)
	{
		vect3i t;
		t[0] = face_ind[i][0];
		t[1] = face_ind[i][1];
		t[2] = face_ind[i][2];

		ind_tex.push_back(t);
	}*/
}

void ParamMethod::sharpen(){}
void ParamMethod::sharpen_alpha(){}
void ParamMethod::reduce_alpha(){}
/*************************************************/
//Functions
/*************************************************/
