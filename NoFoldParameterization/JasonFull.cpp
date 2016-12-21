#include "JasonFull.h"
#include "lbfgs.h"
#include <omp.h>
#include "global.h"

#include <iostream>


int max_iters_t = MAX_ITERS;

#ifdef ZERO_MAX_TRIS
int max_param_tri;
vector<int> zero_tris;
#endif
/*************************************************/
//Constructors
/*************************************************/



#define PARALLEL_GRAD

JasonFull::JasonFull(string n, int num_uv, vector<vector<int>> bo, vector<vector<vect2d>> it)
	: BarrierMethod(n,num_uv, bo,it)
{
#ifdef ZERO_MAX_TRIS
	max_param_tri = 0;
	zero_tris.clear();
#endif
	firstTime = true;

	temp_pos = new double[tex_size*2];
	time_sA_int = time_line_search = time_golden = time_func_int = time_grad_int = 0;
#ifdef JASONS_VARIANT
	e_tri_func = jason_function_tri;
	e_gradient = jason_function_grad;
#endif

#ifdef YARONS
	e_tri_func = yaron_function_tri2;
	e_gradient = yaron_function_grad2;
#endif

#ifdef YARONSSQRT
	e_tri_func = yaron_function_tri_sqrt_int;
	e_gradient = yaron_function_grad_sqrt_int;
#endif

#ifdef SYARONS
	e_tri_func = syaron_function_tri;
	e_gradient = syaron_function_grad;
#endif

#ifdef CONFORMAL
	e_tri_func = confromal_function_tri;
	e_gradient = conformal_function_grad;
#endif

#ifdef MIPS
	e_tri_func = mips_function_tri;
	e_gradient = mips_function_grad;
#endif

#ifdef LSCM
	e_tri_func = lscm_function_tri;
	e_gradient = lscm_function_grad;
#endif

#ifdef OLGAS
	e_tri_func = olga_function_tri;
	e_gradient = olga_function_grad;
#endif
}

/*************************************************/
//RUN Functions
/*************************************************/

double ggR = 2 - (1 + sqrt((double)5))/(double)2.0;
double JasonFull::golden_ratio(double x1, double y1, double x2, double y2, double x3, double y3, double *x, double *q)
{
	double x4 = 0;
	if(x3 - x2 > x2 - x1)
		x4 = x2 + ggR*(x3-x2);
	else
		x4 = x2 - ggR*(x2-x1);

	if(x3-x1 < GOLDEN_EPS*(x2+x4))
		return (x1+x3)/(double)2.0;

	double y4 = 0;
	y4 = calc_error(x,-1*x4,q);

	if(y4 < y2)
	{
		if(x3 - x2 > x2 - x1)
			return golden_ratio(x2,y2, x4,y4, x3,y3,x,q);
		else
			return golden_ratio(x1,y1, x4,y4, x2,y2,x,q);
	}
	else
	{
		if(x3 - x2 > x2 - x1)
			return golden_ratio(x1,y1, x2,y2, x4,y4,x,q);
		else
			return golden_ratio(x4,y4, x2,y2, x3,y3,x,q);

	}
}				

void JasonFull::print_timings()
{
	
	printf("Interior\n");
	printf("time for max param = %f\n", time_sA_int);
	printf("\t interior %f\n", time_sA_int);
	//printf("\t line %f\n", time_line_search);
	//printf("\t golden %f\n", time_golden);
	printf("time for func and grad eval = %f\n",time_func_int+time_grad_int);
	printf("\t func %f\n", time_func_int);
	printf("\t grad %f\n", time_grad_int);

	BarrierMethod::print_timings();

	printf("Misc timings\n");
	printf("\t line search %f\n", time_line_search);
	printf("\t golden %f\n", time_golden);
}

static double _evaluate(
    void *instance,
    /*const*/ double *x,
    double *gg,
    /*const*/ int n,
    /*const*/ double step
    )
{
    return reinterpret_cast<JasonFull*>(instance)->evaluate(x, gg, n, step);
}

static void _is_movement(void *instance, double dis, double grad)
{
	reinterpret_cast<JasonFull*>(instance)->is_movement(dis, grad);
}

#include <fstream>
double JasonFull::evaluate(
    /*const*/ double *x,
    double *gg,
    /*const*/ int n,
    /*const*/ double step
    )
{
	//printf("debug out step %.10g\n", step);
	if(step == -2)
	{
		int i0, i1, i2;

		FILE* fptr = fopen("uvs.txt", "wt");
		fprintf(fptr, "{");
		for (int i = 0; i < tex_size; i++)
		{
			fprintf(fptr, "{%f,%f}", x[2 * i], x[2 * i + 1]);
			if (i != tex_size - 1)
			{
				fprintf(fptr, ",");
			}
		}
		fprintf(fptr, "}");
		fclose(fptr);

		fptr = fopen("isos.txt", "wt");
		FILE *ftris = fopen("tris.txt", "wt");
		fprintf(fptr, "{");
		fprintf(ftris, "{");
		for (int i = 0; i < this->ind_tex[0].size(); i++)
		{
			i0 = indchange[0][ind_tex[0][i][0]];
			i1 = indchange[0][ind_tex[0][i][1]];
			i2 = indchange[0][ind_tex[0][i][2]];

			fprintf(fptr, "{%f,%f,%f}", iso_tris[0][i][1][0], iso_tris[0][i][2][0], iso_tris[0][i][2][1]);
			if (i != this->ind_tex[0].size() - 1)
			{
				fprintf(fptr, ",");
			}
			fprintf(ftris, "{%d,%d,%d}", i0+1, i1+1, i2+1);
			if (i != this->ind_tex[0].size() - 1)
			{
				fprintf(ftris, ",");
			}
		}
		fprintf(fptr, "}");
		fprintf(ftris, "}");
		fclose(fptr);
		fclose(ftris);

		fstream out;
		out.setf(ios::fixed, ios::floatfield);
		out.precision(20);
		out.open("outtris.txt", fstream::out  );

		out << "{";
		double u0,u1,u2,v0,v1,v2, gx0,gy0,gx1,gy1,gx2,gy2;
		for(int i = 0; i < this->ind_tex[0].size();i++)
		{
			i0 = indchange[0][ind_tex[0][i][0]];
			i1 = indchange[0][ind_tex[0][i][1]];
			i2 = indchange[0][ind_tex[0][i][2]];
			
			/*u0 = i0 == -1 ? full_mesh[2* ind_tex[0][i][0]] : x[2*i0];
			v0 = i0 == -1 ? full_mesh[2* ind_tex[0][i][0]+1] : x[2*i0+1];

			u1 = i1 == -1 ? full_mesh[2* ind_tex[0][i][1]] : x[2*i1];
			v1 = i1 == -1 ? full_mesh[2* ind_tex[0][i][1]+1] : x[2*i1+1];

			u2 = i2 == -1 ? full_mesh[2* ind_tex[0][i][2]] : x[2*i2];
			v2 = i2 == -1 ? full_mesh[2* ind_tex[0][i][2]+1] : x[2*i2+1];

			gx0 = i0 == -1 ? 0 : gg[2*i0];
			gy0 = i0 == -1 ? 0 : gg[2*i0+1];

			gx1 = i1 == -1 ? 0 : gg[2*i1];
			gy1 = i1 == -1 ? 0 : gg[2*i1+1];

			gx2 = i2 == -1 ? 0 : gg[2*i2];
			gy2 = i2 == -1 ? 0 : gg[2*i2+1];*/

			u0 =  x[2*i0];
			v0 =  x[2*i0+1];

			u1 =  x[2*i1];
			v1 =  x[2*i1+1];

			u2 =  x[2*i2];
			v2 =  x[2*i2+1];
			/*
			gx0 =  gg[2*i0];
			gy0 = gg[2*i0+1];

			gx1 =  gg[2*i1];
			gy1 =  gg[2*i1+1];

			gx2 =  gg[2*i2];
			gy2 =  gg[2*i2+1];
			*/

			/*fx +=e_tri_func(iso_tris[0][ii][1][0],
							iso_tris[0][ii][2][0],
							iso_tris[0][ii][2][1],
							i0 == -1 ? full_mesh[2* ind_tex[0][ii][0]] : pos[2*i0],
							i0 == -1 ? full_mesh[2* ind_tex[0][ii][0] + 1] : pos[2*i0 + 1],
							i1 == -1 ? full_mesh[2* ind_tex[0][ii][1]] : pos[2*i1],
							i1 == -1 ? full_mesh[2* ind_tex[0][ii][1] + 1] :pos[2*i1+1],
							i2 == -1 ? full_mesh[2* ind_tex[0][ii][2]] : pos[2*i2],
							i2 == -1 ? full_mesh[2* ind_tex[0][ii][2] + 1] :pos[2*i2+1]);*/

			out <<std::fixed<< "{ {{" << u0 << ',' << v0 << "}, " 
				               << '{' << u1 << ',' << v1 << "}, " 
							   << '{' << u2 << ',' << v2 << "}}, " 
							   << iso_tris[0][i][1][0] << ", " << iso_tris[0][i][2][0] << ", " <<  iso_tris[0][i][2][1] << "}";
			if(i < this->ind_tex[0].size() - 1)out << ", ";
		}
		out << "}";

		out.close();
	}

    double fx = 0.0;
		
	if (step != -1 || step != -2)
	{
		fx = calc_error(x, 0, gg);
	}

	if (step == -1 || step == 0)
	{
		calc_full_gradient(x, gg);
	}
	
	//g.barrier_jason_gradient(x,gg);
	
    return fx;
}

static int _progress(
    void *instance,
    double *x,
    double *g,
    double fx,
    double xnorm,
    double gnorm,
    double step,
    int n,
    int k,
    int ls
    )
{
    return reinterpret_cast<JasonFull*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

int JasonFull::progress(
    double *x,
    double *g,
    double fx,
    double xnorm,
    double gnorm,
    double step,
    int n,
    int k,
    int ls
    )
{
    /*printf("barrier 3 Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f, %f, %f\n", fx, x[0], x[1], x[2], x[3]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");*/
    return 0;
}

static double _max_param(
     void *instance,
		double *f,
        double *x,
		double *gg,
        double *q,
		int n
    )
{
    return reinterpret_cast<JasonFull*>(instance)->max_param(f,x,gg,q,n);
}


#ifdef USE_LBFGS_LINESEARCH

double JasonFull::max_parameter(double *f,
    double *x,
	double *gg,
    double *q,
	int n)
{
	return max_param_interior(f,
    x,
	gg,
    q,
	n);
}

void JasonFull::is_movement(double dis, double grad)
{
#ifdef ZERO_MAX_TRIS
//	printf("dis = %.10g, grad mag = %.10g\n", dis, grad);
	if(grad > GRAD_MAG && dis < DIS_MAG )
	{
		//if magnitude big and no movement disable triangle
		zero_tris.push_back( max_param_tri );
	}
	else
	{
		zero_tris.clear();
	}
#endif
}


double cross2(vect2f a, vect2f b)
{
	return -1*b[0]*a[1] + a[0]*b[1];
}

double privSA[NUM_PROCS];
int privSAI[NUM_PROCS];
double JasonFull::max_param(
	double *f,
    double *x,
	double *gg,
    double *q,
	int n
	)
{
		double sA = 9999999;
		//printf("!\n");
#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif
	//INTERIOR MAX PARAM CHECK
//	printf("max param debug\n");
	

#ifdef ZERO_MAX_TRIS
	int priv_min_tri [ NUM_PROCS ];
	for ( int i = 0; i < NUM_PROCS; i++ )
	{
		priv_min_tri[i] = 0;
	}
#endif
	for (int i = 0; i < NUM_PROCS; i++)
	{
		privSA[i] = sA;
		privSAI[i] = 0;
		
	}

	//	printf("size in opt %d",ind_tex[0].size());
	double A;
	int index = 0;

#pragma omp parallel for //removed to get abc values in order
	for(int i = 0; i < ind_tex[0].size(); i++)
	{
		double u0,v0,u1,v1,u2,v2,gx0,gy0,gx1,gy1,gx2,gy2,a,b,c, disc, ta;
		int i0, i1, i2;
		i0 = indchange[0][ind_tex[0][i][0]];
		i1 = indchange[0][ind_tex[0][i][1]];
		i2 = indchange[0][ind_tex[0][i][2]];
		
		/*if(i == 499)
			printf("499 -> %d, %d, %d ....   %d, %d, %d\n", i0,i1,i2,ind_tex[0][i][0],ind_tex[0][i][1],ind_tex[0][i][2] );*/

	//	if(!(i0 == -1 && i1 == -1 && i2 == -1))
		{

		/*if(i0 == -1)
		{
			u0 = full_mesh[ ind_tex[0][i][0]*2 ];
			v0 = full_mesh[ ind_tex[0][i][0]*2 +1];

			gx0 = 0;
			gy0 = 0;
		}
		else*/
		{
			u0 = x[2*i0];
			v0 = x[2*i0+1];

			gx0 = q[i0 * 2];
			gy0 = q[i0 * 2 + 1];
		}

		/*if(i1 == -1)
		{
			u1 = full_mesh[ ind_tex[0][i][1]*2 ];
			v1 = full_mesh[ ind_tex[0][i][1]*2 +1];

			gx1 = 0;
			gy1 = 0;
		}
		else*/
		{
			u1 = x[2*i1];
			v1 = x[2*i1+1];

			gx1 = q[i1 * 2];
			gy1 = q[i1 * 2 + 1];
		}

		/*if(i2 == -1)
		{
			u2 = full_mesh[ ind_tex[0][i][2]*2 ];
			v2 = full_mesh[ ind_tex[0][i][2]*2 +1];

			gx2 = 0;
			gy2 = 0;
		}
		else*/
		{
			u2 = x[2*i2];
			v2 = x[2*i2+1];

			gx2 = q[i2 * 2];
			gy2 = q[i2 * 2 + 1];
		}
		//SWITCH GRADIENT WITH Q
		//Q is flipped
		/*if(271297)
			A = (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2));*/

		c = -(1e-20) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
		b = gy2 *(u0 - u1) + gy0 *(u1 - u2) + gy1 *(u2 -u0) + 
					gx2 *(v1 -v0) + gx1 *(v0 - v2) + gx0 *(v2 -v1);
		a = gx2* (-gy0 + gy1) + gx0* (-gy1 + gy2) + gx1* (gy0 - gy2);


		//Check disc
		disc = b*b - 4*a*c;

		if( disc < 0 )
		{
			//Do nothing
			//printf("neg disc\n");
		}
		else
		{
			if(a == 0)
			{
				ta = -c / b;

				/*if(i == 870 || i == 872 || i == 868)
				{
					printf("\nlin %d = %.10g\n", i, ta);
				}*/

				if(ta > 0 && ta < privSA [ omp_get_thread_num ( ) ])
				{
					/*if(i == 870 || i == 872 || i == 868)
				{
					printf("stop me\n");
				}*/
	#ifdef ZERO_MAX_TRIS			
					priv_min_tri[ omp_get_thread_num ( ) ] = i;
	#endif
					privSA [ omp_get_thread_num ( ) ] = ta;
					privSAI[omp_get_thread_num()] = i;
					//index =i;
				
				}
			}
			else
			{
			//Determine which case
			ta = (-1*b - sqrt( disc )) / (2*a);
			/*if(i == 870 || i == 872 || i == 868)
				{
					printf("\nq1 %d = %.18g\n", i, ta);
				}*/

			if(ta > 0 && ta < privSA [ omp_get_thread_num ( ) ] )
			{
				/*if(i == 870 || i == 872 || i == 868)
				{
					printf("stop me\n");
				}*/
#ifdef ZERO_MAX_TRIS
				priv_min_tri[ omp_get_thread_num ( ) ] = i;
#endif
				privSA [ omp_get_thread_num ( ) ] = ta;
				privSAI[omp_get_thread_num()] = i;
			//	index =i;
			}

			ta = (-1*b+ sqrt( disc )) / (2*a);
			/*if(i == 870 || i == 872 || i ==868)
				{
					printf("\nq2 %d = %.18g\n", i, ta);
				}*/
			if(ta > 0 && ta < privSA [ omp_get_thread_num ( ) ])
			{
				/*if(i == 870 || i == 872|| i == 868)
				{
					printf("stop me\n");
				}*/
#ifdef ZERO_MAX_TRIS			
				priv_min_tri[ omp_get_thread_num ( ) ] = i;
#endif
				privSA [ omp_get_thread_num ( ) ] = ta;
				privSAI[omp_get_thread_num()] = i;
			//	index =i;
				
			}

			}
		}
		}

		/*if(ind_tex[0][i][0] == 930 && ind_tex[0][i][1] == 928 && ind_tex[0][i][2] == 934)
		{
			printf("shit %d, %d, %d\n", ind_tex[0][i][0],ind_tex[0][i][1],ind_tex[0][i][2]);
			printf("ind change %d, %d, %d\n", i0, i1,i2);
			printf("a %.10g, b %.10g, c %.10g", a,b,c);
			printf("sa = %.10g, %.10g\n", (-1*b - sqrt( disc )) / (2*a), (-1*b + sqrt( disc )) / (2*a));
		}*/
	}

	sA = privSA [ 0 ];
	index = privSAI[0];
#ifdef ZERO_MAX_TRIS
	max_param_tri = priv_min_tri[0];
#endif
	for ( int i = 1; i < omp_get_max_threads ( ); i++ )
	{
#ifdef ZERO_MAX_TRIS
		if (privSA[i] < sA)
		{
			max_param_tri = priv_min_tri[i];
		}
#endif
		sA = min ( sA, privSA [ i ] );
		if (privSA[i] < sA)
		{
			index = privSAI[i];
		}
	}
	::g.max_param_ind = index;
	//printf("max param = %.15g index = %d\n",sA, index);
//	printf("%d,%d,%d\n", indchange[0][ind_tex[0][index][0]],indchange[0][ind_tex[0][index][0]],indchange[0][ind_tex[0][index][0]]);


#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_sA_int += endT - startT;
#endif
	

//////////////////////////////////////////////////////////////////////////////
	//Boundary Max Param Check
	currentsA = sA;

#ifdef TURN_OFF_BOUNDARY
#else

	if(g.do_boundary)
	{
#ifdef SEAMLESS_WITHOUT_BOUND
#else
	if(boundaryorder[0].size() > 0)
		sA = BarrierMethod::max_param(f,x,gg,q,n);
#endif
	}

#endif
	
	//printf("max param = %.18g\n",sA);

	//delete[] privSA;

	
	//for(int i = index; i < ind_tex[0].size(); i++)
	//{
	////	if((ind_tex[0][i][0] == 137135 && ind_tex[0][i][1] == 137073 &&ind_tex[0][i][2] == 137079) ||
	////		(ind_tex[0][i][1] == 137070 && ind_tex[0][i][2] == 137069 &&ind_tex[0][i][0] == 137061)||
	////		(ind_tex[0][i][2] == 137070 && ind_tex[0][i][0] == 137069 &&ind_tex[0][i][1] == 137061))
	//	{
	//	double u0,v0,u1,v1,u2,v2,gx0,gy0,gx1,gy1,gx2,gy2,a,b,c, disc, ta;
	//	int i0, i1, i2;
	//	i0 = indchange[0][ind_tex[0][i][0]];
	//	i1 = indchange[0][ind_tex[0][i][1]];
	//	i2 = indchange[0][ind_tex[0][i][2]];
	//	//printf("%d,%d,%d\n",i0,i1,i2);
	//	//printf("%d,%d,%d\n",ind_tex[0][i][0],ind_tex[0][i][1],ind_tex[0][i][2]);
	//	
	//	if(i0 == -1)
	//	{
	//		u0 = full_mesh[ ind_tex[0][i][0]*2 ];
	//		v0 = full_mesh[ ind_tex[0][i][0]*2 +1];

	//		gx0 = 0;
	//		gy0 = 0;
	//	}
	//	else
	//	{
	//		u0 = x[2*i0];
	//		v0 = x[2*i0+1];

	//		gx0 = q[i0 * 2];
	//		gy0 = q[i0 * 2 + 1];
	//	}

	//	if(i1 == -1)
	//	{
	//		u1 = full_mesh[ ind_tex[0][i][1]*2 ];
	//		v1 = full_mesh[ ind_tex[0][i][1]*2 +1];

	//		gx1 = 0;
	//		gy1 = 0;
	//	}
	//	else
	//	{
	//		u1 = x[2*i1];
	//		v1 = x[2*i1+1];

	//		gx1 = q[i1 * 2];
	//		gy1 = q[i1 * 2 + 1];
	//	}

	//	if(i2 == -1)
	//	{
	//		u2 = full_mesh[ ind_tex[0][i][2]*2 ];
	//		v2 = full_mesh[ ind_tex[0][i][2]*2 +1];

	//		gx2 = 0;
	//		gy2 = 0;
	//	}
	//	else
	//	{
	//		u2 = x[2*i2];
	//		v2 = x[2*i2+1];

	//		gx2 = q[i2 * 2];
	//		gy2 = q[i2 * 2 + 1];
	//	}
	//	double tu0,tu1,tu2,tv0,tv1,tv2;
	//	A = (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2));


	//	c = -(1e-15) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
	//	//c =  (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
	//	b = gy2 *(u0 - u1) + gy0 *(u1 - u2) + gy1 *(u2 -u0) + 
	//				gx2 *(v1 -v0) + gx1 *(v0 - v2) + gx0 *(v2 -v1);
	//	a = gx2* (-gy0 + gy1) + gx0* (-gy1 + gy2) + gx1* (gy0 - gy2);

	//	disc = b*b - 4*a*c;

	//	vect2f aa;
	//	aa[0] = u1-u0;
	//	aa[1] = v1-v0;
	//	vect2f bb;
	//	bb[0] = u2 - u0;
	//	bb[1] = v2-v0;
	//	double c1 = cross2(aa,bb);
	//	tu0 = u0;
	//	tu1 = u1;
	//	tu2 = u2;
	//	tv0 = v0; 
	//	tv1 = v1;
	//	tv2 = v2;

	//	u0 += gx0*sA;
	//	u1 += gx1*sA;
	//	u2 += gx2*sA;
	//	v0 += gy0*sA;
	//	v1 += gy1*sA;
	//	v2 += gy2*sA;

	//	aa[0] = u1-u0;
	//	aa[1] = v1-v0;
	//	bb[0] = u2 - u0;
	//	bb[1] = v2-v0;
	//	double c2 = cross2(aa,bb);
	//	//printf("grad = {{%.15g,%.15g},{%.15g, %.15g},{%.15g,%.15g}};\n",gx0,gy0,gx1,gy1,gx2,gy2);
	//	//printf("uv = {{%.15g,%.15g},{%.15g, %.15g},{%.15g,%.15g}};\n",tu0,tv0,tu1,tv1,tu2,tv2);
	//	

	//	//SWITCH GRADIENT WITH Q
	//	//Q is flipped
	//		//if((A < 0 && (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2)) > 0 ) ||  (A > 0 && (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2)) < 0 )   );
	//		if( (c1 < 0 && c2 > 0) || (c1 > 0 && c2 < 0) );
	//		{
	//			//printf("i = %d, sA = %.10g\n a = %.10g, b = %.10g, c = %.10g, disc = %.10g\n",i, sA,a,b,c,disc);
	//			//printf("before: c1 = %.10g\n", c1);
	//			//printf("uv = {{%.15g,%.15g},{%.15g, %.15g},{%.15g,%.15g}};\n",tu0,tv0,tu1,tv1,tu2,tv2);
	//			//printf("grad = {{%.15g,%.15g},{%.15g, %.15g},{%.15g,%.15g}};\n",gx0,gy0,gx1,gy1,gx2,gy2);
	//			//printf("gradsa = {{%.15g,%.15g},{%.15g, %.15g},{%.15g,%.15g}};\n",gx0*sA,gy0*sA,gx1*sA,gy1*sA,gx2*sA,gy2*sA);

	//			////printf("after: c2 = %.15g\n",c2);
	//			//printf("uv2 = {{%.15g,%.15g},{%.15g, %.15g},{%.15g,%.15g}};\n",u0,v0,u1,v1,u2,v2);

	//			int t;
	//			//scanf("%d",t);

	//		}
	//	break;
	//	}
	//}


	//printf("max param %.10g\n", sA);

	::g.remove_tri_index = index;
	return sA;

	}

double JasonFull::max_param_interior(
	double *f,
    double *x,
	double *gg,
    double *q,
	int n
	)
{
	printf("max param interior ?  debug\n");
		double sA = 9999999;
	double u0,v0,u1,v1,u2,v2,gx0,gy0,gx1,gy1,gx2,gy2,a,b,c, disc, ta;

#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif
	//INTERIOR MAX PARAM CHECK

	for (int i = 0; i < ind_tex[0].size(); i++)
	{
		int i0, i1, i2;
		i0 = indchange[0][ind_tex[0][i][0]];
		i1 = indchange[0][ind_tex[0][i][1]];
		i2 = indchange[0][ind_tex[0][i][2]];

		if(i0 == -1 && i1 == -1 && i2 == -1)
			continue;
		//float u0 = faceData[i].he->v->uvpos[0];
		//	float v0 = faceData[i].he->v->uvpos[1];

		//	float u1 = faceData[i].he->next->v->uvpos[0];
		//	float v1 = faceData[i].he->next->v->uvpos[1];

		//	float u2 = faceData[i].he->next->next->v->uvpos[0];
		//	float v2 = faceData[i].he->next->next->v->uvpos[1];

		//	float A = (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2));
		//	//if (area >= 0)
		//	if (A < 0)
		//	{
		//		printf(" flipped = %d\n", i);
		//	}
		/*if(i0 == -1)
		{
			u0 = full_mesh[ ind_tex[0][i][0]*2 ];
			v0 = full_mesh[ ind_tex[0][i][0]*2 +1];

			gx0 = 0;
			gy0 = 0;
		}
		else*/
		{
			u0 = x[2*i0];
			v0 = x[2*i0+1];

			gx0 = q[i0 * 2];
			gy0 = q[i0 * 2 + 1];
		}

		/*if(i1 == -1)
		{
			u1 = full_mesh[ ind_tex[0][i][1]*2 ];
			v1 = full_mesh[ ind_tex[0][i][1]*2 +1];

			gx1 = 0;
			gy1 = 0;
		}
		else*/
		{
			u1 = x[2*i1];
			v1 = x[2*i1+1];

			gx1 = q[i1 * 2];
			gy1 = q[i1 * 2 + 1];
		}

		/*if(i2 == -1)
		{
			u2 = full_mesh[ ind_tex[0][i][2]*2 ];
			v2 = full_mesh[ ind_tex[0][i][2]*2 +1];

			gx2 = 0;
			gy2 = 0;
		}
		else*/
		{
			u2 = x[2*i2];
			v2 = x[2*i2+1];

			gx2 = q[i2 * 2];
			gy2 = q[i2 * 2 + 1];
		}
		//u0 = x[2*i0];
		//v0 = x[2*i0+1];

		//u1 = x[2*i1];
		//v1 = x[2*i1+1];

		//u2 = x[2*i2];
		//v2 = x[2*i2+1];

		////SWITCH GRADIENT WITH Q
		////Q is flipped
		//gx0 = q[i0 * 2];
		//gy0 = q[i0 * 2 + 1];
		//
		//gx1 = q[i1 * 2];
		//gy1 = q[i1 * 2 + 1];

		//gx2 = q[i2 * 2];
		//gy2 = q[i2 * 2 + 1];

		c = -(1e-5) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
		b = gy2 *(u0 - u1) + gy0 *(u1 - u2) + gy1 *(u2 -u0) + 
					gx2 *(v1 -v0) + gx1 *(v0 - v2) + gx0 *(v2 -v1);
		a = gx2* (-gy0 + gy1) + gx0* (-gy1 + gy2) + gx1* (gy0 - gy2);

		//Check disc
		disc = b*b - 4*a*c;

		if( disc < 0 )
		{
			//Do nothing
			//printf("neg disc\n");
		}
		else
		{
			//Determine which case
			ta = (-1*b - sqrt( disc )) / (2*a);

			if(ta > 0 && ta < sA)
			{
				sA = ta;
			}

			ta = (-1*b+ sqrt( disc )) / (2*a);
			if(ta > 0 && ta < sA)
			{
				sA = ta;
			}
		}

	}	
	return sA;

	}
#else
double JasonFull::max_param(
	double *f,
    double *x,
	double *gg,
    double *q,
	int n
	)
{

	double sA = 9999999;
	double u0,v0,u1,v1,u2,v2,gx0,gy0,gx1,gy1,gx2,gy2,a,b,c, disc, ta;

#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif
	//INTERIOR MAX PARAM CHECK
	for(int i = 0; i < ind_tex[0].s; i++)
	{
		u0 = x[2*ind_tex[0][i][0]];
		v0 = x[2*ind_tex[0][i][0]+1];

		u1 = x[2*ind_tex[0][i][1]];
		v1 = x[2*ind_tex[0][i][1]+1];

		u2 = x[2*ind_tex[0][i][2]];
		v2 = x[2*ind_tex[0][i][2]+1];

		//SWITCH GRADIENT WITH Q
		//Q is flipped
		gx0 = q[ind_tex[0][i][0] * 2];
		gy0 = q[ind_tex[0][i][0] * 2 + 1];
		
		gx1 = q[ind_tex[0][i][1] * 2];
		gy1 = q[ind_tex[0][i][1] * 2 + 1];

		gx2 = q[ind_tex[0][i][2] * 2];
		gy2 = q[ind_tex[0][i][2] * 2 + 1];

		c = -(1e-10) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
		b = gy2 *(u0 - u1) + gy0 *(u1 - u2) + gy1 *(u2 -u0) + 
					gx2 *(v1 -v0) + gx1 *(v0 - v2) + gx0 *(v2 -v1);
		a = gx2* (-gy0 + gy1) + gx0* (-gy1 + gy2) + gx1* (gy0 - gy2);

		//Check disc
		disc = b*b - 4*a*c;

		if( disc < 0 )
		{
			//Do nothing
			//printf("neg disc\n");
		}
		else
		{
			//Determine which case
			ta = (-1*b - sqrt( disc )) / (2*a);

			if(ta > 0 && ta < sA)
			{
				sA = ta;
			}

			ta = (-1*b+ sqrt( disc )) / (2*a);
			if(ta > 0 && ta < sA)
			{
				sA = ta;
			}
		}
	}
#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_sA_int += endT - startT;
#endif
	

//////////////////////////////////////////////////////////////////////////////
	//Boundary Max Param Check
	currentsA = sA;

	sA = BarrierMethod::max_param(f,x,gg,q,n);
//	printf("sA = %.10g\n", sA);

//////////////////////////////////////////////////////////////////////////////
	double fx = 0;
	fx = calc_error(x, 0, x);

	//perform linesearch
	#ifdef DO_FULL_TIMINGS
		startT = get_time();
	#endif
	double x1 = 0;
	double x2 = sA * SA_COEFF;
	double x3 = sA;
	double y3 = 99999999;
	double y1 = fx;
	double y2 = 0;

	y2 = calc_error(x,-1*x2,q);
	//printf("%d -> y2 = %.10g\n", 0, y2);
	if(!_finite(y2) )
	{
		printf("SHIT\n");
		y2 = 99999999;
		//y2 = g.barrier_scott_unsymmetric_function(x,-1*x2,q);
	}

	int iters = 0;
	//Find Point 3
	while(iters < 100 && y2 > y1)
	{
		x2 = x2 * .5;
		y2 = 0;

		y2 = calc_error(x,-1*x2,q);
	//	printf("%d -> y2 = %.10g\n", iters+1, y2);
		if(!_finite(y2) )
		{
			printf("SHIT\n");
			y2 = 99999999;
			//y2 = g.barrier_scott_unsymmetric_function(x,-1*x2,q);
		}
		iters++;
	}

	#ifdef DO_FULL_TIMINGS
		endT = get_time();
		time_line_search += endT - startT;
	#endif

	if(iters == 100)
		printf("shit\n");

//	printf("sA2 = %.10g\n", x2);
	#ifdef DO_FULL_TIMINGS
	startT = get_time();
	#endif
	if(sA != 0)
	sA = golden_ratio(x1,y1,x2,y2,x3,y3,x,q);
//	printf("sA3 = %.10g\n", sA);
	#ifdef DO_FULL_TIMINGS
	endT = get_time();
	time_golden += endT- startT;
	#endif


	for(int i = 0; i < tex_size*2; i++)
	{
		x[i] = x[i] + sA*q[i];
	}

	*f = calc_error(x, 0, x);
	calc_full_gradient(x, gg);

	return sA;

}
#endif
void JasonFull::init_lbfgs(int n)
{
	int nn = n;
	 nn += 7;
    nn /= 8;
    nn *= 8;
   
	m_x = lbfgs_malloc(nn*2 ) ;
	for(int i = 0; i < nn; i++)
		m_x[i] = 0;

	lbfgs_init_mem(nn*2, MPARAM, 0);
}



void JasonFull::runHessian(int iters, double *pos, double *ans)
{
	omp_set_num_threads(NUM_PROCS);
	BarrierMethod::run(iters, pos, ans);

	double *quad = new double[2 * tex_size];
	double *nonlin = new double[2 * tex_size];
	Eigen::VectorXd grad, dir, nonlingrad, nonlindir;
	grad.resize(2 * tex_size);
	nonlingrad.resize(2 * tex_size);
	dir.resize(2 * tex_size);
	nonlindir.resize(2 * tex_size);
	
//	evaluate(pos, NULL, 2 * tex_size, -2);

	jason_precompute_quad_hessian();

	double lastFx = calc_error(pos, 0, NULL);
	double fx;

	for (int j = 0; j < 2 * tex_size; j++)
	{
		ans[j] = pos[j];
	}

	memset(nonlin, 0, sizeof(double)* 2 * tex_size);

	ofstream outputs("hessianoutput.txt");
	outputs << "{";
	double start = get_time();

	int i;
	for (i = 0; i < iters; i++)
	{
		calc_full_gradient(ans, grad.data());

		dir = hessianSolver.solve(grad);

		double step = .999 * max_param(NULL, ans, NULL, dir.data(), 2 * tex_size);
		int k;
		for (k = 0; k < 32; k++)
		{
			fx = calc_error(ans, step, dir.data());
			if (fx < lastFx)
			{
				for (int j = 0; j < 2 * tex_size; j++)
				{
					ans[j] -= step * dir.data()[j];
				}
				break;
			}
			step /= 2.0;
		}
		if (k >= 32)
		{
			printf("line search failed.\n");
			break;
		}

		/*
		Convergence test.
		The criterion is given by the following formula:
		|gg(x)| / \max(1, |x|) < \epsilon
		*/
/*		if (xnorm < 1.0) xnorm = 1.0;
		if (gnorm / xnorm <= g.opt_epsilon) {
			printf("success gnorm:%f xnorm:%f\n",gnorm, xnorm);
			break;
		}
	*/
		/*
		Test for stopping criterion.
		The criterion is given by the following formula:
		(f(past_x) - f(x)) / f(x) < \delta
		*/
		/* We don't test the stopping criterion while i < 3. */
		if (3 <= i) {
			double rate = (lastFx - fx) / fx;

			if ( fx <= 256917.078058 ) {
			//if (rate < 1e-6) {
				printf("stop %g\n",fx);
				break;
			}
		}
		
		/* Store the current value of the objective function. */
		lastFx = fx;
		char buf [ 512 ];

		sprintf ( buf, "%.9g", calc_error(ans, 0, NULL ) );

		outputs << "{";
		outputs << i + 1;
		outputs << ",";
		outputs << get_time() - start;
		outputs << ",";
//		outputs << calc_error(ans, 0, NULL);
		outputs << buf;
		outputs << "},";
	}

	outputs << "}";
	outputs.close();

	printf("took %d iterations!\n", i);

	printf("Hessian error: %f\n", calc_error(ans, 0, NULL));

	delete[] quad;
	delete[] nonlin;
}


void JasonFull::runKVF(int iters, double *pos, double *ans)
{
	omp_set_num_threads(NUM_PROCS);
	BarrierMethod::run(iters, pos, ans);

	double *quad = new double[2 * tex_size];
	double *nonlin = new double[2 * tex_size];
	Eigen::VectorXd dir, grad;
	dir.resize(2 * tex_size);
	grad.resize(2 * tex_size);
	
	double lastFx = calc_error(pos, 0, NULL);
	double fx;

	for (int j = 0; j < 2 * tex_size; j++)
	{
		ans[j] = pos[j];
	}
	memset(nonlin, 0, sizeof(double)* 2 * tex_size);

	//innerFlapOptimization();
	ofstream outputs("kvfoutput.txt");
	outputs << "{";

	double start = get_time();

	int i;
	for (i = 0; i < 200; i++)
	{
		calc_full_gradient(ans, grad.data());

		if ( i == 0 )
		{
			addressCounter = 0;
			computeKVF_First_Time ( ans );
			storeAddresses();
		}
		else
		{
			computeKVF ( ans );
		}
		dir = kvfSolver.solve(grad);
		//double scottTest = 0;
		//for ( int kk = 0; kk < 2 * tex_size; kk++ )
		//{
		//	scottTest += dir.data() [ kk ] * dir.data() [ kk ];
		//}
		//scottTest = sqrt ( scottTest );
		//for ( int kk = 0; kk < 2 * tex_size; kk++ )
		//{
		//	dir.data() [ kk ] /= scottTest;
		//}

		double step = .999 * max_param(NULL, ans, NULL, dir.data(), 2 * tex_size);
		printf ( "%g\n", step );
		int k;
		for (k = 0; k < 32; k++)
		{
			fx = calc_error(ans, step, dir.data());
			if (fx < lastFx)
			{
				for (int j = 0; j < 2 * tex_size; j++)
				{
					ans[j] -= step * dir.data()[j];
				}
				break;
			}
			step /= 2.0;
		}
		if (k >= 32)
		{
			printf("line search failed.\n");
			break;
		}

		//if(i == 27)
		//{
		//	ofstream outputs("horseData.txt");
		//	outputs << "{";

		//	for (int j = 0; j < ind_tex[0].size(); j++)
		//	{
		//		if(j == 14296 || j == 14393 || j == 14561)
		//		{
		//			outputs << "{";
		//			int i1, i2, i3;

		//			i1 = ind_tex[0][j][0];
		//			i2 = ind_tex[0][j][1];
		//			i3 = ind_tex[0][j][2];

		//			double i1x = ans [ i1 * 2 ];
		//			double i1y = ans[ i1 * 2 + 1 ];
		//			double i2x = ans[ i2 * 2 ];
		//			double i2y = ans[ i2 * 2 + 1 ];
		//			double i3x = ans[ i3 * 2 ];
		//			double i3y = ans[ i3 * 2 + 1 ];

		//			outputs << "{";

		//			outputs << "{"; //current values

		//			outputs << "{";
		//			outputs << i1x;
		//			outputs << ",";
		//			outputs << i1y;
		//			outputs << "},";
		//			outputs << "{";
		//			outputs << i2x;
		//			outputs << ",";
		//			outputs << i2y;
		//			outputs << "},";
		//			outputs << "{";
		//			outputs << i3x;
		//			outputs << ",";
		//			outputs << i3y;
		//			outputs << "}";

		//			outputs << "},";

		//			double d1x = dir.data() [ i1 * 2 ];
		//			double d1y = dir.data()[ i1 * 2 + 1 ];
		//			double d2x = dir.data()[ i2 * 2 ];
		//			double d2y = dir.data()[ i2 * 2 + 1 ];
		//			double d3x = dir.data()[ i3 * 2 ];
		//			double d3y = dir.data()[ i3 * 2 + 1 ];

		//			outputs << "{"; //direction values

		//			outputs << "{";
		//			outputs << d1x;
		//			outputs << ",";
		//			outputs << d1y;
		//			outputs << "},";
		//			outputs << "{";
		//			outputs << d2x;
		//			outputs << ",";
		//			outputs << d2y;
		//			outputs << "},";
		//			outputs << "{";
		//			outputs << d3x;
		//			outputs << ",";
		//			outputs << d3y;
		//			outputs << "}";

		//			outputs << "},";

		//			double u0,v0,u1,v1,u2,v2,gx0,gy0,gx1,gy1,gx2,gy2,a,b,c, disc, ta;
		//			int i0;
		//			i0 = indchange[0][ind_tex[0][i][0]];
		//			i1 = indchange[0][ind_tex[0][i][1]];
		//			i2 = indchange[0][ind_tex[0][i][2]];
		//			
		//			u0 = ans[2*i0];
		//			v0 = ans[2*i0+1];

		//			gx0 = dir.data()[i0 * 2];
		//			gy0 = dir.data()[i0 * 2 + 1];

		//			u1 = ans[2*i1];
		//			v1 = ans[2*i1+1];

		//			gx1 = dir.data()[i1 * 2];
		//			gy1 = dir.data()[i1 * 2 + 1];

		//			u2 = ans[2*i2];
		//			v2 = ans[2*i2+1];

		//			gx2 = dir.data()[i2 * 2];
		//			gy2 = dir.data()[i2 * 2 + 1];

		//			c = -(1e-20) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
		//			b = gy2 *(u0 - u1) + gy0 *(u1 - u2) + gy1 *(u2 -u0) + 
		//						gx2 *(v1 -v0) + gx1 *(v0 - v2) + gx0 *(v2 -v1);
		//			a = gx2* (-gy0 + gy1) + gx0* (-gy1 + gy2) + gx1* (gy0 - gy2);

		//			outputs << "{"; //max param

		//			outputs << "{";
		//			outputs << a;
		//			outputs << ",";
		//			outputs << b;
		//			outputs << ",";
		//			outputs << c;
		//			outputs << "}";

		//			outputs << "}";

		//			outputs << "},";
		//		}
		//	}
		//	outputs << "}";
		//	outputs.close();
		//}
		//
		for (int j = 0; j < ind_tex[0].size(); j++)
		{
			int i1, i2, i3;

			i1 = ind_tex[0][j][0];
			i2 = ind_tex[0][j][1];
			i3 = ind_tex[0][j][2];
			
			double x12 = ans [ i1 * 2 ] - ans [ i2 * 2 ];
			double y12 = ans[ i1 * 2 + 1 ] - ans [ i2 * 2 + 1 ];
			double x23 = ans[ i2 * 2 ] - ans [ i3 * 2 ];
			double y23 = ans[ i2 * 2 + 1 ] - ans [ i3 * 2 + 1 ];
			double x31 = ans[ i3 * 2 ] - ans[ i1 * 2 ];
			double y31 = ans[ i3 * 2 + 1 ] - ans[ i1 * 2 + 1 ];

			double area = ( 0.5 *  ( y12 * x31 -x12 * y31 ) ); // really the sqrt of area

			if(area < 0)
			{
				printf("Triangle %d: %d, %d, %d on iteration %d\n", j, i1, i2, i3, i);
			}
		}

		/*
		Convergence test.
		The criterion is given by the following formula:
		|gg(x)| / \max(1, |x|) < \epsilon
		*/
/*		if (xnorm < 1.0) xnorm = 1.0;
		if (gnorm / xnorm <= g.opt_epsilon) {
			printf("success gnorm:%f xnorm:%f\n",gnorm, xnorm);
			break;
		}
	*/
		/*
		Test for stopping criterion.
		The criterion is given by the following formula:
		(f(past_x) - f(x)) / f(x) < \delta
		*/
		/* We don't test the stopping criterion while i < 3. */
 		if (3 <= i) {
			double rate = (lastFx - fx) / fx;

			//if ( fx <= 256917.078058 ) {
			if (rate < 1e-200) {
				printf("stop %g\n",fx);
				break;
			}
		}
		
		/* Store the current value of the objective function. */
		lastFx = fx;
/*
		char buf [ 512 ];
		sprintf ( buf, "%.9g", calc_error(ans, 0, NULL ) );

		outputs << "{";
		outputs << i + 1;
		outputs << ",";
		outputs << get_time() - start;
		outputs << ",";
		outputs << buf;
		outputs << "},"; */
	}

	outputs << "}";
	outputs.close();

	printf("took %d iterations!\n", i);

	printf("KVF error: %f\n", calc_error(ans, 0, NULL));

	delete[] quad;
	delete[] nonlin;
}

int JasonFull::run(int iters, double *pos, double *ans)
{
	//ind_tex[0] = indices_tex;
	omp_set_num_threads(NUM_PROCS);
	BarrierMethod::run(iters, pos, ans);
	//cout << "Starting Jason Full\n";
	//omp_set_num_threads(1);

	double fx;
	printf("run lbfgs\n");
	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
	
	param.delta =1e-3;//1e-6;
	
	/*param.m = MPARAM;
	param.past = 5;
*/
	param.m = 5;
	param.past = 3;
	
	param.gtol = 1e-3;
	//param.epsilon = 1e-4;
	//printf("Convergence Epsilon = %d\n", param.epsilon);
	/*switch( g.opt_epsilon ){
	}*/
	param.epsilon = g.opt_epsilon ;
	printf("Convergence Epsilon = %.10g\n", param.epsilon);
	param.ftol = 1e-4;
	param.max_iterations = max_iters_t;
	param.wolfe = .9;
	param.min_step = 1e-30;
	param.max_linesearch = 100;

	//Big
	/*param.delta =1e-15;
	param.past = 0;
	param.gtol = 1e-15;
	param.epsilon = 1e-3;
	param.ftol = 1e-4;
	param.max_iterations = max_iters_t;
	param.wolfe = .9;
	param.min_step = 1e-30;
	param.max_linesearch = 100;
*/
	//param.
	
	
	for (int i = 0;i < 2*tex_size; i++) 
	{
		m_x[i] = pos[i];
    }
	//
	//Precompute
	//
	


#ifdef USE_LBFGS_LINESEARCH
	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
	printf("debug run start\n");
#ifdef ZERO_MAX_TRIS
	int ret = lbfgs3( tex_size*2 , m_x, &fx, _evaluate, _progress, _max_param,  _is_movement, this, &param);
#else
	int ret = lbfgs3( tex_size*2 , m_x, &fx, _evaluate, _progress, _max_param, this, &param);
#endif
#else
	int ret = lbfgs2( tex_size*2 , m_x, &fx, _evaluate, _progress, _max_param, this, &param);
#endif

    /* Report the result. */
	printf("debug end\n");
	if(::g.printInfo)
	{
		printf("L-BFGS optimization terminated with status code = %d\n", ret);
		printf("  fx = %f\n", fx);
	}

	for (int i = 0;i < 2*tex_size; i++) 
	{
		ans[i] = m_x[i];
    }

	//lbfgs_free ( m_x );
	return ret;
}

double JasonFull::barrier_jason_function_tri(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2)
{
	//return (1.0/(((L)))*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*(((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx))+(((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*(((xx)*(xx))+((y)*(y))))));
/*	
	double den = 2 * L * y;
	double a = (-(u0*y) + u1*y - L*v0 + xx*v0 - xx*v1 + L*v2)/den;
	double b = (L*u0 - xx*u0 + xx*u1 - L*u2 - y*v0 + y*v1)/den;
	double c = (-(u0*y) + u1*y + L*v0 - xx*v0 + xx*v1 - L*v2)/den;
	double d = (xx*(u0 - u1) + L*(-u0 + u2) + y*(-v0 + v1))/den;
	
	return L*y*(2*(pow(a,2) + pow(b,2) + pow(c,2) + pow(d,2))*(1 + pow(pow(a,2) + pow(b,2) - pow(c,2) - pow(d,2),2)))/
   pow(pow(a,2) + pow(b,2) - pow(c,2) - pow(d,2),2);
   */
	double area = u1*v0 - u2*v0 - u0*v1 + u2*v1 + u0*v2 - u1*v2;
	double u01 = u0 - u1;
	double v01 = v0 - v1;
	double vecU = L*(u0 - u2) - u01*xx;
	double vecV = L*(v0 - v2) - v01*xx;
	return ( area * area + L * L * y * y ) * ( vecU * vecU + vecV * vecV + (u01 * u01 + v01 * v01) * y * y ) / ( L * y * area * area );

	
	//	return ( area * area + L * L * y * y ) * ( vecU * vecU + vecV * vecV + ((u0 - u1)*(u0 - u1) + (v0 - v1) * (v0 - v1)) * y * y ) / ( L * y * area * area );
}

double JasonFull::barrier_jason_function(double *pos)
{
	double fx = 0;

#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif
#pragma omp parallel for reduction(+:fx)
	for (int ii = 0; ii < ind_tex[0].size(); ii++)
	{

		int i0, i1, i2;
		i0 = indchange[0][ind_tex[0][ii][0]];
		i1 = indchange[0][ind_tex[0][ii][1]];
		i2 = indchange[0][ind_tex[0][ii][2]];

		if(i0 == -1 && i1 == -1 && i2 == -1)
			continue;

#ifdef JASON_PRE_COMP

		//printf("%f,%f,%f\n",iso_tris[0][ii][1][0],iso_tris[0][ii][2][0],iso_tris[0][ii][2][1]);
		fx += paper_jason_function_tri( &this->precomp[ii], 
			 pos[2*i0],
			pos[2*i0 + 1],
			pos[2*i1],
			pos[2*i1+1],
			pos[2*i2],
			pos[2*i2+1]
			);

	
#else
		fx +=e_tri_func(iso_tris[0][ii][1][0],
						iso_tris[0][ii][2][0],
						iso_tris[0][ii][2][1],
						i0 == -1 ? full_mesh[2* ind_tex[0][ii][0]] : pos[2*i0],
						i0 == -1 ? full_mesh[2* ind_tex[0][ii][0] + 1] : pos[2*i0 + 1],
						i1 == -1 ? full_mesh[2* ind_tex[0][ii][1]] : pos[2*i1],
						i1 == -1 ? full_mesh[2* ind_tex[0][ii][1] + 1] :pos[2*i1+1],
						i2 == -1 ? full_mesh[2* ind_tex[0][ii][2]] : pos[2*i2],
						i2 == -1 ? full_mesh[2* ind_tex[0][ii][2] + 1] :pos[2*i2+1]);
		//if(ii == 0)
			//	printf("first func = %f\n", fx);
#endif
		/*if(!_finite(e_tri_func(iso_tris[0][ii][1][0],
						iso_tris[0][ii][2][0],
						iso_tris[0][ii][2][1],
						i0 == -1 ? full_mesh[2* ind_tex[0][ii][0]] : pos[2*i0],
						i0 == -1 ? full_mesh[2* ind_tex[0][ii][0] + 1] : pos[2*i0 + 1],
						i1 == -1 ? full_mesh[2* ind_tex[0][ii][1]] : pos[2*i1],
						i1 == -1 ? full_mesh[2* ind_tex[0][ii][1] + 1] :pos[2*i1+1],
						i2 == -1 ? full_mesh[2* ind_tex[0][ii][2]] : pos[2*i2],
						i2 == -1 ? full_mesh[2* ind_tex[0][ii][2] + 1] :pos[2*i2+1])))
		{
			printf("%d, shit %d, %d, %d\n", ii , i0, i1,i2);
			printf("%.15g, %.15g, %.15g, %.15g, %.15g, %.15g\n",i0 == -1 ? full_mesh[2* ind_tex[0][ii][0]] : pos[2*i0],
						i0 == -1 ? full_mesh[2* ind_tex[0][ii][0] + 1] : pos[2*i0 + 1],
						i1 == -1 ? full_mesh[2* ind_tex[0][ii][1]] : pos[2*i1],
						i1 == -1 ? full_mesh[2* ind_tex[0][ii][1] + 1] :pos[2*i1+1],
						i2 == -1 ? full_mesh[2* ind_tex[0][ii][2]] : pos[2*i2],
						i2 == -1 ? full_mesh[2* ind_tex[0][ii][2] + 1] :pos[2*i2+1]);
		}*/
	}
	//int tttt;
//	scanf("%d", &tttt);
#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_func_int += endT - startT;
#endif
	return fx;
}

double JasonFull::calc_error(double* pos, double step, double* grad)
{
	double fx = 0;
	if(step == 0)
	{
		//Calc Interior
		fx = barrier_jason_function(pos);
		//printf("interior %.10g\n", fx);
#ifdef TURN_OFF_BOUNDARY
#else
		if(g.do_boundary)
		{
		//Perform Global
		fx += BarrierMethod::calc_error(pos, step, grad);
		}
#endif
		//printf("total %.10g\n", fx);
	}
	else
	{
		for(int i = 0; i < tex_size*2; i++)
			temp_pos[i] = pos[i] - step*grad[i];

		//Calc Interior
		fx = barrier_jason_function(temp_pos);
#ifdef TURN_OFF_BOUNDARY
#else
		if(g.do_boundary)
		{
		//Perform Global
		fx += BarrierMethod::calc_error(temp_pos, step, grad);
		}
#endif
	}

	
	//printf("Jason Error\n");
	return fx;
}
void JasonFull::calc_full_gradient(double *pos,  double *gg )
{
	calc_interior_gradient(pos, gg);
//	printf("%.10g, %.10g, %.10g, %.10g, %.10g, %.10g\n\n",g[0],g[1],g[2],g[3],g[4],g[5] );
#ifdef TURN_OFF_BOUNDARY
#else

	if(g.do_boundary)
	BarrierMethod::calc_full_gradient(pos, gg);

//	printf("%.10g, %.10g, %.10g, %.10g, %.10g, %.10g\n\n",g[0],g[1],g[2],g[3],g[4],g[5] );
#endif
//	exit(1);
}

void JasonFull::jason_gradient_quad_nonlin(double *pos, double *quad, double *nonlin)
{
#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif

	int i, j;
#pragma omp parallel private ( i, j )
	{
		memset(threadGrad[omp_get_thread_num()], 0, sizeof(double)*tex_size * 2 * 2);
		//for (i = 0; i < tex_size * 2; i++)
		//   {
		//          threadGrad [ omp_get_thread_num ( ) ] [ i ] = 0;
		//   }

#pragma omp for
		for (i = 0; i < ind_tex[0].size(); i++)
		{
			int i0, i1, i2;

			i0 = ind_tex[0][i][0];
			i1 = ind_tex[0][i][1];
			i2 = ind_tex[0][i][2];

			if (i0 == -1 && i1 == -1 && i2 == -1)
				continue;

			paper_jason_function_grad_quad_nonlin(threadGrad[omp_get_thread_num()], &(threadGrad[omp_get_thread_num()][tex_size * 2]), &this->precomp[i],
				pos[2 * i0],
				pos[2 * i0 + 1],
				pos[2 * i1],
				pos[2 * i1 + 1],
				pos[2 * i2],
				pos[2 * i2 + 1],
				i0, i1, i2
				);
		}

#pragma omp for
		for (i = 0; i < tex_size * 2; i++)
		{
			quad[i] = threadGrad[0][i];
			nonlin[i] = threadGrad[0][i + tex_size * 2];

			for (j = 1; j < NUM_PROCS; j++)
			{
				quad[i] += threadGrad[j][i];
				nonlin[i] += threadGrad[j][i + tex_size * 2];
			}
		}
	}


#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_grad_int += endT - startT;
#endif
}

void JasonFull::calc_quad_nonlin_gradient(double *pos, double *quad, double *nonlin)
{
	jason_gradient_quad_nonlin(pos, quad, nonlin);
#ifdef TURN_OFF_BOUNDARY
#else

	if (g.do_boundary)
		BarrierMethod::calc_full_gradient(pos, nonlin);
#endif
}


void JasonFull::jason_precompute_quad_hessian(void)
{
	double start = get_time();

	vector<Eigen::Triplet<double> > trips;
	trips.reserve(ind_tex[0].size() * 18+2*tex_size);

	for (int i = 0; i < ind_tex[0].size(); i++)
	{
		int i0, i1, i2;

		i0 = ind_tex[0][i][0];
		i1 = ind_tex[0][i][1];
		i2 = ind_tex[0][i][2];

		if (i0 == -1 && i1 == -1 && i2 == -1)
			continue;

		paper_jason_function_hessian_quad( trips, &(precomp[i]), i0, i1, i2 );
	}

	for (int i = 0; i < 2 * tex_size; i++)
	{
		trips.push_back(Eigen::Triplet<double>(i, i, 0.0001));
	}

	hessian.resize(2 * tex_size, 2 * tex_size);
	hessian.setFromTriplets(trips.begin(), trips.end());

/*	FILE *fptr = fopen("hessian.txt", "wt");
	fprintf(fptr, "{");
	for (int k = 0; k < hessian.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(hessian, k); it; ++it)
		{
			fprintf(fptr,"{%d,%d,%f},", it.row(), k, it.value());
		}
	}
	fprintf(fptr, "}");
	fclose(fptr);*/
	double startT = get_time();

	hessianSolver.analyzePattern(hessian);
	double analyzeT = get_time();
	hessianSolver.factorize(hessian);
	double end = get_time();
	printf("Took %f to compute hessian, %f to analyze, and %f to factorization with %f total\n", startT - start, analyzeT - startT, end - analyzeT, end-start);
}

#ifdef PARALLEL_GRAD
void JasonFull::barrier_jason_gradient(double *pos, double *newgradient)
{
#ifdef DO_FULL_TIMINGS
       double startT = get_time();
#endif
/*
       for(int i = 0; i < tex_size*2; i++)
       {
              newgradient[i] = 0;
       }
 
       for(int i = 0; i < ind_tex[0].s; i++)
       {
              e_gradient(newgradient, iso_tris[0][i][1][0], iso_tris[0][i][2][0], iso_tris[0][i][2][1],
                     pos[ind_tex[0][i][0]*2], pos[ind_tex[0][i][0]*2+1],
                     pos[ind_tex[0][i][1]*2], pos[ind_tex[0][i][1]*2+1],
                     pos[ind_tex[0][i][2]*2], pos[ind_tex[0][i][2]*2+1],
                     ind_tex[0][i][0],ind_tex[0][i][1],ind_tex[0][i][2]
                     );
 
       }
*/

       int i, j;
#pragma omp parallel private ( i, j )
{
		memset(threadGrad[omp_get_thread_num()], 0, sizeof(double)*tex_size * 2);
	   //for (i = 0; i < tex_size * 2; i++)
    //   {
    //          threadGrad [ omp_get_thread_num ( ) ] [ i ] = 0;
    //   }
 
#pragma omp for
	   for (i = 0; i < ind_tex[0].size(); i++)
       {
		   int i0, i1, i2;
		   /*i0 = indchange[0][ind_tex[0][i][0]];
		   i1 = indchange[0][ind_tex[0][i][1]];
		   i2 = indchange[0][ind_tex[0][i][2]];
*/
		   i0 = ind_tex[0][i][0];
		   i1 = ind_tex[0][i][1];
		   i2 = ind_tex[0][i][2];

		   if(i0 == -1 && i1 == -1 && i2 == -1)
			continue;
              /*e_gradient(threadGrad [ omp_get_thread_num ( ) ], iso_tris[0][i][1][0], iso_tris[0][i][2][0], iso_tris[0][i][2][1],
                     pos[i0*2], 
					 pos[i0*2+1],
                     pos[i1*2], 
					 pos[i1*2+1],
                     pos[i2*2], 
					 pos[i2*2+1],
                     i0,i1,i2
                     );
*/
#ifdef JASON_PRE_COMP
		   /*
		   paper_jason_function_grad(threadGrad [ omp_get_thread_num ( ) ], &this->precomp[i],
                         pos[2*i0],
						 pos[2*i0 + 1],
						 pos[2*i1],
						 pos[2*i1+1],
						 pos[2*i2],
						 pos[2*i2+1],
                     i0,i1,i2
                     );
*/
		   jason_function_grad(threadGrad[omp_get_thread_num()], this->precomp[i].L, this->precomp[i].x, this->precomp[i].y,
			   pos[2*i0],
			   pos[2*i0 + 1],
			   pos[2*i1],
			   pos[2*i1+1],
			   pos[2*i2],
			   pos[2*i2+1],
			   i0,i1,i2
			   );
#else
			  e_gradient(threadGrad [ omp_get_thread_num ( ) ], iso_tris[0][i][1][0], iso_tris[0][i][2][0], iso_tris[0][i][2][1],
                     i0 == -1 ? full_mesh[2* ind_tex[0][i][0]] : pos[2*i0],
						i0 == -1 ? full_mesh[2* ind_tex[0][i][0] + 1] : pos[2*i0 + 1],
						i1 == -1 ? full_mesh[2* ind_tex[0][i][1]] : pos[2*i1],
						i1 == -1 ? full_mesh[2* ind_tex[0][i][1] + 1] :pos[2*i1+1],
						i2 == -1 ? full_mesh[2* ind_tex[0][i][2]] : pos[2*i2],
						i2 == -1 ? full_mesh[2* ind_tex[0][i][2] + 1] :pos[2*i2+1],
                     i0,i1,i2
                     );
#endif
       }
 
#pragma omp for
	//   printf("max threads %d\n", omp_get_max_threads());
       for ( i = 0; i < tex_size*2; i++ )
       {
              newgradient [ i ] = threadGrad [ 0 ] [ i ];
			  
			 /* if (i == 0)
			  {
				  printf("thread %d = %f\n", 0, threadGrad[0][i]);
			  }*/
			  for (j = 1; j < NUM_PROCS; j++)
              {
				 /* if (i == 0)
				  {
					  printf("thread %d = %f\n", j, threadGrad[j][i]);
				  }*/
                     newgradient [ i ] += threadGrad [ j ] [ i ];
              }
			  //printf("%g\n", newgradient[i]);
       }
	  // printf ( "test\n" );
}

#ifdef ZERO_MAX_TRIS
	if(zero_tris.size() > 0)
		printf("cut triangles: ");
	for(int i  = 0; i < zero_tris.size(); i++)
	{
		printf("%d, ", zero_tris[i]);	

		if( indchange[0][ind_tex[0][zero_tris[i]][0]] != -1 )
		{
			newgradient[ind_tex[0][zero_tris[i]][0]*2] = 0;
			newgradient[ind_tex[0][zero_tris[i]][0]*2 + 1] = 0;
		}

		if( indchange[0][ind_tex[0][zero_tris[i]][1]] != -1 )
		{
			newgradient[ind_tex[0][zero_tris[i]][1]*2] = 0;
			newgradient[ind_tex[0][zero_tris[i]][1]*2 + 1] = 0;
		}

		if( indchange[0][ind_tex[0][zero_tris[i]][2]] != -1 )
		{
			newgradient[ind_tex[0][zero_tris[i]][2]*2] = 0;
			newgradient[ind_tex[0][zero_tris[i]][2]*2 + 1] = 0;
		}
	}
	if(zero_tris.size() > 0)
		printf("\n");
#endif

#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_grad_int += endT - startT;
#endif
}
#else
void JasonFull::barrier_jason_gradient(double *pos, double *newgradient)
{
#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif
	for(int i = 0; i < tex_size*2; i++)
	{
		newgradient[i] = 0;
	}

	for(int i = 0; i < ind_tex[0].s; i++)
	{
		e_gradient(newgradient, iso_tris[0][i][1][0], iso_tris[0][i][2][0], iso_tris[0][i][2][1], 
			pos[ind_tex[0][i][0]*2], pos[ind_tex[0][i][0]*2+1], 
			pos[ind_tex[0][i][1]*2], pos[ind_tex[0][i][1]*2+1], 
			pos[ind_tex[0][i][2]*2], pos[ind_tex[0][i][2]*2+1],
			ind_tex[0][i][0],ind_tex[0][i][1],ind_tex[0][i][2]
			);
/*
		double L = ;
		double xx = ;
		double y = ;

		double u0 = ;
		double v0 = ;

		double u1 = ;
		double v1 = ;

		double u2 = ;
		double v2 = ;*/

		/*newgradient[ind_tex[0][i][0]*2] += ((2.0)*1.0/(((L)))*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*(((L)*(((2.0)*(u0))+((-1.0)*(u1))+((-1.0)*(u2)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((((u0)+((-1.0)*(u2)))*((u1)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))*((v0)+((-1.0)*(v2)))*((y)*(y)))+(((L)*(L)*(L))*(((-1.0)*((v0)+((-1.0)*(v1)))*(((u2)*(u2))+((2.0)*((v0)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2))))))+(((u1)*(u1))*((v0)+((-1.0)*(v2))))+((-1.0)*(u0)*((u1)+((-1.0)*(u2)))*(((2.0)*(v0))+((-1.0)*(v1))+((-1.0)*(v2))))+((u1)*(u2)*(((-1.0)*(v1))+(v2))))*(xx)*((y)*(y)))+((-1.0)*((u0)+((-1.0)*(u1)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*(((-1.0)*((u0)+((-1.0)*(u2)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))))+((((-1.0)*(v0))+(v1))*(((-1.0)*((u0)+((-1.0)*(u1)))*((u1)+((-1.0)*(u2))))+((-1.0)*((v0)+((-1.0)*(v1)))*((v1)+((-1.0)*(v2)))))*((xx)*(xx))*((y)*(y)))+((((-1.0)*(v0))+(v1))*(((-1.0)*((u0)+((-1.0)*(u1)))*((u1)+((-1.0)*(u2))))+((-1.0)*((v0)+((-1.0)*(v1)))*((v1)+((-1.0)*(v2)))))*((y)*(y)*(y)*(y)))))));
		newgradient[ind_tex[0][i][0]*2+1] += (1.0/(((L)))*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*(((2.0)*((u1)+((-1.0)*(u2)))*((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))*((((L)*(L))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx))+(((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((u1)+((-1.0)*(u2)))*(((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx))+(((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*(((xx)*(xx))+((y)*(y))))))+((2.0)*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((v0)+((-1.0)*(v2))))+((L)*(((-2.0)*(v0))+(v1)+(v2))*(xx))+(((v0)+((-1.0)*(v1)))*(((xx)*(xx))+((y)*(y))))))));

		newgradient[ind_tex[0][i][1]*2] +=(1.0/(((L)))*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*(((-2.0)*(L)*((u0)+((-1.0)*(u2)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2)))))*(xx))+((-2.0)*((L)*(L)*(L)*(L))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2)))))*((v0)+((-1.0)*(v2)))*((y)*(y)))+((2.0)*((L)*(L)*(L))*((((v0)+((-1.0)*(v1)))*(((u2)*(u2))+((2.0)*(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2)))))))+((u1)*(u2)*((v0)+((-1.0)*(v2))))+(((u0)*(u0))*(((2.0)*(v0))+((-1.0)*(v1))+((-1.0)*(v2))))+((u0)*(((u1)*(((-1.0)*(v0))+(v2)))+((u2)*(((-3.0)*(v0))+((2.0)*(v1))+(v2))))))*(xx)*((y)*(y)))+((2.0)*((u0)+((-1.0)*(u1)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+((-2.0)*((L)*(L))*((v0)+((-1.0)*(v1)))*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))));
		newgradient[ind_tex[0][i][1]*2+1] +=((2.0)*1.0/(((L)))*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*(((L)*((v0)+((-1.0)*(v2)))*((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((u0)+((-1.0)*(u2)))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2)))))*((y)*(y)))+(((L)*(L)*(L))*(((-2.0)*((u0)*(u0)*(u0)))+((2.0)*((u0)*(u0))*((u1)+((2.0)*(u2))))+((u1)*(((2.0)*((u2)*(u2)))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2))))))+((-1.0)*(u0)*(((2.0)*(u2)*(((2.0)*(u1))+(u2)))+(((v0)+((-1.0)*(v2)))*(((2.0)*(v0))+((-1.0)*(v1))+((-1.0)*(v2))))))+((u2)*((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx)*((y)*(y)))+(((v0)+((-1.0)*(v1)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((u0)+((-1.0)*(u1)))*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))));

		newgradient[ind_tex[0][i][2]*2] +=(1.0/(((L)))*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*(((-2.0)*(L)*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((L)*((u0)+((-1.0)*(u2))))+((((-1.0)*(u0))+(u1))*(xx)))*(((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))+(((L)*(L))*((y)*(y)))))+((2.0)*(((-1.0)*(v0))+(v1))*((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))*((((L)*(L))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx))+(((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*(((-1.0)*(v0))+(v1))*(((((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2)))*(((u1)*(v0))+((-1.0)*(u2)*(v0))+((-1.0)*(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+((-1.0)*(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+((-1.0)*(u2)))*((u0)+((-1.0)*(u2))))+(((v0)+((-1.0)*(v2)))*((v0)+((-1.0)*(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx))+(((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*(((xx)*(xx))+((y)*(y))))))));
		newgradient[ind_tex[0][i][2]*2+1] += ((2.0)*1.0/(((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))))*1.0/(((y)))*((((v0)+((-1.0)*(v1)))*((((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2))))*(((u2)*(((-1.0)*(v0))+(v1)))+((u1)*((v0)+((-1.0)*(v2))))+((u0)*(((-1.0)*(v1))+(v2)))))*(xx))+((-1.0)*((L)*(L)*(L))*((u0)+((-1.0)*(u2)))*(((u0)*(u0))+((u1)*(u2))+((-1.0)*(u0)*((u1)+(u2)))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*((y)*(y)))+(((L)*(L))*(((2.0)*((u0)*(u0)*(u0)))+((-2.0)*((u1)*(u1))*(u2))+((-2.0)*((u0)*(u0))*(((2.0)*(u1))+(u2)))+((-1.0)*(u2)*(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))+((u0)*(((2.0)*(u1)*((u1)+((2.0)*(u2))))+(((v0)+((-1.0)*(v1)))*(((2.0)*(v0))+((-1.0)*(v1))+((-1.0)*(v2))))))+((-1.0)*(u1)*((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v2)))))*(xx)*((y)*(y)))+((L)*((((v0)+((-1.0)*(v2)))*((((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))*(((u2)*((v0)+((-1.0)*(v1))))+((u0)*((v1)+((-1.0)*(v2))))+((u1)*(((-1.0)*(v0))+(v2))))))+((((-1.0)*(u0))+(u1))*((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*((xx)*(xx))*((y)*(y)))+((((-1.0)*(u0))+(u1))*((((u0)+((-1.0)*(u1)))*((u0)+((-1.0)*(u1))))+(((v0)+((-1.0)*(v1)))*((v0)+((-1.0)*(v1)))))*((y)*(y)*(y)*(y)))))));
*/

		//double L = iso_tris[0][i][1][0];
		//double x = iso_tris[0][i][2][0];
		//double y = iso_tris[0][i][2][1];

		//double x1 = pos[ind_tex[0][i][0]*2];
		//double y1 = pos[ind_tex[0][i][0]*2+1];

		//double x2 = pos[ind_tex[0][i][1]*2];
		//double y2 = pos[ind_tex[0][i][1]*2+1];

		//double x3 = pos[ind_tex[0][i][2]*2];
		//double y3 = pos[ind_tex[0][i][2]*2+1];

		//newgradient[ind_tex[0][i][0]*2] += (2*((pow(x,2) + pow(y,2))*(pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*(-y2 + y3)*
  //    pow(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3,2) - 
  //   2*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*(-y2 + y3)*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)) + 
  //   2*(pow(L,2)*(x1 - x3) + L*x*(-2*x1 + x2 + x3) + 
  //      (x1 - x2)*(pow(x,2) + pow(y,2)))*
  //    (x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)))/
  // (L*y*pow(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3),3));

		//newgradient[ind_tex[0][i][0]*2+1] += (2*(x2 - x3)*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*
  //    pow(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3,2) - 
  //   2*(x2 - x3)*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)) + 
  //   2*(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))*
  //    ((pow(x,2) + pow(y,2))*(y1 - y2) + pow(L,2)*(y1 - y3) + 
  //      L*x*(-2*y1 + y2 + y3))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)))/
  // (L*y*pow(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3),3));
		//	newgradient[ind_tex[0][i][1]*2] +=(-2*pow(L,2)*pow(y,2)*(pow(x,2) + pow(y,2))*(y1 - y2)*
  //    (pow(x1,2) + x2*x3 - x1*(x2 + x3) + (y1 - y2)*(y1 - y3)) - 
  //   2*pow(L,4)*pow(y,2)*(pow(x1 - x3,2) + pow(y1 - y3,2))*
  //    (y1 - y3) - 2*L*x*(x1 - x3)*
  //    pow(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3),3) + 
  //   2*(x1 - x2)*(pow(x,2) + pow(y,2))*
  //    pow(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3),3) + 
  //   2*pow(L,3)*x*pow(y,2)*
  //    ((y1 - y2)*(pow(x3,2) + 2*pow(y1 - y3,2)) + 
  //      x2*x3*(y1 - y3) + pow(x1,2)*(2*y1 - y2 - y3) + 
  //      x1*(x2*(-y1 + y3) + x3*(-3*y1 + 2*y2 + y3))))/
  // (L*y*pow(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3),3));
		//	newgradient[ind_tex[0][i][1]*2+1] +=(2*(pow(L,2)*(x1 - x2)*pow(y,2)*(pow(x,2) + pow(y,2))*
  //      (pow(x1,2) + x2*x3 - x1*(x2 + x3) + (y1 - y2)*(y1 - y3)) + 
  //     pow(L,3)*x*pow(y,2)*
  //      (-2*pow(x1,3) + 2*pow(x1,2)*(x2 + 2*x3) + 
  //        x2*(2*pow(x3,2) + pow(y1 - y3,2)) - 
  //        x1*(2*x3*(2*x2 + x3) + (y1 - y3)*(2*y1 - y2 - y3)) + 
  //        x3*(y1 - y2)*(y1 - y3)) + 
  //     pow(L,4)*(x1 - x3)*pow(y,2)*
  //      (pow(x1 - x3,2) + pow(y1 - y3,2)) + 
  //     (pow(x,2) + pow(y,2))*(y1 - y2)*
  //      pow(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3),3) + 
  //     L*x*(y1 - y3)*pow(x3*(-y1 + y2) + x2*(y1 - y3) + 
  //        x1*(-y2 + y3),3)))/
  // (L*y*pow(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3),3));
		//	newgradient[ind_tex[0][i][2]*2] +=(2*(-y1 + y2)*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*
  //    pow(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3,2) - 
  //   2*(-y1 + y2)*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)) - 
  //   2*L*(x*(-x1 + x2) + L*(x1 - x3))*
  //    (x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)))/
  // (L*y*pow(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3),3));
		//	newgradient[ind_tex[0][i][2]*2+1] +=(2*(x1 - x2)*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*
  //    pow(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3,2) - 
  //   2*(x1 - x2)*((pow(x,2) + pow(y,2))*
  //       (pow(x1 - x2,2) + pow(y1 - y2,2)) - 
  //      2*L*x*(pow(x1,2) + x2*x3 - x1*(x2 + x3) + 
  //         (y1 - y2)*(y1 - y3)) + 
  //      pow(L,2)*(pow(x1 - x3,2) + pow(y1 - y3,2)))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)) - 
  //   2*L*(x*(-y1 + y2) + L*(y1 - y3))*
  //    (x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))*
  //    (pow(L,2)*pow(y,2) + 
  //      pow(-(x3*y1) - x1*y2 + x3*y2 + x2*(y1 - y3) + x1*y3,2)))/
  // (L*y*pow(x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3),3));


		//double area = u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2);
		//double area3 = pow ( area, 3 );
		//newgradient[ind_tex[0][i][0]*2] += (2*(area3*(L - xx)*(L*(u0 - u2) + (-u0 + u1)*xx) + (area3*(u0 - u1) + 
  //        pow(L,2)*(L*(v0 - v2) + (-v0 + v1)*xx)*(L*((u0 - u2)*(u1 - u2) + (v0 - v2)*(v1 - v2)) + ((-u0 + u1)*(u1 - u2) - (v0 - v1)*(v1 - v2))*xx))*pow(y,2) + 
  //     pow(L,2)*(v0 - v1)*((u0 - u1)*(u1 - u2) + (v0 - v1)*(v1 - v2))*pow(y,4)))/(area3*L*y);
		//newgradient[ind_tex[0][i][0]*2+1] += (-2*pow(L,2)*(u1 - u2)*pow(y,2)*(pow(L,2)*(pow(u0 - u2,2) + pow(v0 - v2,2)) - 2*L*((u0 - u1)*(u0 - u2) + (v0 - v1)*(v0 - v2))*xx + 
  //      (pow(u0 - u1,2) + pow(v0 - v1,2))*(pow(xx,2) + pow(y,2))) + 2*area3*(pow(L,2)*(v0 - v2) + L*(-2*v0 + v1 + v2)*xx + (v0 - v1)*(pow(xx,2) + pow(y,2))) + 
  //   2*area*pow(L,2)*pow(y,2)*(pow(L,2)*(v0 - v2) + L*(-2*v0 + v1 + v2)*xx + (v0 - v1)*(pow(xx,2) + pow(y,2))))/(area3*L*y);
		//
		//newgradient[ind_tex[0][i][1]*2] += (2*area3*xx*(L*(u0 - u2) + (-u0 + u1)*xx) - 2*(area3*(u0 - u1) + 
  //      pow(L,2)*(L*(v0 - v2) + (-v0 + v1)*xx)*(L*(pow(u0 - u2,2) + pow(v0 - v2,2)) - ((u0 - u1)*(u0 - u2) + (v0 - v1)*(v0 - v2))*xx))*pow(y,2) - 
  //   2*pow(L,2)*(v0 - v1)*(pow(u0,2) + u1*u2 - u0*(u1 + u2) + (v0 - v1)*(v0 - v2))*pow(y,4))/(area3*L*y);
		//newgradient[ind_tex[0][i][1]*2+1] += (2*(area3*xx*(L*(v0 - v2) + (-v0 + v1)*xx) + (area3*(-v0 + v1) + 
  //        pow(L,2)*(L*(u0 - u2) + (-u0 + u1)*xx)*(L*(pow(u0 - u2,2) + pow(v0 - v2,2)) - (pow(u0,2) + u1*u2 - u0*(u1 + u2) + (v0 - v1)*(v0 - v2))*xx))*pow(y,2) + 
  //     pow(L,2)*(u0 - u1)*(pow(u0,2) + u1*u2 - u0*(u1 + u2) + (v0 - v1)*(v0 - v2))*pow(y,4)))/(area3*L*y);
	
		//newgradient[ind_tex[0][i][2]*2] += (2*(area3*(L*(-u0 + u2) + (u0 - u1)*xx) + area*pow(L,2)*(L*(-u0 + u2) + (u0 - u1)*xx)*pow(y,2) + 
  //     L*(v0 - v1)*pow(y,2)*(pow(L,2)*(pow(u0 - u2,2) + pow(v0 - v2,2)) - 2*L*((u0 - u1)*(u0 - u2) + (v0 - v1)*(v0 - v2))*xx + 
  //        (pow(u0 - u1,2) + pow(v0 - v1,2))*(pow(xx,2) + pow(y,2)))))/(area3*y);
		//newgradient[ind_tex[0][i][2]*2+1] += (2*(area3*(L*(-v0 + v2) + (v0 - v1)*xx) - L*(L*(u0 - u2) + (-u0 + u1)*xx)*(L*(pow(u0,2) + u1*u2 - u0*(u1 + u2) + (v0 - v1)*(v0 - v2)) - (pow(u0 - u1,2) + pow(v0 - v1,2))*xx)*
  //      pow(y,2) - L*(u0 - u1)*(pow(u0 - u1,2) + pow(v0 - v1,2))*pow(y,4)))/(area3*y);
	}

	double sum = 0;// fabs ( newgradient[0]);
	for(int i = 0; i < tex_size*2; i++)
	{
//		sum = max ( sum, fabs(newgradient[i]) );
		sum += fabs(newgradient[i]);
	}
	for(int i = 0; i < tex_size*2; i++)
	{
		newgradient[i] /= sum;
	}
#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_grad_int += endT - startT;
#endif
}
#endif

void JasonFull::doQuadrantAddresses(int i, int j)
{
	if(i == j) //is on the diag
	{
		auto item = diagIndexes[i];
		auto p = item.first;
		auto startplus = kvf.m_outerIndex[j + 1];
		auto pplus = startplus + item.second;

		addresses[addressCounter++] = &kvf.m_data.value(p);
		addresses[addressCounter++] = &kvf.m_data.value(pplus);
		addresses[addressCounter++] = &kvf.m_data.value(p+1);
		addresses[addressCounter++] = &kvf.m_data.value(pplus + 1);
	}
	else
	{
		auto start = kvf.m_outerIndex[j];
		auto end = kvf.m_innerNonZeros ? kvf.m_outerIndex[j] + kvf.m_innerNonZeros[j] : kvf.m_outerIndex[j + 1];
		auto p = kvf.m_data.searchLowerIndex(start,end-1, i);

		auto diff = p - start;
		auto startplus = kvf.m_outerIndex[j+1];
		auto pplus = startplus + diff;
		
		addresses[addressCounter++] = &kvf.m_data.value(p);
		addresses[addressCounter++] = &kvf.m_data.value(pplus);
		addresses[addressCounter++] = &kvf.m_data.value(p+1);
		addresses[addressCounter++] = &kvf.m_data.value(pplus + 1);
	}
}

void JasonFull::storeDiagIndexes()
{
	diagIndexes.reserve(tex_size);

	for (int i = 0; i < 2 * tex_size; i++)
	{
		auto start = kvf.m_outerIndex[i];
		auto end = kvf.m_innerNonZeros ? kvf.m_outerIndex[i] + kvf.m_innerNonZeros[i] : kvf.m_outerIndex[i + 1];
		auto p = kvf.m_data.searchLowerIndex(start,end-1, i);

		diagIndexes.push_back(make_pair(p, p-start));
	}
}

void JasonFull::doStoreAddresses (int i1, int i2, int i3 )
{
	doQuadrantAddresses(i1 * 2, i1 * 2);
	//addresses[addressCounter++] = &(kvf.coeffRef(i1 * 2,i1 * 2));
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i1 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i1 * 2 + 1);

	doQuadrantAddresses(i1 * 2, i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i2 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i2 * 2  + 1);

	doQuadrantAddresses(i1 * 2, i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i3 * 2 + 1);

	doQuadrantAddresses(i2 * 2, i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i1 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i1 * 2  + 1);
	//
	doQuadrantAddresses(i2 * 2, i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i2 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i2 * 2 + 1);

	doQuadrantAddresses(i2 * 2, i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i3 * 2 + 1);

	doQuadrantAddresses(i3 * 2, i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i1 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i1 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i1 * 2 + 1);

	doQuadrantAddresses(i3 * 2, i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i2 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i2 * 2 + 1);

	doQuadrantAddresses(i3 * 2, i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i3 * 2 + 1);

	//addresses[addressCounter++] = &(kvf.coeffRef(i1 * 2,i1 * 2));

	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i1 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i1 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i1 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i2 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i1 * 2);
	//
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i1 * 2);
	//
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i1 * 2);
	//
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i1 * 2 + 1);
	//
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i2 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i1 * 2 + 1);
	//
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i2 * 2  + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i1 * 2  + 1);
	//
	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i1 * 2 + 1);

	//addresses[addressCounter++] = &kvf.coeffRef(i1 * 2 + 1,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i1 * 2 + 1);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i2 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i2 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i2 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i2 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i2 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i2 * 2 + 1);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i3 * 2);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i2 * 2 + 1);

	//addresses[addressCounter++] = &kvf.coeffRef(i2 * 2 + 1,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i2 * 2 + 1);

	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i3 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2,i3 * 2 + 1);
	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i3 * 2);

	//addresses[addressCounter++] = &kvf.coeffRef(i3 * 2 + 1,i3 * 2 + 1);
}

void JasonFull::storeAddresses()
{
	double start = get_time();

	addresses = new Eigen::SparseMatrix<double, 0, int>::Scalar* [ 36 * ind_tex[0].size() + 2 * tex_size ];
	//addresses.reserve ( 36 * ind_tex[0].size() + 2 * tex_size );

	storeDiagIndexes();

	for (int i = 0; i < ind_tex[0].size(); i++)
	{
		int i0, i1, i2;

		i0 = ind_tex[0][i][0];
		i1 = ind_tex[0][i][1];
		i2 = ind_tex[0][i][2];

		doStoreAddresses(i0, i1, i2);
	}

	for (int i = 0; i < 2 * tex_size; i++)
	{
		addresses[addressCounter++] = &kvf.m_data.value(diagIndexes[i].first);
		//addresses[addressCounter++] = &(kvf.coeffRef ( i, i ));
	}

	double end = get_time();
	printf("Took %f to store all kvf addresses\n", end-start);
}

void JasonFull::localKVF ( vector<Eigen::Triplet<double> > &trips, double *x, int i1, int i2, int i3 )
{
	double x12 = x [ i1 * 2 ] - x [ i2 * 2 ];
	double y12 = x [ i1 * 2 + 1 ] - x [ i2 * 2 + 1 ];
	double x23 = x [ i2 * 2 ] - x [ i3 * 2 ];
	double y23 = x [ i2 * 2 + 1 ] - x [ i3 * 2 + 1 ];
	double x31 = x [ i3 * 2 ] - x [ i1 * 2 ];
	double y31 = x [ i3 * 2 + 1 ] - x [ i1 * 2 + 1 ];
	double area = sqrt ( 0.5 * fabs ( y12 * x31 -x12 * y31 ) ); // really the sqrt of area
	x12 /= area;
	y12 /= area;
	x23 /= area;
	y23 /= area;
	x31 /= area;
	y31 /= area;

	trips.push_back(Eigen::Triplet<double>(i1 * 2, i1 * 2, x23 * x23 + 2 * y23 * y23));

	trips.push_back(Eigen::Triplet<double>(i1 * 2, i1 * 2 + 1, -x23 * y23));
	trips.push_back(Eigen::Triplet<double>(i1 * 2 + 1, i1 * 2, -x23 * y23));

	trips.push_back(Eigen::Triplet<double>(i1 * 2, i2 * 2, x31*x23+2*y31*y23));
	trips.push_back(Eigen::Triplet<double>(i2 * 2, i1 * 2, x31*x23+2*y31*y23));

	trips.push_back(Eigen::Triplet<double>(i1 * 2, i2 * 2 + 1, -x23*y31));
	trips.push_back(Eigen::Triplet<double>(i2 * 2 + 1, i1 * 2, -x23*y31));
	
	trips.push_back(Eigen::Triplet<double>(i1 * 2, i3 * 2, x12*x23+2*y12*y23));
	trips.push_back(Eigen::Triplet<double>(i3 * 2, i1 * 2, x12*x23+2*y12*y23));
	
	trips.push_back(Eigen::Triplet<double>(i1 * 2, i3 * 2 + 1, -x23 * y12));
	trips.push_back(Eigen::Triplet<double>(i3 * 2 + 1, i1 * 2, -x23 * y12));
	
	trips.push_back(Eigen::Triplet<double>(i1 * 2 + 1, i1 * 2 + 1, 2 * x23 * x23 + y23 * y23));
	
	trips.push_back(Eigen::Triplet<double>(i1 * 2 + 1, i2 * 2, -x31 * y23));
	trips.push_back(Eigen::Triplet<double>(i2 * 2, i1 * 2 + 1, -x31 * y23));
	
	trips.push_back(Eigen::Triplet<double>(i1 * 2 + 1, i2 * 2  + 1, 2 * x23 * x31 + y23 * y31));
	trips.push_back(Eigen::Triplet<double>(i2 * 2 + 1, i1 * 2  + 1, 2 * x23 * x31 + y23 * y31));
	
	trips.push_back(Eigen::Triplet<double>(i1 * 2 + 1, i3 * 2, -x12 * y23));
	trips.push_back(Eigen::Triplet<double>(i3 * 2, i1 * 2 + 1, -x12 * y23));

	trips.push_back(Eigen::Triplet<double>(i1 * 2 + 1, i3 * 2 + 1, 2 * x12 * x23 + y12 * y23));
	trips.push_back(Eigen::Triplet<double>(i3 * 2 + 1, i1 * 2 + 1, 2 * x12 * x23 + y12 * y23));

	trips.push_back(Eigen::Triplet<double>(i2 * 2, i2 * 2, x31 * x31 + 2 * y31 * y31));

	trips.push_back(Eigen::Triplet<double>(i2 * 2, i2 * 2 + 1, -x31 * y31));
	trips.push_back(Eigen::Triplet<double>(i2 * 2 + 1, i2 * 2, -x31 * y31));

	trips.push_back(Eigen::Triplet<double>(i2 * 2, i3 * 2, x12 * x31 + 2 * y12 * y31));
	trips.push_back(Eigen::Triplet<double>(i3 * 2, i2 * 2, x12 * x31 + 2 * y12 * y31));

	trips.push_back(Eigen::Triplet<double>(i2 * 2, i3 * 2 + 1, -x31 * y12));
	trips.push_back(Eigen::Triplet<double>(i3 * 2 + 1, i2 * 2, -x31 * y12));

	trips.push_back(Eigen::Triplet<double>(i2 * 2 + 1, i2 * 2 + 1, 2 * x31 * x31 + y31 * y31));

	trips.push_back(Eigen::Triplet<double>(i2 * 2 + 1, i3 * 2, -x12 * y31));
	trips.push_back(Eigen::Triplet<double>(i3 * 2, i2 * 2 + 1, -x12 * y31));

	trips.push_back(Eigen::Triplet<double>(i2 * 2 + 1, i3 * 2 + 1, 2 * x12 * x31 + y12 * y31));
	trips.push_back(Eigen::Triplet<double>(i3 * 2 + 1, i2 * 2 + 1, 2 * x12 * x31 + y12 * y31));

	trips.push_back(Eigen::Triplet<double>(i3 * 2, i3 * 2, x12 * x12 + 2 * y12 * y12));

	trips.push_back(Eigen::Triplet<double>(i3 * 2, i3 * 2 + 1, -x12 * y12));
	trips.push_back(Eigen::Triplet<double>(i3 * 2 + 1, i3 * 2, -x12 * y12));

	trips.push_back(Eigen::Triplet<double>(i3 * 2 + 1, i3 * 2 + 1, 2 * x12 * x12 + y12 * y12));
}

//void JasonFull::localKVFAdd ( double *x, int i1, int i2, int i3 )
//{
//	double x12 = x [ i1 * 2 ] - x [ i2 * 2 ];
//	double y12 = x [ i1 * 2 + 1 ] - x [ i2 * 2 + 1 ];
//	double x23 = x [ i2 * 2 ] - x [ i3 * 2 ];
//	double y23 = x [ i2 * 2 + 1 ] - x [ i3 * 2 + 1 ];
//	double x31 = x [ i3 * 2 ] - x [ i1 * 2 ];
//	double y31 = x [ i3 * 2 + 1 ] - x [ i1 * 2 + 1 ];
//	double area = sqrt ( 0.5 * fabs ( y12 * x31 -x12 * y31 ) ); // really the sqrt of area
//	x12 /= area;
//	y12 /= area;
//	x23 /= area;
//	y23 /= area;
//	x31 /= area;
//	y31 /= area;
//
//	kvf.coeffRef(i1 * 2,i1 * 2) +=  (x23 * x23 + 2 * y23 * y23 );
//
//	kvf.coeffRef(i1 * 2,i1 * 2 + 1) += -x23 * y23;
//	kvf.coeffRef(i1 * 2 + 1,i1 * 2) +=  -x23 * y23;
//
//	kvf.coeffRef(i1 * 2,i2 * 2) += x31*x23+2*y31*y23;
//	kvf.coeffRef(i2 * 2,i1 * 2) += x31*x23+2*y31*y23;
//
//	kvf.coeffRef(i1 * 2,i2 * 2 + 1) += -x23*y31;
//	kvf.coeffRef(i2 * 2 + 1,i1 * 2) += -x23*y31;
//	
//	kvf.coeffRef(i1 * 2,i3 * 2) += x12*x23+2*y12*y23;
//	kvf.coeffRef(i3 * 2,i1 * 2) += x12*x23+2*y12*y23;
//	
//	kvf.coeffRef(i1 * 2,i3 * 2 + 1) += -x23 * y12;
//	kvf.coeffRef(i3 * 2 + 1,i1 * 2) += -x23 * y12;
//	
//	kvf.coeffRef(i1 * 2 + 1,i1 * 2 + 1) +=  (2 * x23 * x23 + y23 * y23);
//	
//	kvf.coeffRef(i1 * 2 + 1,i2 * 2) += -x31 * y23;
//	kvf.coeffRef(i2 * 2,i1 * 2 + 1) += -x31 * y23;
//	
//	kvf.coeffRef(i1 * 2 + 1,i2 * 2  + 1) += 2 * x23 * x31 + y23 * y31;
//	kvf.coeffRef(i2 * 2 + 1,i1 * 2  + 1) += 2 * x23 * x31 + y23 * y31;
//	
//	kvf.coeffRef(i1 * 2 + 1,i3 * 2) += -x12 * y23;
//	kvf.coeffRef(i3 * 2,i1 * 2 + 1) += -x12 * y23;
//
//	kvf.coeffRef(i1 * 2 + 1,i3 * 2 + 1) += 2 * x12 * x23 + y12 * y23;
//	kvf.coeffRef(i3 * 2 + 1,i1 * 2 + 1) += 2 * x12 * x23 + y12 * y23;
//
//	kvf.coeffRef(i2 * 2,i2 * 2) +=  (x31 * x31 + 2 * y31 * y31);
//
//	kvf.coeffRef(i2 * 2,i2 * 2 + 1) += -x31 * y31;
//	kvf.coeffRef(i2 * 2 + 1,i2 * 2) += -x31 * y31;
//
//	kvf.coeffRef(i2 * 2,i3 * 2) += x12 * x31 + 2 * y12 * y31;
//	kvf.coeffRef(i3 * 2,i2 * 2) += x12 * x31 + 2 * y12 * y31;
//
//	kvf.coeffRef(i2 * 2,i3 * 2 + 1) += -x31 * y12;
//	kvf.coeffRef(i3 * 2 + 1,i2 * 2) += -x31 * y12;
//
//	kvf.coeffRef(i2 * 2 + 1,i2 * 2 + 1) +=  (2 * x31 * x31 + y31 * y31);
//
//	kvf.coeffRef(i2 * 2 + 1,i3 * 2) += -x12 * y31;
//	kvf.coeffRef(i3 * 2,i2 * 2 + 1) += -x12 * y31;
//
//	kvf.coeffRef(i2 * 2 + 1,i3 * 2 + 1) += 2 * x12 * x31 + y12 * y31;
//	kvf.coeffRef(i3 * 2 + 1,i2 * 2 + 1) += 2 * x12 * x31 + y12 * y31;
//
//	kvf.coeffRef(i3 * 2,i3 * 2) +=  (x12 * x12 + 2 * y12 * y12);
//
//	kvf.coeffRef(i3 * 2,i3 * 2 + 1) += -x12 * y12;
//	kvf.coeffRef(i3 * 2 + 1,i3 * 2) += -x12 * y12;
//
//	kvf.coeffRef(i3 * 2 + 1,i3 * 2 + 1) +=  (2 * x12 * x12 + y12 * y12);
//}

void JasonFull::localKVFAdd ( double *x, int i1, int i2, int i3 )
{
	double x12 = x [ i1 * 2 ] - x [ i2 * 2 ];
	double y12 = x [ i1 * 2 + 1 ] - x [ i2 * 2 + 1 ];
	double x23 = x [ i2 * 2 ] - x [ i3 * 2 ];
	double y23 = x [ i2 * 2 + 1 ] - x [ i3 * 2 + 1 ];
	double x31 = x [ i3 * 2 ] - x [ i1 * 2 ];
	double y31 = x [ i3 * 2 + 1 ] - x [ i1 * 2 + 1 ];
	double area = sqrt ( 0.5 *  ( y12 * x31 -x12 * y31 ) ); // really the sqrt of area
	x12 /= area;
	y12 /= area;
	x23 /= area;
	y23 /= area;
	x31 /= area;
	y31 /= area;

	*(addresses[addressCounter++]) +=  (x23 * x23 + 2 * y23 * y23 ); //i1, i1
	*(addresses[addressCounter++]) += -x23 * y23; //i1, i1 + 1
	*(addresses[addressCounter++]) +=  -x23 * y23; //i1 + 1, i1
	*(addresses[addressCounter++]) +=  (2 * x23 * x23 + y23 * y23); //i1 + 1, i1 + 1

	*(addresses[addressCounter++]) += x31*x23+2*y31*y23; //i1, i2
	*(addresses[addressCounter++]) += -x23*y31; //i1, i2 + 1
	*(addresses[addressCounter++]) += -x31 * y23; //i1 + 1, i2 
	*(addresses[addressCounter++]) += 2 * x23 * x31 + y23 * y31; //i1 + 1, i2 + 1

	*(addresses[addressCounter++]) += x12*x23+2*y12*y23; //i1, i3
	*(addresses[addressCounter++]) += -x23 * y12; //i1, i3 + 1
	*(addresses[addressCounter++]) += -x12 * y23; //i1 + 1, i3
	*(addresses[addressCounter++]) += 2 * x12 * x23 + y12 * y23; //i1 + 1, i3 + 1

	*(addresses[addressCounter++]) += x31*x23+2*y31*y23; //i2, i1
	*(addresses[addressCounter++]) += -x31 * y23; //i2, i1 + 1
	*(addresses[addressCounter++]) += -x23*y31; //i2 + 1, i1
	*(addresses[addressCounter++]) += 2 * x23 * x31 + y23 * y31; //i2 + 1, i1 + 1

	*(addresses[addressCounter++]) +=  (x31 * x31 + 2 * y31 * y31); //i2, i2
	*(addresses[addressCounter++]) += -x31 * y31; //i2, i2 + 1
	*(addresses[addressCounter++]) += -x31 * y31; //i2 + 1, i2
	*(addresses[addressCounter++]) +=  (2 * x31 * x31 + y31 * y31); //i2 + 1, i2 + 1

	*(addresses[addressCounter++]) += x12 * x31 + 2 * y12 * y31; //i2, i3
	*(addresses[addressCounter++]) += -x31 * y12; //i2, i3 + 1
	*(addresses[addressCounter++]) += -x12 * y31; //i2 + 1, i3
	*(addresses[addressCounter++]) += 2 * x12 * x31 + y12 * y31; //i2 + 1, i3 + 1

	*(addresses[addressCounter++]) += x12*x23+2*y12*y23; //i3, i1
	*(addresses[addressCounter++]) += -x12 * y23; //i3, i1 + 1
	*(addresses[addressCounter++]) += -x23 * y12; //i3 + 1, i1
	*(addresses[addressCounter++]) += 2 * x12 * x23 + y12 * y23; //i3 + 1, i1 + 1

	*(addresses[addressCounter++]) += x12 * x31 + 2 * y12 * y31; //i3, i2
	*(addresses[addressCounter++]) += -x12 * y31; //i3, i2 + 1
	*(addresses[addressCounter++]) += -x31 * y12; //i3 + 1, i2
	*(addresses[addressCounter++]) += 2 * x12 * x31 + y12 * y31; //i3 + 1, i2 + 1

	*(addresses[addressCounter++]) +=  (x12 * x12 + 2 * y12 * y12); //i3, i3
	*(addresses[addressCounter++]) += -x12 * y12; //i3, i3 + 1
	*(addresses[addressCounter++]) += -x12 * y12; //i3 + 1, i3
	*(addresses[addressCounter++]) +=  (2 * x12 * x12 + y12 * y12); //i3 + 1, i3 + 1

	/**(addresses[addressCounter++]) +=  (x23 * x23 + 2 * y23 * y23 );

	*(addresses[addressCounter++]) += -x23 * y23;
	*(addresses[addressCounter++]) +=  -x23 * y23;

	*(addresses[addressCounter++]) += x31*x23+2*y31*y23;
	*(addresses[addressCounter++]) += x31*x23+2*y31*y23;

	*(addresses[addressCounter++]) += -x23*y31;
	*(addresses[addressCounter++]) += -x23*y31;
	
	*(addresses[addressCounter++]) += x12*x23+2*y12*y23;
	*(addresses[addressCounter++]) += x12*x23+2*y12*y23;
	
	*(addresses[addressCounter++]) += -x23 * y12;
	*(addresses[addressCounter++]) += -x23 * y12;
	
	*(addresses[addressCounter++]) +=  (2 * x23 * x23 + y23 * y23);
	
	*(addresses[addressCounter++]) += -x31 * y23;
	*(addresses[addressCounter++]) += -x31 * y23;
	
	*(addresses[addressCounter++]) += 2 * x23 * x31 + y23 * y31;
	*(addresses[addressCounter++]) += 2 * x23 * x31 + y23 * y31;
	
	*(addresses[addressCounter++]) += -x12 * y23;
	*(addresses[addressCounter++]) += -x12 * y23;

	*(addresses[addressCounter++]) += 2 * x12 * x23 + y12 * y23;
	*(addresses[addressCounter++]) += 2 * x12 * x23 + y12 * y23;

	*(addresses[addressCounter++]) +=  (x31 * x31 + 2 * y31 * y31);

	*(addresses[addressCounter++]) += -x31 * y31;
	*(addresses[addressCounter++]) += -x31 * y31;

	*(addresses[addressCounter++]) += x12 * x31 + 2 * y12 * y31;
	*(addresses[addressCounter++]) += x12 * x31 + 2 * y12 * y31;

	*(addresses[addressCounter++]) += -x31 * y12;
	*(addresses[addressCounter++]) += -x31 * y12;

	*(addresses[addressCounter++]) +=  (2 * x31 * x31 + y31 * y31);

	*(addresses[addressCounter++]) += -x12 * y31;
	*(addresses[addressCounter++]) += -x12 * y31;

	*(addresses[addressCounter++]) += 2 * x12 * x31 + y12 * y31;
	*(addresses[addressCounter++]) += 2 * x12 * x31 + y12 * y31;

	*(addresses[addressCounter++]) +=  (x12 * x12 + 2 * y12 * y12);

	*(addresses[addressCounter++]) += -x12 * y12;
	*(addresses[addressCounter++]) += -x12 * y12;

	*(addresses[addressCounter++]) +=  (2 * x12 * x12 + y12 * y12);*/
}

void JasonFull::computeKVF_First_Time(double *x)
{
	double start = get_time();

	vector<Eigen::Triplet<double> > trips;
	trips.reserve(ind_tex[0].size() * 36+2*tex_size);

	//FILE *fptr = fopen ( "mesh.txt", "wt" );
	//fprintf(fptr, "{{");
	//for ( int i = 0; i < ind_tex[0].size() - 1; i++ )
	//{
	//	fprintf ( fptr, "{%d,%d,%d},", ind_tex[0][i][0]+1, ind_tex[0][i][1]+1, ind_tex[0][i][2]+1 );
	//}
	//fprintf ( fptr, "{%d,%d,%d}},", ind_tex[0][ind_tex[0].size() - 1][0]+1, ind_tex[0][ind_tex[0].size() - 1][1]+1, ind_tex[0][ind_tex[0].size() - 1][2]+1 );
	//fprintf(fptr,"{");
	//for ( int i = 0; i < 2 * tex_size - 2; i+= 2 )
	//{
	//	fprintf(fptr,"{%f,%f},", x[i], x[i+1]);
	//}
	//fprintf(fptr,"{%f,%f}}}\n", x[2*tex_size-2], x[2*tex_size-1]);
	//fclose ( fptr);

	for (int i = 0; i < ind_tex[0].size(); i++)
	{
		int i0, i1, i2;

		i0 = ind_tex[0][i][0];
		i1 = ind_tex[0][i][1];
		i2 = ind_tex[0][i][2];

		if (i0 == -1 && i1 == -1 && i2 == -1)
			continue;

		localKVF( trips, x, i0, i1, i2 );
	}

	for (int i = 0; i < 2 * tex_size; i++)
	{
		trips.push_back(Eigen::Triplet<double>(i, i, 0.0001));
	}

	kvf.resize(2 * tex_size, 2 * tex_size);
	kvf.setFromTriplets(trips.begin(), trips.end());

	//fptr = fopen("kvf.txt", "wt");
	//fprintf(fptr, "{");
	//for (int k = 0; k < kvf.outerSize(); ++k)
	//{
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(kvf, k); it; ++it)
	//	{
	//		fprintf(fptr,"{%d,%d,%f},", it.row(), k, it.value());
	//	}
	//}
	//fprintf(fptr, "}");
	//fclose(fptr);

	double startT = get_time();

	kvfSolver.analyzePattern(kvf);
	double analyzeT = get_time();
	kvfSolver.factorize(kvf);
	double end = get_time();
	printf("Took %f to compute kvf, %f to analyze, and %f to factorization with %f total\n", startT - start, analyzeT - startT, end - analyzeT, end-start);
}

void JasonFull::innerFlapOptimization(double *x)
{
	for(unsigned int i = 0; i < g.flapLoopup.size(); i++)
	{
		auto alpha = g.flapLoopup[i].alpha;
		auto L = g.flapLoopup[i].L;

		auto verts = g.flapLoopup[i].verts;

		auto p1x = x[verts[1] * 2];
		auto p1y = x[verts[1] * 2 + 1];

		auto p2x = x[verts[2] * 2];
		auto p2y = x[verts[2] * 2 + 1];

		auto length = sqrt((p2x - p1x) * (p2x - p1x) + (p2y - p1y) * (p2y - p1y) );

		auto p2p1perpx = -(p2y - p1y);
		auto p2p1perpy = p2x - p1x;

		//auto p2p1perp = p2.uvpos - p1.uvpos;
		//auto temp = p2p1perp[0];
		//p2p1perp[0] = -p2p1perp[1];
		//p2p1perp[1] = temp;
		
		//p3.uvpos = (1 - alpha)*p1.uvpos + alpha * p2.uvpos + (p2p1perp * L / (p2.uvpos - p1.uvpos).length());

		x[verts[0] * 2] = (1 - alpha) * p1x + alpha * p2x + (p2p1perpx * L / length);
		x[verts[0] * 2 + 1] = (1 - alpha) * p1y + alpha * p2y + (p2p1perpy * L / length);
	}
}


void JasonFull::computeKVF(double *x)
{
	double start = get_time();

	// does this work?
	//kvf *= 0; // ?
//	kvf.setZero();
	memset ( &(kvf.data().value(0)), 0, sizeof ( double ) * kvf.data().size ( ) );
//	memset ( kvf.valuePtr ( ), 0, kvf.data ( ).value.size ( ) * sizeof ( double ) );
//	kvf.setZero ( );
	addressCounter = 0;

	innerFlapOptimization(x);

	for (int i = 0; i < ind_tex[0].size(); i++)
	{
		int i0, i1, i2;

		i0 = ind_tex[0][i][0];
		i1 = ind_tex[0][i][1];
		i2 = ind_tex[0][i][2];

		localKVFAdd( x, i0, i1, i2 );
	}
//	kvf = kvf + kvf.transpose ( );

	for (int i = 0; i < 2 * tex_size; i++)
	{
		//kvf.coeffRef ( i, i ) += 0.0001;
		*(addresses[addressCounter++]) += 0.0001;
	}

	double startT = get_time();

	kvfSolver.factorize(kvf);
	double end = get_time();
	printf("Took %f to compute kvf, %f to factorization with %f total\n", startT - start, end - startT, end-start);
}

void JasonFull::calc_interior_gradient( double *pos, double *g )
{
	//BarrierMethod::calc_interior_gradient(pos, g);

	//Calc interior gradient here
	barrier_jason_gradient(pos, g);
}
void JasonFull::display_error( double *e, double *pos )
{
	BarrierMethod::display_error(e, pos);

	int size = this->tex_size / 2.0f;

	for(int i = 0; i < size; i++)
	{
		e[i] = 0;
	}

	for(int ii = 0; ii < ind_tex[0].size(); ii++)
	{
		int i0, i1, i2;
		i0 = indchange[0][ind_tex[0][ii][0]];
		i1 = indchange[0][ind_tex[0][ii][1]];
		i2 = indchange[0][ind_tex[0][ii][2]];

		if(i0 == -1 && i1 == -1 && i2 == -1)
			continue;

		double fx =e_tri_func(iso_tris[0][ii][1][0],
							iso_tris[0][ii][2][0],
							iso_tris[0][ii][2][1],
							pos[2*i0],
							pos[2*i0+1],
							pos[2*i1],
							pos[2*i1+1],
							pos[2*i2],
							pos[2*i2+1]);

		e[ii] += fx;
	}
}

/*************************************************/
//Functions
/*************************************************/

void JasonFull::set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it)
{
	BarrierMethod::set_collapse_level(totalTex, curTex, face_ind, ic, bo, it);

	if (firstTime)
	{
		threadGrad = new double*[omp_get_max_threads()];

		firstTime = false;
		for (int i = 0; i < omp_get_max_threads(); i++)
		{
			threadGrad[i] = new double[totalTex * 2 * 2]; // additional factor of 2 for quad and nonlin gradient calculation
		}
	}

}