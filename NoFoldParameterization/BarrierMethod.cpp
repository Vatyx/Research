#include "BarrierMethod.h"
#include <iostream>

#include "vect.h"
#include "global.h"

#include <omp.h>
//Scott method for Max Param
//


void BarrierMethod::init_lbfgs(int n)
{
}
pair<double,double> findBezRoot(double a, double b, double c, double tmin, double tmax)
{
	//printf("\tf %.10g, %.10g, %.10g\n",a,b,c);
	if( (a > 0 && b > 0 && c > 0) || (a < 0 && b < 0 && c < 0) )
	{// No possible root
		return pair<double, double>(-1,-1);
	}

	if( (abs(b - (a+c)/2.0f) < 10e-10 && (a*c <= 0 || a*b >= 0))  || (tmax - tmin < 10e-10) ) 
	{
		if(abs(a-c) < 10e-10)
		{
			return pair<double, double> (tmin, tmax);
		}
		double t = a / (a-c);

		/*if( t < tmin || t > tmax )
		{
			return pair<double, double>(-1,-1);
		}*/

		if( t < 0 || t > 1 )
		{
			return pair<double, double>(-1,-1);
		}

		return pair<double, double>(t*tmax + (1-t)*tmin, -1);
	}

	//printf("\t%.10g, %.10g, %.10g\n",a,b,c);
	//subd
	pair<double, double> p1 = findBezRoot(a, (a+b) / 2.0f, (a + 2*b + c)/ 4.0f, tmin, (tmax+tmin)/2.0f );
	pair<double, double> p2 = findBezRoot((a+2*b+c)/4.0f, (b+c)/2.0f, c, (tmax+tmin)/2.0f, tmax);

	if(p1.first == -1)return p2;
	if(p2.first == -1)return p1;

	return pair<double,double>(p1.first, p2.first);

}

pair<double, double> findIntersect(vect2d e1, vect2d v1, vect2d e2, vect2d v2,
	vect2d p, vect2d vp, double tmax)
{
	double a = -1*e2[1]*p[0] + e1[1]* (-1*e2[0] + p[0]) + e1[0]* (e2[1] - p[1]) + e2[0] *p[1];

	double b = -1*e1[1]* e2[0] + e1[0]* e2[1] + e1[1]* p[0] - e2[1] *p[0] - e1[0] *p[1] + e2[0]* p[1] + 
		.5f* tmax* (-1*e2[0]* v1[1] + p[0]* v1[1] - e1[1] *v2[0] + p[1] *(-1*v1[0] + v2[0]) + e1[0]* v2[1] - p[0] *v2[1] + 
         e2[1] *(v1[0] - vp[0]) + e1[1]* vp[0] - e1[0]* vp[1] +  e2[0]* vp[1]);

	double c = -1*e1[1]* e2[0] + e1[0]* e2[1] + e1[1] *p[0] - e2[1]* p[0] - e1[0]* p[1] + e2[0]* p[1] + 
		tmax* (-1*e2[0]* v1[1] + p[0]* v1[1] - e1[1]* v2[0] + p[1]* (-1*v1[0] + v2[0]) + e1[0] *v2[1] - 
		p[0] *v2[1] +  e2[1] *(v1[0] - vp[0]) + e1[1] *vp[0] - e1[0] *vp[1] + e2[0]* vp[1]) + 
    tmax*tmax* (-1*v2[1]* vp[0] + v1[1] *(-1*v2[0] + vp[0]) + v1[0]* (v2[1] - vp[1]) + v2[0] *vp[1]);

	pair<double, double> roots = findBezRoot(a,b,c, 0,tmax);
	//printf("roots... %.10g, %.10g\n", roots.first, roots.second);

	if(roots.first != -1 && roots.second != -1)
	{//Range of roots
		/*double temp = (p-e1).dot(e2-e1) / (e2-e1).dot(e2-e1);
		if( 0 <= temp && temp <= 1)
		{
			printf("!!!\n");
			return pair<double, double> (0,-1);
		}*/
			pair<double,double> tt1 = findBezRoot((e1-p).dot(e1-p), 
												 (p-e1).dot(p-e1) + tmax*(v1-vp).dot(e1-p),
												 (p-e1).dot(p-e1) + 2*tmax*(v1-vp).dot(e1-p) + tmax*tmax*(v1-vp).dot(v1-vp),
												 0,
												 tmax);

			pair<double,double> tt2 = findBezRoot((e2-p).dot(e2-p), 
												 (p-e2).dot(p-e2) + tmax*(v2-vp).dot(e2-p),
												 (p-e2).dot(p-e2) + 2*tmax*(v2-vp).dot(e2-p) + tmax*tmax*(v2-vp).dot(v2-vp),
												 0,
												 tmax);

			if(tt1.first == -1)
				return tt2;
			if(tt2.first == -1)
				return tt1;

			return pair<double, double>(  min((double) tt1.first, (double)tt2.first ) , -1  );
	}

	int size = 0;
	if(roots.first == -1)
		size = 0;
	else if(roots.second == -1)
		size = 1;
	else
		size = 2;

	//printf("size = %d\n", size);
	for(int i = 0; i < size; i++ )
	{
		if(i == 0)
		{
		vect2d ee1 = e1 + v1*roots.first;
		vect2d ee2 = e2 + v2*roots.first;
		vect2d pp = p+vp*roots.first;

		double temp = (pp-ee1).dot(ee2-ee1) / (ee2-ee1).dot(ee2-ee1);
		//printf("t1 = %.10g\n", temp);
		if(0 <= temp && temp <= 1)
			return pair<double,double> (roots.first, -1);
		}
		else
		{
			vect2d ee1 = e1 + v1*roots.second;
		vect2d ee2 = e2 + v2*roots.second;
		vect2d pp = p+vp*roots.second;

		double temp = (pp-ee1).dot(ee2-ee1) / (ee2-ee1).dot(ee2-ee1);
		//printf("t2 = %.10g\n", temp);
		if(0 <= temp && temp <= 1)
			return pair<double,double> (roots.second, -1);
		}
	}

	//printf("fuck\n");
	return pair<double, double>(-1,-1);
}
//////////////////////////////////////////////////////////////////////
//Jason method for Max Param

double BarrierMethod::max_parameter(double *f,
    double *x,
	double *gg,
    double *q,
	int n)
{
	return 0;
}
void BarrierMethod::is_movement(double dis, double grad)
{
}

void BarrierMethod::change_sA(double &sA, double *x, double *q, int i, int ip1, int j)
{
	if(i != j && ip1 != j)
	{
		double u0,u1,u2,v0,v1,v2,gx0,gx1,gx2,gy0,gy1,gy2;

		int i0, i1, i2;
		i0 = indchange[0][i];
		i1 = indchange[0][ip1];
		i2 = indchange[0][j];
		
		//if(i0 == -1)
		/*{
			u0 = full_mesh[2*i];
			v0 = full_mesh[2*i+1];
			gx0 = 0;
			gy0 = 0;	
		}
		else*/
		{
			u0 = x[2*i0];
			v0 = x[2*i0+1];
			gx0 = q[2*i0];
			gy0 = q[2*i0 + 1];
		}

		//if(i1 == -1)
		/*{
			u1 = full_mesh[2*ip1];
			v1 = full_mesh[2*ip1+1];
			gx1 = 0;
			gy1 = 0;
		}
		else*/
		{
			u1 = x[2*i1];
			v1 = x[2*i1+1];
			gx1 = q[2*i1];
			gy1 = q[2*i1 + 1];
		}

		//if(i2 == -1)
		/*{
			u2 = full_mesh[2*j];
			v2 = full_mesh[2*j+1];
			gx2 = 0;
			gy2 = 0;
		}
		else*/
		{
			u2 = x[2*i2];
			v2 = x[2*i2+1];
			gx2 = q[2*i2];
			gy2 = q[2*i2 + 1];
		}

			
		/*u0 = x[2*i];
		v0 = x[2*i+1];
		gx0 = q[2*i];
		gy0 = q[2*i + 1];

		u1 = x[2*ip1];
		v1 = x[2*ip1+1];
		gx1 = q[2*ip1];
		gy1 = q[2*ip1 + 1];

		u2 = x[2*j];
		v2 = x[2*j+1];
		gx2 = q[2*j];
		gy2 = q[2*j + 1];*/

		//SWITCH GRADIENT WITH Q
		//Q is flipped
		
		
		

		

		double c = -(1e-20) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1*v0 + v2));
		double b = gy2 *(u0 - u1) + gy0 *(u1 - u2) + gy1 *(u2 -u0) + 
					gx2 *(v1 -v0) + gx1 *(v0 - v2) + gx0 *(v2 -v1);
		double a = gx2* (-gy0 + gy1) + gx0* (-gy1 + gy2) + gx1* (gy0 - gy2);

		/*if ( g.printInfo && i == 762 && ip1 == 767 && j == 2916)
		{
			printf ( "a = {%f, %f}; b = {%f, %f}; p = {%f, %f};\n",u0,v0,u1,v1,u2,v2);
			printf("ga = {%f, %f}; gb = {%f, %f}; gp = {%f, %f};\n", gx0, gy0, gx1, gy1, gx2,gy2);
		}*/

		//Check disc
		double disc = b*b - 4*a*c;

		if( disc < 0 )
		{
			//Do nothing
			//printf("neg disc\n");
		}
		else
		{
			if(a == 0)
			{
				double ta = -c / b;
				if(ta > 0 && ta < sA)
				{
					vect2d a;
					a[0] = u0  + ta*gx0;
					a[1] = v0  + ta*gy0;

					vect2d b;
					b[0] = u1 + ta*gx1;
					b[1] = v1  + ta*gy1;

					vect2d p;
					p[0] = u2  + ta*gx2;
					p[1] = v2  + ta*gy2;

					double t = ((p - a).dot(b-a)) / (a - b).length2();

					if (t >= 0.0 && t <= 1 /*&& t2 >= 0.0 && t2 <= 1*/) 
					{
						sA = ta;
					}
				}
			}
			else
			{
			//Determine which case
			double ta = (-1*b - sqrt( disc )) / (2*a);

			if(ta > 0 && ta < sA)
			{
				/*vect2d a;
				a[0] = x[i*2]  + ta*q[i*2];
				a[1] = x[i*2+1]  + ta*q[i*2+1];

				vect2d b;
				b[0] = x[ip1*2]  + ta*q[ip1*2];
				b[1] = x[ip1*2+1]  + ta*q[ip1*2+1];

				vect2d p;
				p[0] = x[j*2]  + ta*q[j*2];
				p[1] = x[j*2+1]  + ta*q[j*2+1];*/

				vect2d a;
				a[0] = u0  + ta*gx0;
				a[1] = v0  + ta*gy0;

				vect2d b;
				b[0] = u1 + ta*gx1;
				b[1] = v1  + ta*gy1;

				vect2d p;
				p[0] = u2  + ta*gx2;
				p[1] = v2  + ta*gy2;

				double t = ((p - a).dot(b-a)) / (a - b).length2();

				if (t >= 0.0 && t <= 1 /*&& t2 >= 0.0 && t2 <= 1*/) 
				{
					sA = ta;
				}
			}

			ta = (-1*b+ sqrt( disc )) / (2*a);
			if(ta > 0 && ta < sA)
			{
				/*vect2d a;
				a[0] = x[i*2]  + ta*q[i*2];
				a[1] = x[i*2+1]  + ta*q[i*2+1];

				vect2d b;
				b[0] = x[ip1*2]  +ta*q[ip1*2];
				b[1] = x[ip1*2+1]   + ta*q[ip1*2+1];

				vect2d p;
				p[0] = x[j*2]  + ta*q[j*2];
				p[1] = x[j*2+1]  + ta*q[j*2+1];*/

				vect2d a;
				a[0] = u0  + ta*gx0;
				a[1] = v0  + ta*gy0;

				vect2d b;
				b[0] = u1 + ta*gx1;
				b[1] = v1  + ta*gy1;

				vect2d p;
				p[0] = u2  + ta*gx2;
				p[1] = v2  + ta*gy2;

				double t = ((p - a).dot(b-a)) / (a - b).length2();

				if (t >= 0.0 && t <= 1 /*&& t2 >= 0.0 && t2 <= 1*/) 
				{
					sA = ta;
				}

			}
		}}
	}
}


void BarrierMethod::set_collapse_level(int totalTex, int curTex, vector<vect3i> *face_ind, vector<int> *ic, vector<vector<int>> *bo, vector<vector<vect2d>> *it)
{
	ParamMethod::set_collapse_level(totalTex, curTex, face_ind, ic, bo, it);

	total_tex_size = totalTex;
	tex_size = curTex;
	//printf("curTex = %d\n", curTex);


	indchange = ic;
	/*indchange.clear();
	for (int i = 0; i < ic.size(); i++)
	{
		int t = ic[i];
		indchange.push_back(t);
	}*/

	boundaryorder = bo;
	/*boundaryorder.clear();
	for (int i = 0; i < bo.size(); i++)
	{
		vector<int> temp;
		for (int j = 0; j < bo[i].size(); j++)
		{
			int t = bo[i][j];
			temp.push_back(t);
		}
		boundaryorder.push_back(temp);
	}*/

	iso_tris = it;
	for (int i = 0; i < it->size(); i++)
	{
		iso_pre t;
		t = iso_pre( it[0][i][1][0],it[0][i][2][0],it[0][i][2][1] );

		precomp.push_back(t);
	}
	printf("precomputation = %d\n", precomp.size());


	

	/*iso_tris.clear();
	for (int i = 0; i < it.size(); i++)
	{
		vector<vect2d> temp;
		for (int j = 0; j < it[i].size(); j++)
		{
			vect2d t;
			t[0] = it[i][j][0];
			t[1] = it[i][j][1];
			temp.push_back(t);
		}
		iso_tris.push_back(temp);
	}*/
}

BarrierMethod::BarrierMethod(string n, int num_uv, vector<vector<int>> bo, vector<vector<vect2d>> it)
	: ParamMethod(n)
{
	func_exp = STARTING_EXP;
	alpha_val = ALPHA;

//	boundaryorder = bo;
//	iso_tris = it;

	tex_size = num_uv;
	time_func_bound = time_grad_bound = time_funcgrad_intervalsetup = time_sA_intervalsetup = time_sA_sort = time_sA_bound = 0;

}


void BarrierMethod::print_timings()
{
	ParamMethod::print_timings();
	printf("Barrier\n");
	printf("time for max param = %f\n", time_sA_intervalsetup+time_sA_sort+time_sA_bound);
	printf("\t interval %f\n", time_sA_intervalsetup);
	printf("\t sort %f\n", time_sA_sort);
	printf("\t param %f\n", time_sA_bound);
	printf("time for func and grad eval = %f\n",time_func_bound+ time_grad_bound+time_funcgrad_intervalsetup);
	printf("\t grid %f\n", time_funcgrad_intervalsetup);
	printf("\t func %f\n", time_func_bound);
	printf("\t grad %f\n", time_grad_bound);
}

int BarrierMethod::run(int iters, double *pos, double *ans)
{
	ParamMethod::run(iters, pos, ans);
//	cout << "Starting Barrier Method\n";

	mins_interval.clear();
	maxs_interval.clear();
	pmins_interval.clear();
	pmaxs_interval.clear();

	int tot_boundary = 0;
	for (int i = 0; i < boundaryorder[0].size(); i++)
		tot_boundary += boundaryorder[0][i].size();

	mins_interval.resize(tot_boundary);
	maxs_interval.resize(tot_boundary);

	pmins_interval.resize(tot_boundary);
	pmaxs_interval.resize(tot_boundary);

	//work here
	return 0;
}

bool interval_sort_func( pair<int, double> i, pair<int, double> j) { return (i.second < j.second); }


#ifdef USE_INTERVAL_TREE
double BarrierMethod::max_param( double *f, double *x, double *gg,  double *q, int n )
{
	printf ( "%d threads\n", omp_get_num_threads ( ) );

	//Get intervals
#ifdef DO_FULL_TIMINGS
		double startT = get_time();
#endif
		double sA = currentsA;

		int temp;
	/*	if(!_finite(gg[0]))
		{
			printf("barrier gg\n");
			scanf("%d",temp);
		}

		if(!_finite(x[0]))
		{
			printf("barrier x\n");
			scanf("%d",temp);
		}

		if(!_finite(q[0]))
		{
			printf("barrier q\n");
			scanf("%d",temp);
		}*/
		//Find bounding box
		vect2d minpoint, maxpoint;
		int i0 = indchange[0][(int)boundaryorder[0][0][0]];

		/*if(i0 == -1)
		{
			int ori = (int)boundaryorder[0][0][0];
			minpoint[0] = this->full_mesh[2 * ori];
			minpoint[1] = this->full_mesh[2 * ori+1];
			maxpoint[0] =this->full_mesh[2 * ori];
			maxpoint[1] = this->full_mesh[2 * ori+1];
		}
		else*/
		{
			minpoint[0] = min(x[2 * i0], x[2 * i0] + sA*q[2 * i0]);
			minpoint[1] = min(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]);
			maxpoint[0] = max(x[2 * i0], x[2 * i0] + sA*q[2 * i0]);
			maxpoint[1] = max(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]);
		}
		
		for (int j = 0; j < boundaryorder[0].size(); j++)
		{
			for (int i = 0; i < boundaryorder[0][j].size(); i++)
			{
				i0 = indchange[0][boundaryorder[0][j][i]];

				/*if(i0 == -1)
				{
					int ori = indchange[0][boundaryorder[0][j][i]];
					minpoint[0] = min(minpoint[0], full_mesh[2*ori] );
					minpoint[1] = min(minpoint[1], full_mesh[2*ori+1] );

					maxpoint[0] = max(maxpoint[0], full_mesh[2*ori] );
					maxpoint[1] = max(maxpoint[1], full_mesh[2*ori+1] );
				}
				else*/
				{
					minpoint[0] = min(minpoint[0], min(x[2 * i0], x[2 * i0] + sA*q[2 *i0]));
					minpoint[1] = min(minpoint[1], min(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]));

					maxpoint[0] = max(maxpoint[0], max(x[2 * i0], x[2 * i0] + sA*q[2 * i0]));
					maxpoint[1] = max(maxpoint[1], max(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]));
				}
			}
		}

		//Generate Tree
		tree.clear();
		//tree.setBoundingBox ( BoundingBox ( minpoint, maxpoint ) );

		tree.setBoundingBox (  minpoint, maxpoint );

		int i1;
		for (int j = 0; j < boundaryorder[0].size(); j++)
		{
			for (int i = 0; i < boundaryorder[0][j].size(); i++)
			{
				i0 = indchange[0][boundaryorder[0][j][i]];
				int ip1 = (i + 1) % boundaryorder[0][j].size();
				i1 = indchange[0][boundaryorder[0][j][ip1]];
				vect2d p1;
				vect2d p2;

				/*if(i0 == -1)
				{
					int ori = boundaryorder[0][j][i];
					p1[0] = min(full_mesh[ori*2],
					min(x[2 * i1], x[2 * i1] + sA*q[2 * i1])
					);

					p1[1] = min(full_mesh[ori*2+1],
						min(x[2 * i1 + 1], x[2 * i1 + 1] + sA*q[2 * i1 + 1])
						);

				
					p2[0] = max(full_mesh[ori*2],
						max(x[2 * i1], x[2 *i1] + sA*q[2 * i1])
						);

					p2[1] = max(full_mesh[ori*2],
						max(x[2 * i1 + 1], x[2 * i1 + 1] + sA*q[2 * i1 + 1])
						);
				}
				else if(i1 == -1)
				{
					int ori = boundaryorder[0][j][ip1];
					p1[0] = min(min(x[2 * i0], x[2 * i0] + sA*q[2 * i0]),
						full_mesh[ori*2]
						);

					p1[1] = min(min(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]),
						full_mesh[ori*2+1]
						);

				
					p2[0] = max(max(x[2 * i0], x[2 * i0] + sA*q[2 * i0]),
						full_mesh[ori*2]
						);

					p2[1] = max(max(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]),
						full_mesh[ori*2+1]
						);
				}
				else if(i0 == -1 && i1 == -1)
				{
					int ori = boundaryorder[0][j][i];
					int ori2 = boundaryorder[0][j][ip1];
					p1[0] = min(full_mesh[ori*2],
					full_mesh[ori2*2]
					);

					p1[1] = min(full_mesh[ori*2+1],
						full_mesh[ori2*2+1]
						);

				
					p2[0] = max(full_mesh[ori*2],
						full_mesh[ori2*2]
						);

					p2[1] = max(full_mesh[ori*2],
						full_mesh[ori2*2+1]
						);
				}
				else*/
				{
					p1[0] = min(min(x[2 * i0], x[2 * i0] + sA*q[2 * i0]),
						min(x[2 * i1], x[2 * i1] + sA*q[2 * i1])
						);

					p1[1] = min(min(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]),
						min(x[2 * i1 + 1], x[2 * i1 + 1] + sA*q[2 * i1 + 1])
						);

				
					p2[0] = max(max(x[2 * i0], x[2 * i0] + sA*q[2 * i0]),
						max(x[2 * i1], x[2 *i1] + sA*q[2 * i1])
						);

					p2[1] = max(max(x[2 * i0 + 1], x[2 * i0 + 1] + sA*q[2 * i0 + 1]),
						max(x[2 * i1 + 1], x[2 * i1 + 1] + sA*q[2 * i1 + 1])
						);
				}

				tree.insert(BoundingBox(p1, p2));
			}
		}
		//		tree.finishedInserting ( );

#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_sA_intervalsetup += endT - startT;
#endif

#ifdef DO_FULL_TIMINGS
		startT = get_time();
#endif
	double privSA[NUM_PROCS];

	//Find intersections and test
		for (int j = 0; j < boundaryorder[0].size(); j++)
		{
			for (int i = 0; i < NUM_PROCS; i++)
			{
				privSA[i] = sA;
			}

#pragma omp parallel for
			for (int i = 0; i < boundaryorder[0][j].size(); i++)
			{
				int_pos[ omp_get_thread_num() ].clear();
				int ind = indchange[0][boundaryorder[0][j][i]];

				vect2d p1;
				vect2d p2;

				/*if(ind == -1)
				{
					int ori = boundaryorder[0][j][i];
					p1[0] = full_mesh[2*ori] - BBOX_INTER_THRES ;
					p1[1] = full_mesh[2*ori+1] - BBOX_INTER_THRES;
				
					p2[0] = full_mesh[2*ori] + BBOX_INTER_THRES;
					p2[1] = full_mesh[2*ori+1] + BBOX_INTER_THRES;
				}
				else*/
				{
					p1[0] = min(x[2 * ind], x[2 * ind] + sA*q[2 * ind]);
					p1[1] = min(x[2 * ind + 1], x[2 * ind + 1] + sA*q[2 * ind + 1]);
				
					p2[0] = max(x[2 * ind], x[2 * ind] + sA*q[2 * ind]);
					p2[1] = max(x[2 * ind + 1], x[2 * ind + 1] + sA*q[2 * ind + 1]);
				}
				tree.findIntersect(BoundingBox(p1, p2), int_pos[omp_get_thread_num()]);
				//vect2d p1, p2;
				//p1[0] = x[2*ind];
				//p1[1] = x[2*ind+1];
				//p2[0] = x[2*ind] + sA*q[2*ind];
				//p2[1] = x[2*ind+1] + sA*q[2*ind+1];
				//tree.findIntersectLineSegment ( p1, p2, int_pos );

				for (int jj = 0; jj < int_pos[omp_get_thread_num()].length(); jj++)
				{
					int indt = int_pos[omp_get_thread_num()].data[jj];
					
					int bi = 0;
					for (int t = 0; t < boundaryorder[0].size(); t++)
					{
						if (indt >= boundaryorder[0][t].size())
							indt -= boundaryorder[0][t].size();
						else break;
						bi++;
					}

					i0 = indchange[0][boundaryorder[0][bi][indt]];
					i1 = indchange[0][boundaryorder[0][bi][(indt + 1) % boundaryorder[0][bi].size()]];

					if(i0 == -1 && i1 == -1 && ind == -1)
						continue;


					/*change_sA(sA, x, q,
						i0,
						i1,
						ind
						);*/

					change_sA(privSA[omp_get_thread_num()], x, q,
						boundaryorder[0][bi][indt],
						boundaryorder[0][bi][(indt + 1) % boundaryorder[0][bi].size()],
						boundaryorder[0][j][i]);
				}
			}
			sA = privSA[0];

			for (int i = 1; i < omp_get_max_threads(); i++)
			{
				sA = min(sA, privSA[i]);
			}

		}
/*
	change_sA( sA, x, q, 
					boundaryorder[i],  
					boundaryorder[ (i + 1)%boundaryorder.size()],
					boundaryorder[ mins_interval[j].first ]
					);*/

#ifdef DO_FULL_TIMINGS
	endT = get_time();
	time_sA_bound += endT- startT;
#endif


	return sA;
}
#else
double BarrierMethod::max_param( double *f, double *x, double *gg,  double *q, int n )
{
	//Get intervals
#ifdef DO_FULL_TIMINGS
		double startT = get_time();
#endif
	double sA = currentsA;

#ifdef SCOTTS_CODE
	for(int i = 0; i < boundaryorder.size(); i++)
	{
		double t;

		t = min( min(x[2*boundaryorder[i]], x[2*boundaryorder[i]] + sA*q[2*boundaryorder[i]]),
			     min(x[2*boundaryorder[(i+1)%boundaryorder.size()]], x[2*boundaryorder[(i+1)%boundaryorder.size()]] + sA*q[2*boundaryorder[(i+1)%boundaryorder.size()]])
			);
		std::pair<int, double> pt(i,t);
		mins_interval[i] = pt;

		t = max( max(x[2*boundaryorder[i]], x[2*boundaryorder[i]] + sA*q[2*boundaryorder[i]]),
			     max(x[2*boundaryorder[(i+1)%boundaryorder.size()]], x[2*boundaryorder[(i+1)%boundaryorder.size()]] + sA*q[2*boundaryorder[(i+1)%boundaryorder.size()]])
				 );

		maxs_interval[i] = t;

		t = min( x[2*boundaryorder[i]], x[2*boundaryorder[i]] + sA*q[2*boundaryorder[i]]);
		std::pair<int, double> ptt(i,t);
		pmins_interval[i] = ptt;

		t = max(x[2*boundaryorder[i]], x[2*boundaryorder[i]] + sA*q[2*boundaryorder[i]]);

		pmaxs_interval[i] = t;
	}

#ifdef DO_FULL_TIMINGS
		double time_sA_intervalsetup += endT - startT;
#endif

#ifdef DO_FULL_TIMINGS
		startT = get_time();
#endif
	//sort
	sort( mins_interval.begin(), mins_interval.end(), interval_sort_func );
	sort( pmins_interval.begin(), pmins_interval.end(), interval_sort_func );

#ifdef DO_FULL_TIMINGS
	endT = get_time();
	time_sA_sort+= endT- startT;
#endif

	int ind = 0;
	vect2d sv, fv, pv, start, finish, pp;
	std::pair<double,double> ttt;

#ifdef DO_FULL_TIMINGS
		startT = get_time();
#endif


	//bool first = false;
	for(int i = 0; i < boundaryorder.size(); i++)
	{
		bool first = false;
		for(int j = ind; j < boundaryorder.size(); j++)
		{

			if( (pmins_interval[j].second <= maxs_interval[mins_interval[i].first] && pmins_interval[j].second >= mins_interval[i].second)
				     ||
				     (pmaxs_interval[pmins_interval[j].first] <= maxs_interval[mins_interval[i].first] && pmaxs_interval[pmins_interval[j].first] >= mins_interval[i].second)
				   ) //point in interval
			{

				if(!first)
				{
					ind = j;
					first = true;
				}
				

//#ifdef TTTT
//				int ii = boundaryorder[ mins_interval[i].first ];
//				int iin = boundaryorder[ (mins_interval[i].first+1)%g.boundaryorder.size() ];
//				int jj = boundaryorder[ pmins_interval[j].first ];
//				change_sA( sA, x, q, 
//					ii,  
//					iin,
//					jj
//					);
//#else

				int ii = boundaryorder[ mins_interval[i].first ];
				int iin = boundaryorder[ (mins_interval[i].first+1)%boundaryorder.size() ];
				int jj = boundaryorder[ pmins_interval[j].first ];

				if(jj != ii && jj != iin)
				{
					start[0] = x[ii*2];
					start[1] = x[ii*2+1];

					finish[0] = x[iin*2];
					finish[1] = x[iin*2+1];

					pp[0] = x[jj*2];
					pp[1] = x[jj*2+1];

					sv[0] = q[ii*2];
					sv[1] = q[ii*2+1];

					fv[0] = q[iin*2];
					fv[1] = q[iin*2+1];

					pv[0] = q[jj*2];
					pv[1] = q[jj*2+1];

					ttt = findIntersect(start, sv, finish, fv, pp, pv, sA);

					if(ttt.first != -1)
					{
						if(ttt.first > 0 && ttt.first < sA)
							sA = ttt.first;
					}

				}
//#endif
			}
			else if( pmins_interval[j].second > maxs_interval[mins_interval[i].first] ) //point past interval
			{
				j = boundaryorder.size();
			}
		}
	}
	#ifdef DO_FULL_TIMINGS
	endT = get_time();
	time_sA_bound += endT- startT;
#endif
#else

	for(int i = 0; i < boundaryorder.size(); i++)
	{
		double t;

		t = min( min(x[2*boundaryorder[i]], x[2*boundaryorder[i]] + sA*q[2*boundaryorder[i]]),
			     min(x[2*boundaryorder[(i+1)%boundaryorder.size()]], x[2*boundaryorder[(i+1)%boundaryorder.size()]] + sA*q[2*boundaryorder[(i+1)%boundaryorder.size()]])
			);
		std::pair<int, double> pt(i,t);
		mins_interval[i] = pt;

		t = max( max(x[2*boundaryorder[i]], x[2*boundaryorder[i]] + sA*q[2*boundaryorder[i]]),
			     max(x[2*boundaryorder[(i+1)%boundaryorder.size()]], x[2*boundaryorder[(i+1)%boundaryorder.size()]] + sA*q[2*boundaryorder[(i+1)%boundaryorder.size()]])
				 );

		maxs_interval[i] = t;
	}

#ifdef DO_FULL_TIMINGS
		double endT = get_time();
		time_sA_intervalsetup += endT - startT;
#endif
#ifdef DO_FULL_TIMINGS
		startT = get_time();
#endif
	//sort
	sort( mins_interval.begin(), mins_interval.end(), interval_sort_func );
//	sort( pmins_interval.begin(), pmins_interval.end(), interval_sort_func );
#ifdef DO_FULL_TIMINGS
	endT = get_time();
	time_sA_sort+= endT- startT;
#endif

#ifdef DO_FULL_TIMINGS
		startT = get_time();
#endif
	//Check collisions
	for(int i = 0; i < boundaryorder.size(); i++)
	{
//		for(int j = i+1; j < g.boundaryorder.size(); j++)
		for(int j = 0; j < boundaryorder.size() && mins_interval [ j ].second <= maxs_interval[mins_interval[i].first]; j++)
		{
			if ( //mins_interval[j].second > maxs_interval[mins_interval[i].first] ||
				maxs_interval[mins_interval[j].first] < mins_interval[i].second )
			{
			}
			else
//			if( mins_interval[j].second <=  maxs_interval[mins_interval[i].first] )
			//if( x[j] <  maxs_interval[mins_interval[i].first] )
			{

				//Checking both points of interval... A little redundant
				change_sA( sA, x, q, 
					boundaryorder[mins_interval[i].first],  
					boundaryorder[ (mins_interval[(i)].first + 1)%boundaryorder.size()],
					boundaryorder[ mins_interval[j].first ]
					);
/*
				change_sA( sA, x, q, 
					boundaryorder[mins_interval[i].first],  
					boundaryorder[ (mins_interval[(i)].first + 1)%boundaryorder.size()],
					boundaryorder[ (mins_interval[(j)].first + 1)%boundaryorder.size() ]
					);
*/
					if(sA < .00001)
					{
					//	printf("sA = %.10g\n", sA);
					}
			}
//			else
//			{
			//	j = g.boundaryorder.size();
//			}
		}
	}
	#ifdef DO_FULL_TIMINGS
	endT = get_time();
	time_sA_bound += endT- startT;
#endif
#endif

	return sA;
}
#endif
void BarrierMethod::calc_interior_gradient(double *pos,  double *g )
{
}
void BarrierMethod::display_error( double *e, double *pos )
{
	ParamMethod::display_error(e,pos);
}


double BarrierMethod::calc_seamless_error(double a, double b)
{
	return seamless_exp * ( 4*a*(a-b)*b*(a+b)  + a*a*a*a - 6*a*a*b*b + b*b*b*b  );
}

double BarrierMethod::calc_seamless_error(double u1, double v1, double u2, double v2)
{
	return seamless_exp * ( (((u1)*(u1)*(u1)*(u1)*(u1)*(u1)*(u1)*(u1))+((-8.0)*((u1)*(u1)*(u1)*(u1)*(u1)*(u1)*(u1))*(u2))+((u2)*(u2)*(u2)*(u2)*(u2)*(u2)*(u2)*(u2))+((4.0)*((u1)*(u1)*(u1)*(u1)*(u1)*(u1))*(((7.0)*((u2)*(u2)))+(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2))))))+((-8.0)*((u1)*(u1)*(u1)*(u1)*(u1))*(u2)*(((7.0)*((u2)*(u2)))+((3.0)*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))))+((-8.0)*((u1)*(u1)*(u1))*(u2)*(((7.0)*((u2)*(u2)*(u2)*(u2)))+((10.0)*((u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+((-3.0)*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))))+((2.0)*((u1)*(u1)*(u1)*(u1))*(((35.0)*((u2)*(u2)*(u2)*(u2)))+((30.0)*((u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+((-3.0)*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))))+((4.0)*((u1)*(u1))*(((7.0)*((u2)*(u2)*(u2)*(u2)*(u2)*(u2)))+((15.0)*((u2)*(u2)*(u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+((-9.0)*((u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+((-1.0)*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))))+((4.0)*((u2)*(u2)*(u2)*(u2)*(u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+((-6.0)*((u2)*(u2)*(u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+((-4.0)*((u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2))))+((-8.0)*(u1)*(u2)*(((u2)*(u2)*(u2)*(u2))+((4.0)*((u2)*(u2))*(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))+(((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))*((v1)+((-1.0)*(v2)))))*((u2)+(v1)+((-1.0)*(v2)))*((u2)+((-1.0)*(v1))+(v2)))) );
}
/*************************************************/
//Functions
/*************************************************/	



/*******************************************************/
//HELPER FUNCS
/*******************************************************/
void BarrierMethod::sharpen()
{
	func_exp += 2;
	printf("Sharpen Function to = %d\n", func_exp);
}

void BarrierMethod::sharpen_alpha()
{
	alpha_val += ALPHA_INC;
	printf("sharpern alpha = %f\n", alpha_val);
}

void BarrierMethod::reduce_alpha()
{
	alpha_val -= ALPHA_INC;
	printf("sharpern alpha = %f\n", alpha_val);
}


void BarrierMethod::fill_point_grid(const double *pos)
{
#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif
	//compute bbox
	pg.clear();
	int i0 = indchange[0][boundaryorder[0][0][0]];
	vect2d minpoint, maxpoint;
	/*minpoint[0] = i0 == -1 ? full_mesh[2*boundaryorder[0][0][0]] : pos[2 * i0]; 
	minpoint[1] = i0 == -1 ? full_mesh[2*boundaryorder[0][0][0]+1] :pos[2 * i0 + 1];
	maxpoint[0] = i0 == -1 ? full_mesh[2*boundaryorder[0][0][0]] : pos[2 * i0]; 
	maxpoint[1] = i0 == -1 ? full_mesh[2*boundaryorder[0][0][0]+1] :pos[2 * i0+ 1];
*/
	minpoint[0] = pos[2 * i0]; 
	minpoint[1] = pos[2 * i0 + 1];
	maxpoint[0] =  pos[2 * i0]; 
	maxpoint[1] = pos[2 * i0+ 1];
	
	for (int j = 0; j < boundaryorder[0].size(); j++)
	{
		for (int i = 0; i < boundaryorder[0][j].size(); i++)
		{
			
			i0 = indchange[0][boundaryorder[0][j][i]];

			/*minpoint[0] = min((double)minpoint[0], i0 == -1 ?  full_mesh[2*boundaryorder[0][j][i]] : (double)pos[2 * i0]);
			minpoint[1] = min((double)minpoint[1], i0 == -1 ?  full_mesh[2*boundaryorder[0][j][i]+1] :(double)pos[2 * i0 + 1]);
			maxpoint[0] = max((double)maxpoint[0], i0 == -1 ?  full_mesh[2*boundaryorder[0][j][i]] :(double)pos[2 * i0]);
			maxpoint[1] = max((double)maxpoint[1], i0 == -1 ?  full_mesh[2*boundaryorder[0][j][i]+1] :(double)pos[2 * i0 + 1]);
*/
			minpoint[0] = min((double)minpoint[0],  (double)pos[2 * i0]);
			minpoint[1] = min((double)minpoint[1], (double)pos[2 * i0 + 1]);
			maxpoint[0] = max((double)maxpoint[0], (double)pos[2 * i0]);
			maxpoint[1] = max((double)maxpoint[1], (double)pos[2 * i0 + 1]);
		}
	}
	pg.setBoundingBox(minpoint, maxpoint);

	//Fill grid
	for (int j = 0; j < boundaryorder[0].size(); j++)
	{
		for (int i = 0; i < boundaryorder[0][j].size(); i++)
		{
			i0 = indchange[0][boundaryorder[0][j][i]];
			vect2d t; 
			/*t[0] = i0 == -1 ? full_mesh[2*boundaryorder[0][j][i]] : pos[2 * i0]; 
			t[1] = i0 == -1 ? full_mesh[2*boundaryorder[0][j][i]+1] : pos[2 * i0 + 1];
*/
			t[0] =  pos[2 * i0]; 
			t[1] =  pos[2 * i0 + 1];
			
			int ex = 0;
			for (int jj = 0; jj < j; jj++)
				ex += boundaryorder[0][jj].size();
			
			pg.addPoint(t, i + ex);
		}
	}

#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_funcgrad_intervalsetup += endT - startT;
#endif
}

double point_to_line_distance(vect2d p, vect2d a, vect2d b)
{
	double t = ((p - a).dot(b-a)) / (a - b).length2();

  if (t < 0.0) 
  {
	  return (p - a).length();      
  }
  else if (t > 1.0) 
  {
	  return (p- b).length();
  }
  return (p - (a + t * (b - a)) ).length();
}

/*******************************************************/
//FUNCTION EVALS
/*******************************************************/
double BarrierMethod::calc_error(double* pos, double step, double* grad)
{
	double fx = ParamMethod::calc_error(pos, step, grad);
#ifdef TURN_OFF_BOUNDARY
#else
	//Get global function error

	if(g.do_boundary)
	{
	#ifdef SEAMLESS_WITHOUT_BOUND
	#else
		if(boundaryorder[0].size() > 0)
		{
	#ifdef GLOBAL_BOUNDARY_SUPPORT
	#else
		fill_point_grid(pos);
	#endif
		double t = fx;
		fx += perform_full_global_barrier_func(fx, pos);//perform_global_fold_barrier_func(fx, pos);
	
		if (fx != t)
		{
			//printf("%bar fx = %.10g\n", fx);
		}
		}
		/*else
		{
			printf("0");
		}*/
	#endif
	}
#endif
	
#ifdef PERFORM_SEAMLESS

	int ind = 0;
	int ind2= 0;
	int next = 0;
	int next2 = 0;
	double seamf = 0;
	double largest = 0;
	int indlarge = 0;
	double u1,v1,u2,v2,u3,v3,u4,v4;
	for(int i = 0; i < this->match_boundary.size(); i++)
	{
		next = match_boundary[i].u1;
		next2 = match_boundary[i].v1;

		u1 = pos[2*next];
		v1 = pos[2*next+1];
		u2 = pos[2*next2];
		v2 = pos[2*next2 + 1];

		u3 = pos[2*match_boundary[i].u2];
		v3 = pos[2*match_boundary[i].u2+1];
		u4 = pos[2*match_boundary[i].v2];
		v4 = pos[2*match_boundary[i].v2+1];
		//printf("%f,%f,%f,%f\n",pos[2*next], pos[2*next2] ,  pos[2*next+1], pos[2*next2+1]);
		match_boundary[i].p1 = ((((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
		match_boundary[i].p2 = ((((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
		match_boundary[i].p3 = ((((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
		match_boundary[i].p4 = ((((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));

		/*match_boundary[i].p1 = ((abs((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(abs((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
		match_boundary[i].p2 = ((abs((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(abs((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
		match_boundary[i].p3 = ((abs((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(abs((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
		match_boundary[i].p4 = ((abs((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(abs((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));*/
		//printf("%.10g\n",match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4);
		seamf += seamless_exp * match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4;

		if(match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4 > largest)
		{
			largest = match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4;
			indlarge = i;
		}
	}

	u1 = pos[2*match_boundary[indlarge].u1];
	v1 = pos[2*match_boundary[indlarge].u1+1];
	u2 = pos[2*match_boundary[indlarge].v1];
	v2 = pos[2*match_boundary[indlarge].v1 + 1];

	u3 = pos[2*match_boundary[indlarge].u2];
	v3 = pos[2*match_boundary[indlarge].u2+1];
	u4 = pos[2*match_boundary[indlarge].v2];
	v4 = pos[2*match_boundary[indlarge].v2+1];
	vect2d temp;
	temp[0] = u1 - u2;
	temp[1] = v1- v2;

	vect2d temp2;
	temp2[0] = u3 - u4;
	temp2[1] = v3- v4;
	printf("\r %f,   %f,  %f,  %f, %d", seamf, largest, temp.length(), temp2.length(), indlarge);

#ifdef SEAMLESS_WITHOUT_BOUND
	return seamf;
#else
	return fx + seamf;
#endif

#else
	return fx;
	#endif
}

double BarrierMethod::perform_full_global_barrier_func(double ff, const double *pos)
{
	double retfx = 0;


#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif

	
	//Find possible cases
	for (int j = 0; j < boundaryorder[0].size(); j++)
	{
		double fx = 0;
#pragma omp parallel for reduction (+:fx)
		for (int i = 0; i < boundaryorder[0][j].size(); i++)
		{
			int t_id = omp_get_thread_num();
			int i0, i1;
			i0 = indchange[0][boundaryorder[0][j][i]]; 
			i1 = indchange[0][boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()]];
		
			//printf("%d,",t_id);

#ifdef GLOBAL_BOUNDARY_SUPPORT
			for(int jj = 0; jj < boundaryorder[0][j].size(); jj++)
			{
				global_barrier_eval(fx, pos, 
					boundaryorder[0][j][i],  
					boundaryorder[0][j][(i + 1) % boundaryorder[0].size()],
					boundaryorder[0][j][jj]
					);
			}
#else
			int_pos[t_id].clear();
			//ResizeArray<int> int_posT;
			/*vect2d t1;
			t1[0] = i0 == -1 ?  full_mesh[2*boundaryorder[0][j][i]] : pos[2 * i0]; 
			t1[1] = i0 == -1 ?  full_mesh[2*boundaryorder[0][j][i]+1] : pos[2 * i0 + 1];
			vect2d t2;
			t2[0] = i1 == -1 ? full_mesh[2*boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()]] : pos[2 * i1];
			t2[1] = i1 == -1 ? full_mesh[2*boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()]+1] : pos[2 * i1 + 1];
*/
			vect2d t1;
			t1[0] =  pos[2 * i0]; 
			t1[1] =  pos[2 * i0 + 1];
			vect2d t2;
			t2[0] =  pos[2 * i1];
			t2[1] =  pos[2 * i1 + 1];

			pg.findPoints(int_pos[t_id], t1, t2, CONST_CLAMP);
			//pg.findPoints ( int_posT, t1, t2, CONST_CLAMP );

			for (int jj = 0; jj < int_pos[t_id].length(); jj++)
				//for(int j = 0; j < int_posT.length(); j++)
			{
				int ind = int_pos[t_id].data[jj];
				int bi = 0;
				for (int t = 0; t < boundaryorder[0].size(); t++)
				{
					if (ind >= boundaryorder[0][t].size())
						ind -= boundaryorder[0][t].size();
					else break;
					bi++;
				}

				global_barrier_eval(fx, pos,
					boundaryorder[0][j][i],
					boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()],
					//boundaryorder[0][int_posT.data[j]]
					boundaryorder[0][bi][ind]
					);
			}
#endif
		}

		retfx += fx;
	}


#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_func_bound += endT - startT;
#endif

#ifdef PERFORM_SEAMLESS

	//int ind = 0;
	//int ind2= 0;
	//int next = 0;
	//int next2 = 0;
	//double seamf = 0;
	//double largest = 0;
	//int indlarge = 0;
	//double u1,v1,u2,v2,u3,v3,u4,v4;
	//for(int i = 0; i < this->match_boundary.size(); i++)
	//{
	//	next = match_boundary[i].u1;
	//	next2 = match_boundary[i].v1;

	//	u1 = pos[2*next];
	//	v1 = pos[2*next+1];
	//	u2 = pos[2*next2];
	//	v2 = pos[2*next2 + 1];

	//	u3 = pos[2*match_boundary[i].u2];
	//	v3 = pos[2*match_boundary[i].u2+1];
	//	u4 = pos[2*match_boundary[i].v2];
	//	v4 = pos[2*match_boundary[i].v2+1];
	//	//printf("%f,%f,%f,%f\n",pos[2*next], pos[2*next2] ,  pos[2*next+1], pos[2*next2+1]);
	//	match_boundary[i].p1 = ((((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
	//	match_boundary[i].p2 = ((((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
	//	match_boundary[i].p3 = ((((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
	//	match_boundary[i].p4 = ((((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));

	//	/*match_boundary[i].p1 = ((abs((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(abs((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
	//	match_boundary[i].p2 = ((abs((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(abs((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
	//	match_boundary[i].p3 = ((abs((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(abs((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
	//	match_boundary[i].p4 = ((abs((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(abs((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));*/
	//	//printf("%.10g\n",match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4);
	//	seamf += seamless_exp * match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4;

	//	if(match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4 > largest)
	//	{
	//		largest = match_boundary[i].p1*match_boundary[i].p2*match_boundary[i].p3*match_boundary[i].p4;
	//		indlarge = i;
	//	}
	//}

	//u1 = pos[2*match_boundary[indlarge].u1];
	//v1 = pos[2*match_boundary[indlarge].u1+1];
	//u2 = pos[2*match_boundary[indlarge].v1];
	//v2 = pos[2*match_boundary[indlarge].v1 + 1];

	//u3 = pos[2*match_boundary[indlarge].u2];
	//v3 = pos[2*match_boundary[indlarge].u2+1];
	//u4 = pos[2*match_boundary[indlarge].v2];
	//v4 = pos[2*match_boundary[indlarge].v2+1];
	//vect2d temp;
	//temp[0] = u1 - u2;
	//temp[1] = v1- v2;

	//vect2d temp2;
	//temp2[0] = u3 - u4;
	//temp2[1] = v3- v4;
	//printf("\r %f,   %f,  %f,  %f, %d", seamf, largest, temp.length(), temp2.length(), indlarge);
	//return retfx + ff + seamf;
	return retfx + ff;
#else
	return retfx + ff;
#endif


	

}

void BarrierMethod::set_seamless_param(double t)
{
	seamless_exp = t;
}
//Setup single vertex call
#ifdef GLOBAL_BOUNDARY_SUPPORT
void BarrierMethod::global_barrier_eval(double &func, const double *pos, int i, int ip1,  int j )
{
	

	if(i != j && ip1 != j )
	{
		vect2d p0;
	p0[0] = pos[i*2];
	 p0[1] = pos[i*2 + 1];

	/*	double u0 = pos[i*2];
	double v0 = pos[i*2 + 1];*/
	 vect2d p1;
	p1[0]= pos[ip1*2];
	p1[1] = pos[ip1*2 + 1];
	/*double u1 = pos[ip1*2];
	double v1 = pos[ip1*2 + 1];*/
	vect2d p;
	p[0] = pos[j*2];
	p[1] = pos[j*2 + 1];
	/*double x= pos[j*2];
	double y = pos[j*2 + 1];*/

	double alpha = alpha_val;

#ifdef ALPHA_TIMES_X
	func += alpha / ((p-p0).length() +  (p-p1).length() - (p0- p1).length());
#else
	func += pow(1.0 / ((p-p0).length() +  (p-p1).length() - (p0- p1).length()), alpha);
#endif
	}
}
#else
void BarrierMethod::global_barrier_eval(double &func, const double *pos, int i, int ip1,  int j )
{
	int i0 = i; // indchange[0][i];
	int i1 = ip1; // indchange[0][ip1];
	int i2 = j; // indchange[0][j];

//	if(i0 == -1 && i1 == -1 && i2 == -1)
//		return;

	if(i != j && ip1 != j )
	{
		/*vect2d a;
		a[0] = i0 == -1 ? full_mesh[i*2] : pos[i0*2];
		a[1] = i0 == -1 ? full_mesh[i*2+1] : pos[i0*2+1];

		vect2d b;
		b[0] = i1 == -1 ? full_mesh[ip1*2] : pos[i1*2];
		b[1] = i1 == -1 ? full_mesh[ip1*2+1] : pos[i1*2+1];

		vect2d x;
		x[0] = i2 == -1 ? full_mesh[j*2] : pos[i2*2];
		x[1] = i2 == -1 ? full_mesh[j*2+1] : pos[i2*2+1];
*/
		vect2d a;
		a[0] =  pos[i0*2];
		a[1] =  pos[i0*2+1];

		vect2d b;
		b[0] =  pos[i1*2];
		b[1] =  pos[i1*2+1];

		vect2d x;
		x[0] =  pos[i2*2];
		x[1] =  pos[i2*2+1];

		func += global_barrier_eval(a,b,x);
	}
}
#endif
//Not used... gradient step pushed ot earlier in calls to remove need to update position constantly
//current position
//segment i -> ip1
//vertex j
void BarrierMethod::global_barrier_eval(double &func, const double *pos, const double *grad, double step, int i, int ip1,  int j )
{
	int i0 = indchange[0][i];
	int i1 = indchange[0][ip1];
	int i2 = indchange[0][j];

	if(i != j && ip1 != j )
	{
		/*vect2d a;
		a[0] = i0 == -1 ? full_mesh[i*2] : pos[i0*2] - step*grad[i0*2];
		a[1] = i0 == -1 ? full_mesh[i*2+1] :  pos[i0*2+1]  - step*grad[i0*2+1];

		vect2d b;
		b[0] = i1 == -1 ? full_mesh[ip1*2] : pos[i1*2]  - step*grad[i1*2];
		b[1] = i1 == -1 ? full_mesh[ip1*2+1] : pos[i1*2+1] - step*grad[i1*2+1];

		vect2d x;
		x[0] = i2 == -1 ? full_mesh[j*2] : pos[i2*2] - step*grad[i2*2];
		x[1] = i2 == -1 ? full_mesh[j*2+1] : pos[i2*2+1] - step*grad[i2*2+1];
*/
		vect2d a;
		a[0] =  pos[i0*2] - step*grad[i0*2];
		a[1] =   pos[i0*2+1]  - step*grad[i0*2+1];

		vect2d b;
		b[0] =  pos[i1*2]  - step*grad[i1*2];
		b[1] =  pos[i1*2+1] - step*grad[i1*2+1];

		vect2d x;
		x[0] =  pos[i2*2] - step*grad[i2*2];
		x[1] =  pos[i2*2+1] - step*grad[i2*2+1];

		func += global_barrier_eval(a,b,x);
	}
}

//Calculates actual boundary function for single vertex vs linesegment
double BarrierMethod::global_barrier_eval(vect2d a, vect2d b, vect2d x)
{
	double clamp_dist = CONST_CLAMP;
	double eps_slope = EPSION_FUNC;

	double ll = ( a - b ).length2 ( );
	double l = sqrt ( ll );
	vect2d ortho;
	ortho [ 0 ] = (b [ 1 ] - a [ 1 ]) / l;
	ortho [ 1 ] = (a [ 0 ] - b [ 0 ]) / l;
	double orthoDist = ( x - b ).dot( ortho );
	double dist = abs ( orthoDist );
	if ( dist > clamp_dist )			
		return 0;

	double t = ((x - a).dot(b-a)) / ll;
	if (t < 0.0) 
	{
		dist = (x - a).length();      
	}
	else if (t > 1.0) 
	{
		dist = (x- b).length();
	}
/*	else
	{
		dist = ( x - (a + t * (b - a)) ).length();
	}*/

	return pow( max( (double)0, (double)(1.0f / ((eps_slope/clamp_dist) * dist)  - 1.0f / eps_slope )), func_exp);

}


/*******************************************************/
//GRADIENT EVALS
/*******************************************************/

void BarrierMethod::calc_full_gradient( double *pos,  double *gg )
{
	ParamMethod::calc_full_gradient(pos, gg);
#ifdef TURN_OFF_BOUNDARY

#else
	if(g.do_boundary)
	{
	#ifdef SEAMLESS_WITHOUT_BOUND
	#else
		if(boundaryorder[0].size() > 0)
		{
	#ifdef GLOBAL_BOUNDARY_SUPPORT
	#else
		fill_point_grid(pos);
	#endif

		calc_boundary_gradient(pos, gg);
		}
	#endif
	}

#endif
	

#ifdef PERFORM_SEAMLESS

	int ind = 0;
	int ind2= 0;
	int next = 0;
	int next2 = 0;
	double seamf = 0;

	double u1,u2,v1,v2,u3,v3,u4,v4;

	for(int i = 0; i < this->match_boundary.size(); i++)
	{
		next = match_boundary[i].u1;
		next2 = match_boundary[i].v1;

		u1 = pos[2*next];
		v1 = pos[2*next+1];
		u2 = pos[2*next2];
		v2 = pos[2*next2+1];

		u3 = pos[2*match_boundary[i].u2];
		v3 = pos[2*match_boundary[i].u2+1];
		u4 = pos[2*match_boundary[i].v2];
		v4 = pos[2*match_boundary[i].v2+1];

		gg[ 2*next ] += seamless_exp*( (((2.0)*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((2.0)*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))) ));

		gg[ 2*next2 ] += seamless_exp*( (((-2.0)*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((2.0)*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((2.0)*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((2.0)*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));
		
		gg[ 2*next+1] +=seamless_exp*( (((-2.0)*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((-2.0)*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((-2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((-2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))));

		gg[ 2*next2+1] += seamless_exp*( (((2.0)*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((-2.0)*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((-2.0)*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((-2.0)*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));

		gg[ 2*match_boundary[i].u2 ] += seamless_exp*( (((2.0)*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((2.0)*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((2.0)*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((-2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))));

		gg[2*match_boundary[i].u2+1 ] += seamless_exp* ((((2.0)*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((2.0)*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((-2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((-2.0)*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));
		
		gg[2*match_boundary[i].v2] +=seamless_exp* ((((-2.0)*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((-2.0)*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((-2.0)*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))));

		gg[ 2*match_boundary[i].v2+1] += seamless_exp* ((((-2.0)*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * ((-2.0)*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * ((2.0)*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * ((2.0)*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));



		/*gg[ 2*next ] += seamless_exp*( ((abs((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))) ));

		gg[ 2*next2 ] += seamless_exp*( ((abs((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));
		
		gg[ 2*next+1] +=seamless_exp*( ((abs((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))));

		gg[ 2*next2+1] += seamless_exp*( ((abs((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));

		gg[ 2*match_boundary[i].u2 ] += seamless_exp*( ((abs((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))));

		gg[2*match_boundary[i].u2+1 ] += seamless_exp* (((abs((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));
		
		gg[2*match_boundary[i].v2] +=seamless_exp* (((abs((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))));

		gg[ 2*match_boundary[i].v2+1] += seamless_exp* (((abs((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p2 * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * (abs((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))) * match_boundary[i].p3 * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * (abs((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))) * match_boundary[i].p4)+
			(match_boundary[i].p1 * match_boundary[i].p2 * match_boundary[i].p3 * (abs((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4)))));*/
	}


#endif
	//calc interior for each method
}

void BarrierMethod::calc_boundary_gradient(double* pos, double *g)
{
	//Calc boundary code here!!!
#ifdef DO_FULL_TIMINGS
	double startT = get_time();
#endif

	for (int j = 0; j < boundaryorder[0].size(); j++)
	{
		//Find possible cases
		for(int i = 0; i < boundaryorder[0][j].size(); i++)
		{
#ifdef GLOBAL_BOUNDARY_SUPPORT
			for(int jj = i + 1; jj < boundaryorder[0][j].size() - 1; jj++)
			{
				global_barrier_grad(g, pos,
					boundaryorder[0][j][jj],  
					boundaryorder[0][j][jj + 1],
					boundaryorder[0][j][i]
					);
			}
			if (i != 0 && i != boundaryorder[0][j].size() - 1)
			{
				global_barrier_grad(g, pos,
					boundaryorder[0][j][boundaryorder[0][j].size() - 1],
					boundaryorder[0][j][0],
					boundaryorder[0][j][i]
					);
			}
			for(int jj = 0; jj < i - 1; jj++)
			{
				global_barrier_grad(g, pos,
					boundaryorder[0][j][jj],
					boundaryorder[0][j][jj + 1],
					boundaryorder[0][j][i]
					);
			}
			//for(int j = 0; j < boundaryorder[0].size(); j++)
			//{
			//	global_barrier_grad(g, pos,
			//				boundaryorder[0][i],  
			//				boundaryorder[0][ (i + 1) % boundaryorder[0].size ( )],
			//				boundaryorder[0][j]
			//			);
			//}

#else
			int i0 = indchange[0][boundaryorder[0][j][i]];
			int i1 = indchange[0][boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()]];
			in_pos.clear();
			vect2d t1, t2;
			/*t1[0] = i0 == -1 ? full_mesh[2*boundaryorder[0][j][i]] : pos[2 * i0];
			t1[1] = i0 == -1 ? full_mesh[2*boundaryorder[0][j][i]+1] : pos[2 * i0 + 1];
			t2[0] = i1 == -1 ? full_mesh[2*boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()]] : pos[2 * i1];
			t2[1] = i1 == -1 ? full_mesh[2*boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()]+1] : pos[2 * i1 + 1];
*/
			t1[0] =  pos[2 * i0];
			t1[1] =  pos[2 * i0 + 1];
			t2[0] = pos[2 * i1];
			t2[1] =  pos[2 * i1 + 1];
			pg.findPoints(in_pos, t1, t2, CONST_CLAMP);

			for (int jj = 0; jj < in_pos.length(); jj++)
			{
				int indt = in_pos.data[jj];
				int bi = 0;
				for (int t = 0; t < boundaryorder[0].size(); t++)
				{
					if (indt >= boundaryorder[0][t].size())
						indt -= boundaryorder[0][t].size();
					else break;
					bi++;
				}


				global_barrier_grad(g, pos,
					boundaryorder[0][j][i],
					boundaryorder[0][j][(i + 1) % boundaryorder[0][j].size()],
					boundaryorder[0][bi][indt]
					);
/*
					global_barrier_grad(g, pos,
					i0,
					i1,
					indchange[0][boundaryorder[0][bi][indt]]
					);*/
			}
#endif
		}
	}

#ifdef DO_FULL_TIMINGS
	double endT = get_time();
	time_grad_bound += endT - startT;
#endif
}

//Calculates actual gradient for vertex vs lines segment
//answer in newgradient
//current position
//segment i -> ip1
//vertex j

#ifdef GLOBAL_BOUNDARY_SUPPORT
void BarrierMethod::global_barrier_grad(double *newgradient, const double *pos, int i, int ip1,  int j )
{
	//calc gradient
	if(i != j && ip1 != j  )
	{
		double u0 = pos[i*2];
		double v0 = pos[i*2+1];

		double u1 = pos[ip1*2];
		double v1 = pos[ip1*2+1];

		double x = pos[j*2];
		double y = pos[j*2+1];
		double alpha = alpha_val;

		double len01 = sqrt( (u0-u1)*(u0-u1) + (v0 - v1)*(v0 - v1) );
		double len0p = sqrt( (u0-x)*(u0-x) + (v0 - y)*(v0 - y) );
		double len1p = sqrt( (u1-x)*(u1-x) + (v1 - y)*(v1 - y) );
#ifdef ALPHA_TIMES_X
		double same = alpha / pow(( -1*len01 + len0p + len1p ) ,2.0);
		
		newgradient[i*2] -= same * ( (u1-u0)/len01 + (u0-x)/len0p );
		newgradient[i*2+1] -= same * ( (v1-v0)/len01 + (v0-y)/len0p );

		newgradient[ip1*2] -= same * ( (u0-u1)/len01 + (u1-x)/len1p );
		newgradient[ip1*2+1] -= same * ( (v0-v1)/len01 + (v1-y)/len1p );

		newgradient[j*2] -= same * ( (x-u0)/len0p + (x-u1)/len1p );
		newgradient[j*2+1] -= same * ( (y-v0)/len0p + (y-v1)/len1p );

#else
		double same = alpha * pow(1.0f /( -1*len01 + len0p + len1p ) , 1+alpha);

		newgradient[i*2] += same * ((u0 - u1) / len01  + (x - u0) / len0p);
		newgradient[i*2+1] += same * ((v0 - v1) / len01  + (y - v0) / len0p);

		newgradient[ip1*2] += same * ((u1 - u0) / len01  + (x -u1) / len1p);
		newgradient[ip1*2+1] += same * ((v1 - v0) / len01  + (y -v1) / len1p);

		newgradient[j*2] += same * ((u0 - x) / len0p  + (u1 - x) / len1p);
		newgradient[j*2+1] += same * ((v0 - y) / len0p  + (v1 - y) / len1p);
#endif


	}
}

#else

#define SIMPLE_GRAD_EVAL

void BarrierMethod::global_barrier_grad(double *newgradient, const double *pos, int i, int ip1, int j)
{
	int i0 = indchange[0][i];
	int i1 = indchange[0][ip1];
	int i2 = indchange[0][j];

	if (i0 == -1 && i1 == -1 && i2 == -1)
		return;

	//calc gradient
	if (i != j && ip1 != j)
	{
		/*vect2d a;
		a[0] = i0 == -1 ? full_mesh[i*2] : pos[i0*2];
		a[1] = i0 == -1 ? full_mesh[i*2+1] : pos[i*2+1];

		vect2d b;
		b[0] = i1 == -1 ? full_mesh[ip1*2] : pos[i1*2];
		b[1] = i1 == -1 ? full_mesh[ip1*2+1] : pos[i1*2+1];

		vect2d p;
		p[0] = i2 == -1 ? full_mesh[j*2] : pos[i2*2];
		p[1] = i2 == -1 ? full_mesh[j*2+1] : pos[i2*2+1];
		*/
		vect2d a;
		a[0] = pos[i0 * 2];
		a[1] = pos[i * 2 + 1];

		vect2d b;
		b[0] = pos[i1 * 2];
		b[1] = pos[i1 * 2 + 1];

		vect2d p;
		p[0] = pos[i2 * 2];
		p[1] = pos[i2 * 2 + 1];

		/*if(i == 762 && ip1 == 767 && j == 2916 && ::g.printInfo )
		{
		printf("{{%f, %f}, {%f, %f}, {%f, %f}};\n", a[0], a[1], b[0], b[1], p[0], p[1] );
		}*/

		double clamp_dist = CONST_CLAMP;
		double eps_slope = EPSION_FUNC;

		double ll = (a - b).length2();
		double l = sqrt(ll);
		vect2d ortho;
		ortho[0] = (b[1] - a[1]) / l;
		ortho[1] = (a[0] - b[0]) / l;
		double distA = (p - a).length2();
		double distB = (p - b).length2();
		double orthoDist;
		
		if (distA > distB)
		{
			orthoDist = (p - a).dot(ortho);
		}
		else
		{
			orthoDist = (p - b).dot(ortho);
		}
//		orthoDist = (p - b).dot(ortho);
		double dist = abs(orthoDist);
		if (dist >= clamp_dist)
			return;

		double t = ((p - a).dot(b - a)) / ll;
		if (t < 0.0)
		{
//			dist = (p - a).length();
			dist = sqrt(distA);
		}
		else if (t > 1.0)
		{
//			dist = (p - b).length();
			dist = sqrt(distB);
		}
		if (dist < clamp_dist)
		{
			if (t < 0.0)
			{
#ifdef SIMPLE_GRAD_EVAL
				double alpha = 2.0 * clamp_dist * (clamp_dist - dist) * pow(dist, -4.0);

				newgradient[i0 * 2] += (p[0] - a[0]) * alpha;
				newgradient[i0 * 2 + 1] += (p[1] - a[1]) * alpha;

				newgradient[i2 * 2] += (a[0] - p[0]) * alpha;
				newgradient[i2 * 2 + 1] += (a[1] - p[1]) * alpha;
#else
				//	if(i0 != -1)
				{
					newgradient[i0 * 2] += ((2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*((a[0]) + ((-1.0)*(p[0])))*(((-1.0)*(clamp_dist)) + pow(((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1])))))*((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1]))))))));
					newgradient[i0 * 2 + 1] += ((2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*(((-1.0)*(clamp_dist)) + pow(((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1])))))*((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1])))))))*((a[1]) + ((-1.0)*(p[1]))));
				}

				//	if(i1 != -1)
				{
					newgradient[i1 * 2] += 0;
					newgradient[i1 * 2 + 1] += 0;
				}

				//	if(i2 != -1)
				{
					newgradient[i2 * 2] += ((-2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*((a[0]) + ((-1.0)*(p[0])))*(((-1.0)*(clamp_dist)) + pow(((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1])))))*((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1]))))))));
					newgradient[i2 * 2 + 1] += ((-2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*(((-1.0)*(clamp_dist)) + pow(((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1])))))*((((a[0]) + ((-1.0)*(p[0])))*((a[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(p[1])))*((a[1]) + ((-1.0)*(p[1])))))))*((a[1]) + ((-1.0)*(p[1]))));
				}
#endif
			}
			else if (t > 1.0)
			{
#ifdef SIMPLE_GRAD_EVAL
				double alpha = 2.0 * clamp_dist * (clamp_dist - dist) * pow(dist, -4.0);

				newgradient[i1 * 2] += (p[0] - b[0]) * alpha;
				newgradient[i1 * 2 + 1] += (p[1] - b[1]) * alpha;

				newgradient[i2 * 2] += (b[0] - p[0]) * alpha;
				newgradient[i2 * 2 + 1] += (b[1] - p[1]) * alpha;
#else
				//	if(i0 != -1)
				{
					newgradient[i0 * 2] += 0;
					newgradient[i0 * 2 + 1] += 0;
				}

				//	if(i1 != -1)
				{
					newgradient[i1 * 2] += ((2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*((b[0]) + ((-1.0)*(p[0])))*(((-1.0)*(clamp_dist)) + pow(((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1])))))*((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1]))))))));
					newgradient[i1 * 2 + 1] += ((2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*(((-1.0)*(clamp_dist)) + pow(((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1])))))*((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1])))))))*((b[1]) + ((-1.0)*(p[1]))));
				}

				//	if(i2 != -1)
				{
					newgradient[i2 * 2] += ((-2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*((b[0]) + ((-1.0)*(p[0])))*(((-1.0)*(clamp_dist)) + pow(((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1])))))*((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1]))))))));
					newgradient[i2 * 2 + 1] += ((-2.0)*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*(((-1.0)*(clamp_dist)) + pow(((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1]))))), (1.0) / (2.0)))*1.0 / ((((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1])))))*((((b[0]) + ((-1.0)*(p[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((b[1]) + ((-1.0)*(p[1])))*((b[1]) + ((-1.0)*(p[1])))))))*((b[1]) + ((-1.0)*(p[1]))));
				}

#endif
			}
			else
			{
#ifdef SIMPLE_GRAD_EVAL
				double alpha = 2.0 * clamp_dist * (clamp_dist - dist) * pow(dist, -3.0) / l;
				if (orthoDist != dist)
				{
					newgradient[i2 * 2] -= (a[1] - b[1]) * alpha;
					newgradient[i2 * 2 + 1] -= (b[0] - a[0]) * alpha;

					newgradient[i1 * 2] -= (b[1] - a[1]) * (t)* alpha;
					newgradient[i1 * 2 + 1] -= (a[0] - b[0]) * (t)* alpha;

					t = ((p - b).dot(b - a)) / ll;
					newgradient[i0 * 2] -= (a[1] - b[1]) * (t)* alpha;
					newgradient[i0 * 2 + 1] -= (b[0] - a[0]) * (t)* alpha;
				}
				else
				{
					newgradient[i2 * 2] += (a[1] - b[1]) * alpha;
					newgradient[i2 * 2 + 1] += (b[0] - a[0]) * alpha;

					newgradient[i1 * 2] += (b[1] - a[1]) * (t)* alpha;
					newgradient[i1 * 2 + 1] += (a[0] - b[0]) * (t)* alpha;

					t = ((p - b).dot(b - a)) / ll; // t-1
					newgradient[i0 * 2] += (a[1] - b[1]) * (t)* alpha;
					newgradient[i0 * 2 + 1] += (b[0] - a[0]) * (t)* alpha;

				}
#else
				//	if(i2 != -1)
				{
					newgradient[i2 * 2] += ((-2.0)*((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))*((a[1]) + ((-1.0)*(b[1])))*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*1.0 / (((((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))))*(((-1.0)*(clamp_dist)) + pow((1.0 / ((((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))))*((((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1])))*(((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1]))))), (1.0) / (2.0))));
					newgradient[i2 * 2 + 1] += ((2.0)*((a[0]) + ((-1.0)*(b[0])))*((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*1.0 / (((((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))))*(((-1.0)*(clamp_dist)) + pow((1.0 / ((((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))))*((((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1])))*(((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1]))))), (1.0) / (2.0))));
				}

				//	if(i0 != -1)
				{
					newgradient[i0 * 2] += ((2.0)*((a[1]) + ((-1.0)*(b[1])))*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*((((a[0]) + ((-1.0)*(b[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((b[1]) + ((-1.0)*(p[1])))))*1.0 / (((((-1.0)*(b[1])*(p[0])) + ((a[1])*(((-1.0)*(b[0])) + (p[0]))) + ((a[0])*((b[1]) + ((-1.0)*(p[1])))) + ((b[0])*(p[1])))*(((-1.0)*(b[1])*(p[0])) + ((a[1])*(((-1.0)*(b[0])) + (p[0]))) + ((a[0])*((b[1]) + ((-1.0)*(p[1])))) + ((b[0])*(p[1])))*(((-1.0)*(b[1])*(p[0])) + ((a[1])*(((-1.0)*(b[0])) + (p[0]))) + ((a[0])*((b[1]) + ((-1.0)*(p[1])))) + ((b[0])*(p[1])))))*(((-1.0)*(clamp_dist)) + pow((1.0 / ((((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))))*((((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1])))*(((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1]))))), (1.0) / (2.0))));
					newgradient[i0 * 2 + 1] += ((2.0)*((a[0]) + ((-1.0)*(b[0])))*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*((((a[0]) + ((-1.0)*(b[0])))*((b[0]) + ((-1.0)*(p[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((b[1]) + ((-1.0)*(p[1])))))*1.0 / (((((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))))*(((-1.0)*(clamp_dist)) + pow((1.0 / ((((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))))*((((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1])))*(((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1]))))), (1.0) / (2.0))));
				}

				//	if(i1 != -1)
				{
					newgradient[i1 * 2] += ((2.0)*((a[1]) + ((-1.0)*(b[1])))*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*1.0 / (((((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))))*(((a[0])*(a[0])) + ((a[1])*(a[1])) + ((b[0])*(p[0])) + ((-1.0)*(a[0])*((b[0]) + (p[0]))) + ((b[1])*(p[1])) + ((-1.0)*(a[1])*((b[1]) + (p[1]))))*(((-1.0)*(clamp_dist)) + pow((1.0 / ((((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))))*((((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1])))*(((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1]))))), (1.0) / (2.0))));
					newgradient[i1 * 2 + 1] += ((-2.0)*((a[0]) + ((-1.0)*(b[0])))*(clamp_dist)*1.0 / (((eps_slope)*(eps_slope)))*1.0 / (((((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))*(((a[1])*((b[0]) + ((-1.0)*(p[0])))) + ((b[1])*(p[0])) + ((-1.0)*(b[0])*(p[1])) + ((a[0])*(((-1.0)*(b[1])) + (p[1]))))))*(((a[0])*(a[0])) + ((a[1])*(a[1])) + ((b[0])*(p[0])) + ((-1.0)*(a[0])*((b[0]) + (p[0]))) + ((b[1])*(p[1])) + ((-1.0)*(a[1])*((b[1]) + (p[1]))))*(((-1.0)*(clamp_dist)) + pow((1.0 / ((((((a[0]) + ((-1.0)*(b[0])))*((a[0]) + ((-1.0)*(b[0])))) + (((a[1]) + ((-1.0)*(b[1])))*((a[1]) + ((-1.0)*(b[1])))))))*((((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1])))*(((a[1])*(b[0])) + ((-1.0)*(a[0])*(b[1])) + ((-1.0)*(a[1])*(p[0])) + ((b[1])*(p[0])) + ((a[0])*(p[1])) + ((-1.0)*(b[0])*(p[1]))))), (1.0) / (2.0))));
				}
#endif
			}
		}
	}
}
#endif