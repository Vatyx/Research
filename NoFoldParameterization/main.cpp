#include "def_common.h"
#include "global.h"
#include "timer.h"
#include <math.h>
#include <string>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include "JasonFull.h"
#include <omp.h>
#include "vect.h"

#include <fstream>

ParamMethod *param;
static bool first_collapse = true;

void flapOptimization()
{
	vector<unsigned int> vertexCount(g.halfmesh.totalVerts);
	vector<HE_Face*> faceTracker(g.halfmesh.totalVerts);

	map<unsigned int, HE_Face*> flaps;

	for(unsigned int i = 0; i < g.halfmesh.totalFaces; i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			vertexCount[g.halfmesh.faceData[i].vi.v[j]]++;
			faceTracker[g.halfmesh.faceData[i].vi.v[j]] = &g.halfmesh.faceData[i];
		}
	}

	for(unsigned int i = 0; i < vertexCount.size(); i++)
	{
		if(vertexCount[i] == 1)
		{
			flapOpSave save;
			vect3i verts;
			verts[0] = i;

			HE_Vertex p1;
			HE_Vertex p2;
			HE_Vertex p3;
			if(i == faceTracker[i]->vi[0])
			{
				p3 = g.halfmesh.vertexData[i];

				p1 = g.halfmesh.vertexData[faceTracker[i]->vi[1]];
				p2 = g.halfmesh.vertexData[faceTracker[i]->vi[2]];

				verts[1] = faceTracker[i]->vi[1];
				verts[2] = faceTracker[i]->vi[2];
			}
			else if(i == faceTracker[i]->vi[1])
			{
				p3 = g.halfmesh.vertexData[i];

				p1 = g.halfmesh.vertexData[faceTracker[i]->vi[2]];
				p2 = g.halfmesh.vertexData[faceTracker[i]->vi[0]];

				verts[1] = faceTracker[i]->vi[2];
				verts[2] = faceTracker[i]->vi[0];
			}
			else
			{
				p3 = g.halfmesh.vertexData[i];

				p1 = g.halfmesh.vertexData[faceTracker[i]->vi[0]];
				p2 = g.halfmesh.vertexData[faceTracker[i]->vi[1]];

				verts[1] = faceTracker[i]->vi[0];
				verts[2] = faceTracker[i]->vi[1];
			}

			auto normal = (p2.pos - p1.pos).cross(p3.pos - p1.pos);
			normal /= normal.length();
			auto p2p1perpr3 = normal.cross(p2.pos - p1.pos);

			auto alpha = ( (p3.pos - p1.pos).dot(p2.pos - p1.pos) ) / ((p2.pos - p1.pos).length2());
			auto L = (p3.pos - p1.pos).dot(p2p1perpr3) / (p2.pos - p1.pos).length();

			save.alpha = alpha;
			save.L = L;

			if(L < 0)
			{
				printf("L is %f\n", L);
			}

			auto p2p1perp = p2.uvpos - p1.uvpos;
			auto temp = p2p1perp[0];
			p2p1perp[0] = -p2p1perp[1];
			p2p1perp[1] = temp;
			
			p3.uvpos = (1 - alpha)*p1.uvpos + alpha * p2.uvpos + (p2p1perp * L / (p2.uvpos - p1.uvpos).length());

			g.halfmesh.vertexData[i] = p3;

			save.verts = verts;
			g.flapLoopup.push_back(save);
		}

	}
}

void run_full_optimization()
{
	printf("Finish collapse and recompute floaters\n");

		g.halfmesh.gen_collapse_lvl();

		vect3d maxbox = g.halfmesh.vertexData[0].pos;
		vect3d minbox = g.halfmesh.vertexData[0].pos;

		for (int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			if(maxbox.v[0] < g.halfmesh.vertexData[i].pos.v[0])
			{
				maxbox.v[0] = g.halfmesh.vertexData[i].pos.v[0];
			}

			if(maxbox.v[1] < g.halfmesh.vertexData[i].pos.v[1])
			{
				maxbox.v[1] = g.halfmesh.vertexData[i].pos.v[1];
			}

			if(maxbox.v[2] < g.halfmesh.vertexData[i].pos.v[2])
			{
				maxbox.v[2] = g.halfmesh.vertexData[i].pos.v[2];
			}

			if(minbox.v[0] > g.halfmesh.vertexData[i].pos.v[0])
			{
				minbox.v[0] = g.halfmesh.vertexData[i].pos.v[0];
			}

			if(minbox.v[1] > g.halfmesh.vertexData[i].pos.v[1])
			{
				minbox.v[1] = g.halfmesh.vertexData[i].pos.v[1];
			}

			if(minbox.v[2] > g.halfmesh.vertexData[i].pos.v[2])
			{
				minbox.v[2] = g.halfmesh.vertexData[i].pos.v[2];
			}
		}

		float maxdistance = 0;
		vect3d totalDistance = maxbox - minbox;

		for(int i = 0; i < 3; i++)
		{
			if(maxdistance < totalDistance.v[i])
			{
				maxdistance = totalDistance.v[i];
			}
		}
		
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			auto& vertex = g.halfmesh.vertexData[i].pos;
			
			vertex -= (maxbox + minbox) / 2;
			vertex /= maxdistance;
		}

#ifdef PERFORM_FLOATERS
		for(int ii = 0; ii < g.halfmesh.boundaries.size(); ii++)
		{
			g.halfmesh.tuttes_embedding_chart(ii);
		}

#endif
		g.halfmesh.scale_halmesh(100);

		for (int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			g.mesh.positions[i] = g.halfmesh.vertexData[i].pos;
			g.mesh.tex_coords[i] = g.halfmesh.vertexData[i].uvpos;
		}
		printf("num levels debug %d\n", g.halfmesh.col_lev.size());

		printf("finding and optimizing flaps\n");
		flapOptimization();
		printf("first triangle: (%f, %f), (%f, %f), (%f, %f)\n", 
			g.halfmesh.vertexData[g.halfmesh.faceData[0].vi[0]].uvpos[0], g.halfmesh.vertexData[g.halfmesh.faceData[0].vi[0]].uvpos[1], 
			g.halfmesh.vertexData[g.halfmesh.faceData[0].vi[1]].uvpos[0], g.halfmesh.vertexData[g.halfmesh.faceData[0].vi[1]].uvpos[1],
			g.halfmesh.vertexData[g.halfmesh.faceData[0].vi[2]].uvpos[0], g.halfmesh.vertexData[g.halfmesh.faceData[0].vi[2]].uvpos[1]);
		//return;

		g.halfmesh.col_lev.clear();
		g.halfmesh.gen_collapse_lvl();

		g.halfmesh.analyze_boundary();
		
		g.halfmesh.err.clear();
		
		double tot = 0;
		double maxi = -99999;
		for(int i = 0; i < g.halfmesh.totalFaces;i++)
		{
			double t = g.tri_error( g.halfmesh.col_lev[0].iso[i][1][0],
				g.halfmesh.col_lev[0].iso[i][2][0],
				g.halfmesh.col_lev[0].iso[i][2][1],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos[0],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos[1],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos[0],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos[1],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos[0],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos[1]
				) ;
				
				t = t > 10000 && false ? 10000 : t;
			g.halfmesh.err.push_back( .2*-1*log( 1.0f / (t+1)));
			tot += t;
			if(t > maxi)
				maxi = t;
		}

		printf("max err = %f.. avg = %f\n", maxi, tot/g.halfmesh.totalFaces);

		g.halfmesh.bin_errors.clear();
		g.halfmesh.bins.clear();
		for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			g.halfmesh.bin_errors.push_back(0);
			g.halfmesh.bins.push_back(0);
		}
		vect2d s;
		s[0] = 1;s[1] = 0;
		for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			g.halfmesh.bins[i] = g.find_ang(s, g.halfmesh.boundaries[0][i]->uvpos);
		}

		max_iters_t =  10000000;
		printf("run optimization\n");
		//run opt 
		
		g.debug_dirs.clear();
		g.debug_p.clear();

		if(!param->is_full_created)
		{
			param->is_full_created = true;
			param->full_mesh = new double[g.halfmesh.totalVerts*2];
			g.perform_isometric_flattening(&g.mesh, &g.iso_tris);
			param->init_lbfgs( g.halfmesh.totalVerts );
		}

		vector<int> indc;
		vector<vector<int>> boo; boo.push_back(g.boundaryorder);

		vector<vector<vect2d>> isos_to_opT;
		vector<vect3i> faces_to_op;
		vector<vector<vect2d>> isos_to_op;

		vector<int> remove_ind;
		vector<int> tri_inds;

		g.remove_these_tris.clear();
		printf("Remove Tris = %d\n", g.remove_these_tris.size());
		for(int i = 0; i < g.remove_these_tris.size(); i++)
		{
			for(int j = 0; j < 3; j++)
			{
				int ind = g.halfmesh.faceData[g.remove_these_tris[i]].vi[j];
			//	printf("%d\n", ind);
				bool isin = false;
				for(int ii = 0; ii < remove_ind.size(); ii++)
				{
					if(ind == remove_ind[ii])
					{
						isin = true; break;
					}
				}
				if(!isin)
				{
					remove_ind.push_back(ind);
				}
			}
		}

		int ind = 0;

		//Get indchange
		ind = 0;
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			bool isin = false;
			for(int ii = 0; ii < remove_ind.size(); ii++)
			{
				if(i == remove_ind[ii])
				{
					isin = true; break;
				}
			}
			if(isin)
			{
				indc.push_back(-1);
			}
			else
			{
				indc.push_back(ind);
				ind++;
			}
		}

		for(int i = 0; i < g.halfmesh.totalFaces; i++)
		{
			bool isin = false;
			for(int ii = 0; ii < g.remove_these_tris.size(); ii++)
			{
				if(i == g.remove_these_tris[ii])
				{
					isin = true; break;
				}
			}
			if(!isin)
			{
				vect3i t = g.halfmesh.faceData[i].vi;
				faces_to_op.push_back(t);
				
			}
		}
		printf("total faces = %d, faces to opt = %d / %d\n", g.halfmesh.totalFaces, faces_to_op.size(), isos_to_opT.size());
		printf("total verts = %d, verts to opt = %d\n", g.halfmesh.totalVerts, ind);

		//set collapse level
		param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
			ind,
			&g.halfmesh.col_lev[g.expand_iter].faces,
			&indc, 
			&boo,
			&g.halfmesh.col_lev[g.expand_iter].iso);

		param->set_seamless_param(g.seamless);

		//Find matching boundaries
		param->match_boundary.clear();
		//param->match_v.clear();
		vector<int> inserted;
		inserted.clear();
		for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			MatchingBoundary t;

			bool skipme = false;
			for(int j = 0; j < inserted.size(); j++)
			{
				if(i == inserted[j])
				{
					skipme = true;
				}
			}
			if(skipme)
				continue;
			int match1 = i;
			int match2 = (i+1)%g.halfmesh.boundaries[0].size();
			for(int j = 0; j < g.halfmesh.boundaries[0].size(); j++)
			{

					double len = (g.halfmesh.boundaries[0][i]->pos - g.halfmesh.boundaries[0][j]->pos).length();
					double len2 = (g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->pos - g.halfmesh.boundaries[0][(j-1)%g.halfmesh.boundaries[0].size()]->pos).length();


					if((i == 20 && len < .0001) || i == 20 && len2 < .0001 )
						printf("checker");
					if(len < 0.00000000001 && len2 < .000000000001)
					{
						if(i == j && (i+1)%g.halfmesh.boundaries[0].size() == (j-1)%g.halfmesh.boundaries[0].size())
							printf("wtf seams\n");

						if(i == (j-1)%g.halfmesh.boundaries[0].size() && (i+1)%g.halfmesh.boundaries[0].size() == j)
							printf("wtf seams\n");

						match1 = j;
						match2 = (j-1)%g.halfmesh.boundaries[0].size();
						break;
					}
			}


			t.u1 = g.halfmesh.boundaries[0][i]->index;
			t.v1 = g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->index;

			t.u2 = g.halfmesh.boundaries[0][ match1 ]->index;
			t.v2 = g.halfmesh.boundaries[0][ match2 ]->index;

			bool insert = true;
			if(t.u1 == t.u2)
				insert = false;
			
			if(insert)
			{
				param->match_boundary.push_back(t);
				inserted.push_back(i);
				inserted.push_back(match2);
			}
		}

		printf("matching boundaries = %d\n", param->match_boundary.size());


		for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
		{
			param->full_mesh[ii*2] = g.halfmesh.vertexData[ii].uvpos[0];
			param->full_mesh[ii*2+1] = g.halfmesh.vertexData[ii].uvpos[1];
		}

		double* temppos = new double[ g.halfmesh.totalVerts * 2];
		double* sol = new double[ g.halfmesh.totalVerts * 2];

		ind = 0;
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			if(indc[i] != -1)
			{
				temppos[2*ind] = g.halfmesh.vertexData[i].uvpos[0];
				temppos[2*ind+1] = g.halfmesh.vertexData[i].uvpos[1];
				ind++;
			}
		}

		double start = get_time();

		int ret = -1;
		//((JasonFull*)param)->runHessian(10000, temppos, sol);
		((JasonFull*)param)->runKVF(1000, temppos, sol);

		for(int i = 0; i < ind*2; i++)
		{
			temppos[i] = sol[i];
		}
		double end = get_time();

		g.timing = end - start;
		g.total_timing += g.timing;
		printf("Total time = %f, iteration = %f\n", g.total_timing, g.timing);
		param->print_timings();

		ind = 0;
		double movement = 0;
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			if(indc[i] != -1)
			{
				movement += (g.halfmesh.vertexData[i].uvpos[0] - temppos[2*ind])*(g.halfmesh.vertexData[i].uvpos[0] - temppos[2*ind]);
				movement += (g.halfmesh.vertexData[i].uvpos[1] - temppos[2*ind+1])*(g.halfmesh.vertexData[i].uvpos[1] - temppos[2*ind+1]);
				g.mesh.tex_coords[i][0] = g.halfmesh.vertexData[i].uvpos[0] = temppos[2*ind];
				g.mesh.tex_coords[i][1] = g.halfmesh.vertexData[i].uvpos[1] = temppos[2*ind+1];
				ind++;
			}
		}

		printf("movement = %.15g\n", movement);
		if(movement < REMOVE_DIS && g.remove_tri_index != -1)
		{
			printf("remove triangle = %d\n", g.remove_tri_index );
			g.remove_these_tris.push_back( g.remove_tri_index );
		}


		double min = 9999999999999;
		double max = -999999;
		tot = 0;

		for(int i = 0; i < param->match_boundary.size(); i++)
		{
			double err = param->match_boundary[i].p1*param->match_boundary[i].p2*param->match_boundary[i].p3*param->match_boundary[i].p4;

			tot += err;
			min = err < min ? err : min;
			max = err > max ? err : max;
		}
		
		printf("seam error = %.10g, min = %.10g, max = %.10g\n", tot, min,max);

		//printf("clean up\n");
		delete[] temppos;
		delete[] sol;
		//printf("clean up\n");
}

void test_distribution()
{
	g.halfmesh.peaks.clear();
	for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			g.halfmesh.bin_errors[i] = 0;
		}



	double maxi = -999999999;
		for(int i = 0; i < g.halfmesh.totalFaces;i++)
		{
			vect2d mid = (g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos +
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos +
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos ) / 3.0f;

			double t = mid.length();
			if(t > maxi)
				maxi = t;

		}

		//g.takefrom = .1;
		printf("rad = %f, take from = %f\n", maxi, g.takefrom);

		g.halfmesh.error_view.clear();
		g.halfmesh.error_ind.clear();
		vect2d s;
		s[0] = 1;
		s[1] = 0;
		for(int i  = 0; i <	g.halfmesh.totalFaces; i++)
		{
			vect2d mid = (g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos +
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos +
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos ) / 3.0f;

			double t = mid.length();

			if(t >= maxi * g.takefrom)
			{
				vect2d tt;
				tt[0] = g.find_ang( s, mid ) * 10;
				tt[1] = g.halfmesh.err[i]* 10;
				g.halfmesh.error_view.push_back(tt);

				//if(g.halfmesh.err[i] > .25)
				g.halfmesh.error_ind.push_back(i);

				float minb = 999999;
				float maxb = -999999;

				float tang = g.find_ang(s, g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos);
				if(minb > tang)
					minb = tang;
				if(maxb < tang)
					maxb = tang;

				tang = g.find_ang(s, g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos);
				if(minb > tang)
					minb = tang;
				if(maxb < tang)
					maxb = tang;

				tang = g.find_ang(s, g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos);
				if(minb > tang)
					minb = tang;
				if(maxb < tang)
					maxb = tang;

				if(maxb - minb > PI)//wrap around
				{
					float te = minb;
					minb = maxb;
					maxb += 2*PI + minb;
				}

				//printf("(%f, %f)", minb, maxb);

				bool found = false;
				for(int j = 0; j < g.halfmesh.bins.size(); j++)
				{
					/*if(tt[0]/10.0f > g.halfmesh.bins[j] && tt[0]/10.0f < g.halfmesh.bins[j+1])
					{
						g.halfmesh.bin_errors[j] += tt[1];
						found = true;
						break;
					}*/

					/*if(j == g.halfmesh.bins.size()-1)
					{
						if((minb <= g.halfmesh.bins[j] && maxb >= g.halfmesh.bins[j]) || (minb <= 2*PI + g.halfmesh.bins[0] && maxb >= 2*PI + g.halfmesh.bins[0]) )
						{
							found = true;
							g.halfmesh.bin_errors[j] += tt[1];
						}

						if( maxb > g.halfmesh.bins[j] )
							break;
					}
					else
					{
						if((minb <= g.halfmesh.bins[j] && maxb >= g.halfmesh.bins[j]) || (minb <= g.halfmesh.bins[j+1] && maxb >= g.halfmesh.bins[j+1]) )
						{
							found = true;
							g.halfmesh.bin_errors[j] += tt[1];
						}

						if( maxb > g.halfmesh.bins[j] )
							break;
					}*/

					if(j == g.halfmesh.bins.size()-1)
					{
						if((minb >= g.halfmesh.bins[j] && minb <= 2*PI + g.halfmesh.bins[0]) || (maxb >= g.halfmesh.bins[j] && maxb <= 2*PI + g.halfmesh.bins[0]) )
						{
							found = true;
							g.halfmesh.bin_errors[j] += tt[1];
						}

						if((minb <= g.halfmesh.bins[j] && maxb >= 2*PI + g.halfmesh.bins[0]) )
						{
							found = true;
							g.halfmesh.bin_errors[j] += tt[1];
						}

						if( maxb < g.halfmesh.bins[j] )
							break;
					}
					else
					{
						if((minb >= g.halfmesh.bins[j] && minb <= g.halfmesh.bins[j+1]) || (maxb >= g.halfmesh.bins[j] && maxb <= g.halfmesh.bins[j+1]) )
						{
							found = true;
							g.halfmesh.bin_errors[j] += tt[1];
						}

						if((minb <= g.halfmesh.bins[j] && maxb >= g.halfmesh.bins[j+1]) )
						{
							found = true;
							g.halfmesh.bin_errors[j] += tt[1];
						}

						if( maxb < g.halfmesh.bins[j] )
							break;
					}
				}
				/*if(!found)
					g.halfmesh.bin_errors[g.halfmesh.bins.size()-1] += tt[1];*/
			}

		}

		double binsum = 0;
		for(int i = 0; i < g.halfmesh.bins.size(); i++)
		{
			binsum += g.halfmesh.bin_errors[i]*g.halfmesh.bin_errors[i];
		}
		printf("error sum = %f\n", sqrt(binsum));

		for(int i = 0; i < g.halfmesh.bins.size(); i++)
		{
			g.halfmesh.bin_errors[i] /= sqrt(binsum);
		}

		binsum = 0;
		for(int i = 0; i < g.halfmesh.bins.size(); i++)
		{
			binsum += g.halfmesh.bin_errors[i]*g.halfmesh.bin_errors[i];
		}
		double mean = binsum / g.halfmesh.bins.size();
		double sd = 0;

		for(int i = 0; i < g.halfmesh.bins.size(); i++)
		{
			sd += (g.halfmesh.bin_errors[i] - mean)*(g.halfmesh.bin_errors[i] - mean);
		}
		sd = sqrt(sd / g.halfmesh.bin_errors.size());

		printf("mean error = %f, sd = %f\n", mean, sd);
		float prev = g.halfmesh.bin_errors[0];
		int poten_peak = 0;
		bool possiblepeak = false;
		for(int i = 0; i < g.halfmesh.bins.size(); i++)
		{
		//	printf("%f ",g.halfmesh.bin_errors[i]);

			if(g.halfmesh.bin_errors[i] > prev )
			{
				prev = g.halfmesh.bin_errors[i];
				poten_peak = i;
				possiblepeak = true;
				continue;
			}

			if( possiblepeak && ((prev - g.halfmesh.bin_errors[i]) > (sd / 2.0f)) )
			{
				g.halfmesh.peaks.push_back(poten_peak);
				possiblepeak = false;
			}

			if(!possiblepeak)
			{
				prev = g.halfmesh.bin_errors[i];
				poten_peak = i;
			}

			/*if(g.halfmesh.bin_errors[i] < prev)
			{
				
			}
			prev = g.halfmesh.bin_errors[i];*/
		}
		printf("num peaks = %d\n", g.halfmesh.peaks.size());
		
		
		vector<PeakInterval> intervals;
		for(int i = 0; i < g.halfmesh.peaks.size(); i++)
		{
			printf("\r %d", i);
			PeakInterval tint;
			tint.start = tint.end = tint.peak = g.halfmesh.peaks[i];
			tint.vs = tint.ve = tint.vp = g.halfmesh.bin_errors[g.halfmesh.peaks[i]];
			tint.circ= false;

			intervals.push_back(tint);

			//grow intervals

			bool notstopped = true;
			double stopme = sd / 5.0f;
			
			//grow start
			int j = (tint.peak - 1 + g.halfmesh.bins.size()) % g.halfmesh.bins.size();
			while(notstopped)
			{
				if(intervals[i].vs - g.halfmesh.bin_errors[j] < stopme )
				{
					j = (j - 1 + g.halfmesh.bins.size()) % g.halfmesh.bins.size();
				}
				else if( intervals[i].vs - g.halfmesh.bin_errors[j] >= stopme )
				{
					notstopped = false;
					intervals[i].vs = g.halfmesh.bin_errors[j];
					intervals[i].start = j;
				}
			}

			notstopped = true;
			
			//grow start
			j = (tint.peak + 1 + g.halfmesh.bins.size()) % g.halfmesh.bins.size();
			while(notstopped)
			{
				if(intervals[i].vs - g.halfmesh.bin_errors[j] < stopme )
				{
					j = (j + 1 + g.halfmesh.bins.size()) % g.halfmesh.bins.size();
				}
				else if( intervals[i].ve - g.halfmesh.bin_errors[j] >= stopme )
				{
					notstopped = false;
					intervals[i].ve = g.halfmesh.bin_errors[j];
					intervals[i].end = j;
				}
			}
		}
		

		for(int i = 0; i  < g.halfmesh.peaks.size(); i++)
		{
			g.halfmesh.pi.push_back( intervals[i] );
		}

		vect2d t1; t1[0] = maxi; t1[1] = 0;
		vect2d t2; t2[0] = 0; t2[1] = maxi;
		vect2d t3; t3[0] = -maxi; t3[1] = 0;
		vect2d t4; t4[0] = 0; t4[1] = -maxi;
		g.halfmesh.error_view.push_back(t1);
		g.halfmesh.error_view.push_back(t2);
		g.halfmesh.error_view.push_back(t3);
		g.halfmesh.error_view.push_back(t4);
		printf("num tris shown = %d / %d \n", g.halfmesh.error_view.size(), g.halfmesh.totalFaces);


}

void output_dist_math()
{
	printf("out dist\n");
	string str = "dist.dis";
	ofstream f(str.c_str());

	f << '{';
	for(int i = 0; i < g.halfmesh.error_view.size()-4; i++)
	{
		f << '{' << g.halfmesh.error_view[i][0] << ','
				 << g.halfmesh.error_view[i][1] << '}';
				 
		if(i < g.halfmesh.error_view.size()-5) 
			f << ',';
	}
	f << '}';
	f << endl;

	f.close();
}

void reshape(int width, int height)
{
	glViewport(0, 0, SCREEN_WDITH, SCREEN_HEIGHT);
}
void display()
{
	glClearColor(g.bg_color[0], g.bg_color[1], g.bg_color[2], 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set up matrices
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, double(SCREEN_WDITH)/double(SCREEN_HEIGHT), .01, 1000);
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	//Camera
	glTranslatef(0, 0, -g.zoom);
	glMultMatrixd(g.rotation);
	glTranslatef(-g.focus[0], -g.focus[1], -g.focus[2]);

	if(g.draw_uv01)
	{ //Draw 0,1 bounding box
		glBegin(GL_LINES);

		glColor3f(1,1,1);
		glVertex2f(0,0);
		glVertex2f(1,0);

		glVertex2f(1,0);
		glVertex2f(1,1);

		glVertex2f(1,1);
		glVertex2f(0,1);

		glVertex2f(0,1);
		glVertex2f(0,0);

		glEnd();
	}

	if(g.display_error)
	{
		glBegin(GL_TRIANGLES);
		//Display error triangles
		for(int i = 0; i < g.mesh.indices_tex.s; i++)
		{
			glColor3f(g.err[i], 0,0);
			glVertex2dv(g.mesh.tex_coords[g.mesh.indices_tex[i][0]].v);
			glVertex2dv(g.mesh.tex_coords[g.mesh.indices_tex[i][1]].v);
			glVertex2dv(g.mesh.tex_coords[g.mesh.indices_tex[i][2]].v);
		}
		glEnd();

	}

	
	else if(g.draw_simplified)
	{
		

	//Draw Tex Coodinates
		glColor3f(1,1,0);
		glBegin(GL_LINES);

		//Draw Non Inverted Triangles
		for(int i = 0; i < g.smesh.indices_tex.s; i++)
		{
			glVertex2dv(g.smesh.tex_coords[g.smesh.indices_tex[i][0]].v);
			glVertex2dv(g.smesh.tex_coords[g.smesh.indices_tex[i][1]].v);

			glVertex2dv(g.smesh.tex_coords[g.smesh.indices_tex[i][1]].v);
			glVertex2dv(g.smesh.tex_coords[g.smesh.indices_tex[i][2]].v);

			glVertex2dv(g.smesh.tex_coords[g.smesh.indices_tex[i][2]].v);
			glVertex2dv(g.smesh.tex_coords[g.smesh.indices_tex[i][0]].v);
		}
		glEnd();

		//Draw Inverted Triangles
		glColor3f(1,0,1);
		glBegin(GL_LINES);
		for (int i = 0; i < g.smesh.indices_tex.s; i++)
		{
			vect2d t[3];
			t[0] = g.smesh.tex_coords[g.smesh.indices_tex[i][0]];
			t[1] = g.smesh.tex_coords[g.smesh.indices_tex[i][1]];
			t[2] = g.smesh.tex_coords[g.smesh.indices_tex[i][2]];

			vect2d e[2];
			e[0] = t[1] - t[0];
			e[1] = t[2] - t[0];

			float area = e[0][0]*e[1][1] - e[0][1]*e[1][0];

			float u0 = g.smesh.tex_coords[g.smesh.indices_tex[i][0]][0];
			float v0 = g.smesh.tex_coords[g.smesh.indices_tex[i][0]][1];

			float u1 = g.smesh.tex_coords[g.smesh.indices_tex[i][1]][0];
			float v1 = g.smesh.tex_coords[g.smesh.indices_tex[i][1]][1];

			float u2 = g.smesh.tex_coords[g.smesh.indices_tex[i][2]][0];
			float v2 = g.smesh.tex_coords[g.smesh.indices_tex[i][2]][1];

			float A = (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2));
			//if (area >= 0)
				if (A < 0)
			{
				printf("%d = {{%f,%f},{%f, %f},{%f, %f}}\n", i, u0,v0,u1,v1,u2,v2);
				glVertex2dv(t[0].v);
				glVertex2dv(t[1].v);

				glVertex2dv(t[1].v);
				glVertex2dv(t[2].v);
				
				glVertex2dv(t[2].v);
				glVertex2dv(t[0].v);
			}
		}
		glEnd();

		glColor3f(1,0,1);
		glBegin(GL_POINTS);
		for (int i = 0; i < g.smesh.indices_tex.s; i++)
		{
			vect2d t[3];
			t[0] = g.smesh.tex_coords[g.smesh.indices_tex[i][0]];
			t[1] = g.smesh.tex_coords[g.smesh.indices_tex[i][1]];
			t[2] = g.smesh.tex_coords[g.smesh.indices_tex[i][2]];

			vect2d e[2];
			e[0] = t[1] - t[0];
			e[1] = t[2] - t[0];

			float area = e[0][0]*e[1][1] - e[0][1]*e[1][0];

			float u0 = g.smesh.tex_coords[g.smesh.indices_tex[i][0]][0];
			float v0 = g.smesh.tex_coords[g.smesh.indices_tex[i][0]][1];

			float u1 = g.smesh.tex_coords[g.smesh.indices_tex[i][1]][0];
			float v1 = g.smesh.tex_coords[g.smesh.indices_tex[i][1]][1];

			float u2 = g.smesh.tex_coords[g.smesh.indices_tex[i][2]][0];
			float v2 = g.smesh.tex_coords[g.smesh.indices_tex[i][2]][1];

			float A = (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2));
			//if (area >= 0)
				if (A < 0)
			{
				//printf("{{%f,%f},{%f, %f},{%f, %f}}\n", u0,v0,u1,v1,u2,v2);
				glVertex2dv(t[0].v);
				glVertex2dv(t[1].v);

				glVertex2dv(t[1].v);
				glVertex2dv(t[2].v);
				
				glVertex2dv(t[2].v);
				glVertex2dv(t[0].v);
			}
		}
		glEnd();

		//Draw Boundary
		glColor3f(1,0,0);
		glBegin(GL_LINES);

		//Draw Non Inverted Triangles
		double tot_length = 0;
		double max_length = -999999;
		double min_length = 99999;
		for(int i = 0; i < g.sbo.size(); i++)
		{
			glVertex2dv(g.smesh.tex_coords[g.sbo[i]].v);
			glVertex2dv(g.smesh.tex_coords[g.sbo[(i+1)%g.sbo.size()]].v);

			double t = (g.smesh.tex_coords[g.sbo[i]] - g.smesh.tex_coords[g.sbo[(i+1)%g.sbo.size()]]).length();
			
			tot_length += t;

			max_length = t > max_length ? t : max_length;
			min_length = t < min_length ? t : min_length;
		}
		glEnd();

	}

	if (g.draw_halfmesh)
	{
		if(g.draw_full_mesh)
		{
			g.halfmesh.draw_3Dfull();
		}
		else
		{
		g.halfmesh.draw();

		if(g.draw_grad)
		{
		//	printf("max param = %.10g\n", g.max_param);
			glColor3f(1,1,0);
			glBegin(GL_LINES);

			
			for(int i = 0; i < g.halfmesh.totalVerts; i++)
			{
				if(g.halfmesh.collapsed[i] == -1)
				{
					glVertex2dv(g.halfmesh.vertexData[i].uvpos.v );

					vect2d t;
					//t[0] = g.max_param * g.gradDraw[2*i];
					//t[1] = g.max_param * g.gradDraw[2*i+1];
					t[0] =  g.gradDraw[2*i];
					t[1] =  g.gradDraw[2*i+1];

					glVertex2dv((g.halfmesh.vertexData[i].uvpos + t).v );
				}
			}
			glEnd();
		}
		}
	

		/*glBegin(GL_LINES);
		glColor3f(0, 1, 0);
		for (int i = 0; i < g.debug_dirs.size(); i++)
		{
			glVertex2dv(g.debug_dirs[i][0].v);
			glVertex2dv(g.debug_dirs[i][1].v);
		}
		glEnd();*/

		/*glBegin(GL_POINTS);
		glColor3f(1, 1, 0);
		for (int i = 0; i < g.debug_p.size(); i++)
		{
			if(i == 1)
				glColor3f(0,0,1);
			else
				glColor3f(1, 1, 0);
			glVertex2dv(g.debug_p[i].v);
		}
		glEnd();

		glBegin(GL_POLYGON);
		glColor3f(1, 0, 0);
		glVertex2dv(g.halfmesh.vertexData[2503].uvpos.v);
		glVertex2dv(g.halfmesh.vertexData[2499].uvpos.v);
		glVertex2dv(g.halfmesh.vertexData[3066].uvpos.v);
		glEnd();*/
	}

	//printf("size = %d\n", g.debug_p.size());
	glBegin(GL_POINTS);
	glColor3f(1, 1, 0);
	for (int i = 0; i < g.debug_p.size(); i++)
	{
		glVertex2dv(g.debug_p[i].v);
	}
	glEnd();

	glBegin(GL_POINTS);
	glColor3f(1, 1, 0);
	/*if(g.remove_these_tris.size() > 0)
		printf("%d, %d, %d\n", g.halfmesh.faceData[g.remove_these_tris[0]].vi[0],g.halfmesh.faceData[g.remove_these_tris[0]].vi[1],g.halfmesh.faceData[g.remove_these_tris[0]].vi[2]);*/
	for (int i = 0; i < g.remove_these_tris.size(); i++)
	{
	
		glVertex2dv(g.halfmesh.vertexData[ g.halfmesh.faceData[g.remove_these_tris[i]].vi[0]].uvpos.v);
		glVertex2dv(g.halfmesh.vertexData[ g.halfmesh.faceData[g.remove_these_tris[i]].vi[1]].uvpos.v);
		glVertex2dv(g.halfmesh.vertexData[ g.halfmesh.faceData[g.remove_these_tris[i]].vi[2]].uvpos.v);
	}
	glEnd();

	glFlush ();
	glutSwapBuffers();
}

void calc_display_error()
{
	printf("Not Working\n");
	//ParamMethod* a;
	//a = new JasonFull("jason method\n", g.mesh.tex_coords.s, g.boundaryorder, g.iso_tris);
	//param->ind_tex = g.mesh.indices_tex;

	//double* temppos = new double[g.mesh.tex_coords.s * 2];
	//for(int i = 0; i < g.mesh.tex_coords.s; i++)
	//{
	//	temppos[i*2] = g.mesh.tex_coords[i][0];
	//	temppos[i*2+1] = g.mesh.tex_coords[i][1];
	//}

	//g.err = new double[g.mesh.indices_tex.s];

	//param->display_error(g.err ,temppos);

	////normalize colors
	//double m = g.err[0];
	////printf("m = %f\n", m);
	//for(int i = 1; i < g.mesh.indices_tex.s; i++)
	//{
	//	if( g.err[i] > m )
	//		m = g.err[i];
	//}
	//printf("m = %f\n", m);
	//for(int i = 0; i < g.mesh.indices_tex.s; i++)
	//{
	//	g.err[i] /= m;
	//}

	//delete[] temppos;
}

double find_intersection(vect2d p, vect2d a, vect2d b, vect2d pg)
{
	if ((pg[1] * (a[0] - b[0]) + pg[0] * (b[1] - a[1])) == 0)
		printf(":den = 0.... during expansion");
	return (-1 * (-(1e-10) - (b[0]* (p[1] - a[1]) + p[0]* (a[1] - b[1]) + a[0]* (-1 * p[1] + b[1])))) / (pg[1] *(a[0] - b[0]) + pg[0] *(b[1] - a[1]));
	//return (-1 * (-(1e-10) - (u2* (v0 - v1) + u0* (v1 - v2) + u1* (-1 * v0 + v2)))) / gy0 *(u1 - u2) + gx0 *(v2 - v1);
}
void input( unsigned char key, int x, int y )
{
	if(key == 'r')
	{
		/*if(g.draw_simplified)
			g.flip_uv(&g.smesh);
		else
			g.flip_uv(&g.mesh);*/
		g.halfmesh.flip_tex();
		printf("Reverse UVs\n");
	}
	else if(key == 'e')
	{
		g.debug_seamless = !g.debug_seamless;

		if(g.debug_seamless)
		{
		printf("debug edge = %d    (%d, %d)  ->  (%d, %d) \n", g.which_seamless,
			param->match_boundary[g.which_seamless].u1,param->match_boundary[g.which_seamless].v1,
			param->match_boundary[g.which_seamless].u2,param->match_boundary[g.which_seamless].v2);

		double u1,v1,u2,v2,u3,v3,u4,v4;
		u1 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u1].uvpos[0];
		v1 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u1].uvpos[1];
		u2 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v1].uvpos[0];
		v2 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v1].uvpos[1];

		u3 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u2].uvpos[0];
		v3 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u2].uvpos[1];
		u4 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v2].uvpos[0];
		v4 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v2].uvpos[1];

		double p1 = ((((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
		double p2 = ((((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
		double p3 = ((((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
		double p4 = ((((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));


		printf("error = %.10g\n", p1*p2*p3*p4);

		double ang;
		vect2d vv1;
		vv1[0] = u1 - u2; vv1[1] = v1 - v2;
		vect2d vv2;
		vv2[0] = u3 - u4; vv2[1] = v3 - v4;

		printf("angle = %.10g,  lengths = %.10g, %.10g\n", 180*acos(vv1.dot(vv2) / ( vv1.length() * vv2.length() )) / MATH_PI, vv1.length(), vv2.length() );
		}

	}
	else if(key == 'w')
	{
		g.which_seamless = (g.which_seamless + 1) % param->match_boundary.size() ;
		printf("debug edge = %d    (%d, %d)  ->  (%d, %d) \n", g.which_seamless,
			param->match_boundary[g.which_seamless].u1,param->match_boundary[g.which_seamless].v1,
			param->match_boundary[g.which_seamless].u2,param->match_boundary[g.which_seamless].v2);

		double u1,v1,u2,v2,u3,v3,u4,v4;
		u1 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u1].uvpos[0];
		v1 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u1].uvpos[1];
		u2 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v1].uvpos[0];
		v2 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v1].uvpos[1];

		u3 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u2].uvpos[0];
		v3 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u2].uvpos[1];
		u4 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v2].uvpos[0];
		v4 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v2].uvpos[1];

		double p1 = ((((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
		double p2 = ((((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
		double p3 = ((((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
		double p4 = ((((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));


		printf("error = %.10g\n", p1*p2*p3*p4);

		double ang;
		vect2d vv1;
		vv1[0] = u1 - u2; vv1[1] = v1 - v2;
		vect2d vv2;
		vv2[0] = u3 - u4; vv2[1] = v3 - v4;

		printf("angle = %.10g,  lengths = %.10g, %.10g\n", 180*acos(vv1.dot(vv2) / ( vv1.length() * vv2.length() )) / MATH_PI, vv1.length(), vv2.length() );



	}
	else if(key == 'W')
	{
		g.which_seamless = (g.which_seamless - 1) % param->match_boundary.size() ;
		printf("debug edge = %d    (%d, %d)  ->  (%d, %d) \n", g.which_seamless,
			param->match_boundary[g.which_seamless].u1,param->match_boundary[g.which_seamless].v1,
			param->match_boundary[g.which_seamless].u2,param->match_boundary[g.which_seamless].v2);

		double u1,v1,u2,v2,u3,v3,u4,v4;
		u1 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u1].uvpos[0];
		v1 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u1].uvpos[1];
		u2 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v1].uvpos[0];
		v2 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v1].uvpos[1];

		u3 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u2].uvpos[0];
		v3 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].u2].uvpos[1];
		u4 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v2].uvpos[0];
		v4 = g.halfmesh.vertexData[param->match_boundary[g.which_seamless].v2].uvpos[1];

		double p1 = ((((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2))*((u3)+((-1.0)*(u4))+((-1.0)*(v1))+(v2)))+(((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))*((u1)+((-1.0)*(u2))+(v3)+((-1.0)*(v4)))));
		double p2 = ((((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4)))*((u1)+((-1.0)*(u2))+(u3)+((-1.0)*(u4))))+(((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))*((v1)+((-1.0)*(v2))+(v3)+((-1.0)*(v4)))));
		double p3 = ((((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2)))*((u3)+((-1.0)*(u4))+(v1)+((-1.0)*(v2))))+(((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))*((u1)+((-1.0)*(u2))+((-1.0)*(v3))+(v4))));
		double p4 = ((((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4))*((u1)+((-1.0)*(u2))+((-1.0)*(u3))+(u4)))+(((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))*((v1)+((-1.0)*(v2))+((-1.0)*(v3))+(v4))));


		printf("error = %.10g\n", p1*p2*p3*p4);

		double ang;
		vect2d vv1;
		vv1[0] = u1 - u2; vv1[1] = v1 - v2;
		vect2d vv2;
		vv2[0] = u3 - u4; vv2[1] = v3 - v4;

		printf("angle = %.10g,  lengths = %.10g, %.10g\n", 180*acos(vv1.dot(vv2) / ( vv1.length() * vv2.length() )) / MATH_PI, vv1.length(), vv2.length() );


	}
	else if(key == 's')
	{
		printf("Saving mesh and uvs\n");

		for(int i =0; i < g.mesh.tex_coords.s; i++)
		{
			g.mesh.tex_coords[i] = g.halfmesh.vertexData[i].uvpos;
		}
		string st = SAVETO;

		string mn = MESHNAME;
//		g.save_mesh2(st + mn + "/" + mn +"_" + ".obj");
		g.save_mesh( mn + ".obj");
		//g.save_mesh(SAVETO + mn + "/" + mn + "_"+ ".obj");
		//g.halfmesh.save_mesh();
	}
	else if(key == 'S')
	{
		g.halfmesh.save_mesh2();
	}
	else if (key == 'a')
	{
		param->reduce_alpha();
	}
	else if ( key == 'A')
	{
		param->sharpen_alpha();
	}
	else if ( key == 'd' )
	{
		g.printInfo = !g.printInfo;
	}
	else if(key == '-')
	{
		output_dist_math();
	}
	else if (key == '0')
	{
		g.takefrom+=.0001;
		test_distribution();
	}
	else if (key == ')')
	{
		g.takefrom-=.1;
		test_distribution();
	}
	else if (key == 'z' && false)
	{
		printf("run ring optimization\n");
		//run opt 
		
		g.debug_dirs.clear();
		g.debug_p.clear();

		vector<int> indc;

		for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
		{
			indc.push_back(-1);
		}

		vector<vector<int>> dummyboundary;dummyboundary.clear();

		if(!param->is_full_created)
		{
			param->is_full_created = true;
			param->full_mesh = new double[g.halfmesh.totalVerts*2];
		}

		for(int i = 0; i < g.halfmesh.rings.size()-1; i++)
		{
			for(int j = 0; j < i+1; j++)
			{
				//setup indchange
				for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
				{
					indc[i] = -1;
					param->full_mesh[ii*2] = g.halfmesh.vertexData[ii].uvpos[0];
					param->full_mesh[ii*2+1] = g.halfmesh.vertexData[ii].uvpos[1];
				}

				int ind = 0;
				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					indc[g.halfmesh.rings[j][ii]] = ind;
					ind++;
				}

				param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
					g.halfmesh.rings[j].size(),
					&g.halfmesh.rings_faces[j],
					&indc, 
					j==0 ? &g.halfmesh.col_lev[0].boundaries : &dummyboundary,
					&g.halfmesh.rings_isos[j]);


				double* temppos = new double[g.halfmesh.rings[j].size() * 2];
				double* sol = new double[g.halfmesh.rings[j].size() * 2];
				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					temppos[ii*2] = g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[0];
					temppos[ii*2+1] = g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[1];
					sol[ii*2] = 0;
					sol[ii*2+1] = 0;

				}

				int ret;
				ret = param->run(1000, temppos, sol );

				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[0] = sol[ii*2];
					g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[1] = sol[ii*2+1];
				}

				delete[] temppos;
				delete[] sol;
				 
				/*i = 10000;
				j = 10000;*/

			}
		}
	}
	/*else if(key == 'w')
	{
		printf("reset removed tris\n");
		g.remove_these_tris.clear();
	}*/
	
	else if(key == 'Q')
	{
		double st = get_time();
		max_iters_t =  50;
		printf("run ring optimization\n");
		//run opt 
		
		g.debug_dirs.clear();
		g.debug_p.clear();

		//g.debug_p.push_back(g.halfmesh.vertexData[482].uvpos);

		vector<int> indcF;
		vector<vect3i> faces_to_opF;
		vector<vector<vect2d>> isos_to_opF;
		vector<vector<int>> boundaryF;
		vector<int> tb;
		for(int i = 0; i < g.boundaryorder.size(); i++)
		{
			tb.push_back(g.boundaryorder[i]);
		}
		boundaryF.push_back(tb);
		printf("boundary size = %d\n", g.boundaryorder.size());
		for(int i= 0 ; i < g.halfmesh.totalVerts; i++)
		{
			indcF.push_back(i);
		}
		for(int i = 0; i < g.halfmesh.totalFaces; i++)
		{
			vect3i t;
			t[0] = g.mesh.indices_tex[i][0];
			t[1] = g.mesh.indices_tex[i][1];
			t[2] = g.mesh.indices_tex[i][2];
			faces_to_opF.push_back(t);
		}
		vector<vector<int>> indc;
		vector<vector<int>> dummyboundary;dummyboundary.clear();
		vector<vector<int>> boo;
		//g.boundaryorder.clear();
		//g.get_boundary_order(g.mesh,g.mesh
		boo.push_back(g.boundaryorder);
		vector<vector<vect3i>> faces_to_op;
		vector<vector<vector<vect2d>>> isos_to_op;

		if(!param->is_full_created)
		{
			param->is_full_created = true;
			param->full_mesh = new double[g.halfmesh.totalVerts*2];
			g.perform_isometric_flattening(&g.mesh, &g.iso_tris);
			param->init_lbfgs( g.halfmesh.totalVerts );
		}

		
		for(int j = 0; j < g.halfmesh.rings.size(); j++)
		{
			printf("\rsetup ring %d", j);
			//for(int j = 0; j < i+1; j++)
			{
				vector<int> tempi;
				for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
				{
					tempi.push_back(-1);
				}
				indc.push_back(tempi);

				int ind = 0;
				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					indc[j][g.halfmesh.rings[j][ii]] = ind;
					ind++;
				}
				
				vector<vect3i> ft;
				vector<vector<vect2d>> it;

				for(int ii = 0; ii < g.halfmesh.rings_faces[j].size(); ii++)
				{
					ft.push_back( g.halfmesh.rings_faces[j][ii] );
					it.push_back( g.iso_tris[g.halfmesh.rings_find[j][ii]] );
				}

				if(j > 0)
				{
					for(int ii = 0; ii < g.halfmesh.rings_faces[j-1].size(); ii++)
					{
						ft.push_back( g.halfmesh.rings_faces[j-1][ii] );
						it.push_back( g.iso_tris[g.halfmesh.rings_find[j-1][ii]] );
					}
				}
				
				faces_to_op.push_back(ft);
				isos_to_op.push_back(it);
			
			}
		}


		//printf("18debug verts = %d, %d, %d\n",faces_to_op[18][7782][0],faces_to_op[18][7782][1],faces_to_op[18][7782][2]);
		//printf("19debug verts = %d, %d, %d\n",faces_to_op[19][7782][0],faces_to_op[19][7782][1],faces_to_op[19][7782][2]);
		////////////////////////////////////////////////////////DEBUG
		//497 482 496
		/*printf("%.10g, %.10g, %.10g\n", g.iso_tris[593][1][0],
						g.iso_tris[593][2][0],
						g.iso_tris[593][2][1]);*/
		
		for(int i = 0; i < g.halfmesh.totalFaces; i++)
		{
			g.halfmesh.faceData[i].inside_rings = false;
		}
		/*printf("vert 497\n");
		g.print_debug_info(497);

		printf("vert 492\n");
		g.print_debug_info(492);

		printf("vert 496\n");
		g.print_debug_info(496);*/
		g.halfmesh.iso_tris = g.iso_tris;


		int ttt;

		//scanf("%d", ttt);
		//int i = g.ringI;
		//int j = g.ringJ;
		double* temppos = new double[g.halfmesh.totalVerts * 2];
		double* storegrad = new double[g.halfmesh.totalVerts * 2];
		double* sol = new double[g.halfmesh.totalVerts * 2];

		//for(int i = 0; i < g.halfmesh.rings.size()-1; i++)
		{
			int i = g.halfmesh.rings.size()-1;
			for(int j = 0; j < i+1; j++)
			{
				printf("\r  %d, Optimize Ring %d ",i, j);
				//setup indchange
				for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
				{
					param->full_mesh[ii*2] = g.halfmesh.vertexData[ii].uvpos[0];
					param->full_mesh[ii*2+1] = g.halfmesh.vertexData[ii].uvpos[1];
				}
				
				param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
					g.halfmesh.rings[j].size(),
					&faces_to_op[j],
					&indc[j], 
					//j==0 ? &g.halfmesh.col_lev[0].boundaries : &dummyboundary,
					j == 0 ? &boo : &dummyboundary,
					&isos_to_op[j]);

				//printf("look at me = %d\n", g.halfmesh.rings[j][0]);
				//for(int ii = 0; ii < g.mesh.indices_tex.size(); ii++)
				//{
				//	if(g.mesh.indices_tex[ii][0] == g.halfmesh.rings[j][0] || g.mesh.indices_tex[ii][1] == g.halfmesh.rings[j][0]|| g.mesh.indices_tex[ii][2] == g.halfmesh.rings[j][0])
				//		printf("tri = (%d, %d, %d)\n", g.mesh.indices_tex[ii][0],g.mesh.indices_tex[ii][1],g.mesh.indices_tex[ii][2]);
				//}


				/*for(int ii = 0; ii < faces_to_op[j].size(); ii++)
				{
					if(faces_to_op[j][ii][0] == g.halfmesh.rings[j][0] || faces_to_op[j][ii][1] == g.halfmesh.rings[j][0]|| faces_to_op[j][ii][2] == g.halfmesh.rings[j][0])
						printf("tri = (%d, %d, %d)\n", faces_to_op[j][ii][0],faces_to_op[j][ii][1],faces_to_op[j][ii][2]);
				}*/

				//param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
				//	g.halfmesh.rings[j].size(),
				//	&faces_to_opF,
				//	&indc[j], 
				//	//j==0 ? &g.halfmesh.col_lev[0].boundaries : &dummyboundary,
				//	//j == 0 ? &boo : &dummyboundary,
				//	j == 0 ? &boo : &dummyboundary,
				//	&g.iso_tris);
				
				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					temppos[ii*2] = g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[0];
					temppos[ii*2+1] = g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[1];
					sol[ii*2] = 0;
					sol[ii*2+1] = 0;

				}

				int ret;
				ret = param->run(1000, temppos, sol );

				if( _finite(sol[0]) )
				{
					for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
					{
						g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[0] = sol[ii*2];
						g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[1] = sol[ii*2+1];
					}
				}

				display();
				//break;
			}
			//break;
		}
		delete[] temppos;
		delete[] sol;

		printf("Rings took = %.10g\n", get_time() - st);

	}
	else if (key == 'q')
	{
		double st = get_time();
		max_iters_t =  50;
		printf("run ring optimization\n");
		//run opt 
		
		g.debug_dirs.clear();
		g.debug_p.clear();

		//g.debug_p.push_back(g.halfmesh.vertexData[482].uvpos);

		vector<int> indcF;
		vector<vect3i> faces_to_opF;
		vector<vector<vect2d>> isos_to_opF;
		vector<vector<int>> boundaryF;
		vector<int> tb;
		for(int i = 0; i < g.boundaryorder.size(); i++)
		{
			tb.push_back(g.boundaryorder[i]);
		}
		boundaryF.push_back(tb);
		printf("boundary size = %d\n", g.boundaryorder.size());
		for(int i= 0 ; i < g.halfmesh.totalVerts; i++)
		{
			indcF.push_back(i);
		}
		for(int i = 0; i < g.halfmesh.totalFaces; i++)
		{
			vect3i t;
			t[0] = g.mesh.indices_tex[i][0];
			t[1] = g.mesh.indices_tex[i][1];
			t[2] = g.mesh.indices_tex[i][2];
			faces_to_opF.push_back(t);
		}
		vector<vector<int>> indc;
		vector<vector<int>> dummyboundary;dummyboundary.clear();
		vector<vector<int>> boo;
		//g.boundaryorder.clear();
		//g.get_boundary_order(g.mesh,g.mesh
		boo.push_back(g.boundaryorder);
		vector<vector<vect3i>> faces_to_op;
		vector<vector<vector<vect2d>>> isos_to_op;

		if(!param->is_full_created)
		{
			param->is_full_created = true;
			param->full_mesh = new double[g.halfmesh.totalVerts*2];
			g.perform_isometric_flattening(&g.mesh, &g.iso_tris);
			param->init_lbfgs( g.halfmesh.totalVerts );
		}

		
		for(int j = 0; j < g.halfmesh.rings.size(); j++)
		{
			printf("\rsetup ring %d", j);
			//for(int j = 0; j < i+1; j++)
			{
				vector<int> tempi;
				for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
				{
					tempi.push_back(-1);
				}
				indc.push_back(tempi);

				int ind = 0;
				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					indc[j][g.halfmesh.rings[j][ii]] = ind;
					ind++;
				}
				
				vector<vect3i> ft;
				vector<vector<vect2d>> it;

				for(int ii = 0; ii < g.halfmesh.rings_faces[j].size(); ii++)
				{
					ft.push_back( g.halfmesh.rings_faces[j][ii] );
					it.push_back( g.iso_tris[g.halfmesh.rings_find[j][ii]] );
				}

				if(j > 0)
				{
					for(int ii = 0; ii < g.halfmesh.rings_faces[j-1].size(); ii++)
					{
						ft.push_back( g.halfmesh.rings_faces[j-1][ii] );
						it.push_back( g.iso_tris[g.halfmesh.rings_find[j-1][ii]] );
					}
				}
				
				faces_to_op.push_back(ft);
				isos_to_op.push_back(it);
			
			}
		}


		//printf("18debug verts = %d, %d, %d\n",faces_to_op[18][7782][0],faces_to_op[18][7782][1],faces_to_op[18][7782][2]);
		//printf("19debug verts = %d, %d, %d\n",faces_to_op[19][7782][0],faces_to_op[19][7782][1],faces_to_op[19][7782][2]);
		////////////////////////////////////////////////////////DEBUG
		//497 482 496
		/*printf("%.10g, %.10g, %.10g\n", g.iso_tris[593][1][0],
						g.iso_tris[593][2][0],
						g.iso_tris[593][2][1]);*/
		
		for(int i = 0; i < g.halfmesh.totalFaces; i++)
		{
			g.halfmesh.faceData[i].inside_rings = false;
		}
		/*printf("vert 497\n");
		g.print_debug_info(497);

		printf("vert 492\n");
		g.print_debug_info(492);

		printf("vert 496\n");
		g.print_debug_info(496);*/
		g.halfmesh.iso_tris = g.iso_tris;


		int ttt;

		//scanf("%d", ttt);
		//int i = g.ringI;
		//int j = g.ringJ;
		double* temppos = new double[g.halfmesh.totalVerts * 2];
		double* storegrad = new double[g.halfmesh.totalVerts * 2];
		double* sol = new double[g.halfmesh.totalVerts * 2];

		for(int i = 0; i < g.halfmesh.rings.size()-1; i++)
		{
			for(int j = 0; j < i+1; j++)
			{
				printf("\r  %d, Optimize Ring %d ",i, j);
				//setup indchange
				for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
				{
					param->full_mesh[ii*2] = g.halfmesh.vertexData[ii].uvpos[0];
					param->full_mesh[ii*2+1] = g.halfmesh.vertexData[ii].uvpos[1];
				}
				
				param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
					g.halfmesh.rings[j].size(),
					&faces_to_op[j],
					&indc[j],
					j == 0 ? &boo : &dummyboundary,
					&isos_to_op[j]);
				
				for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
				{
					temppos[ii*2] = g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[0];
					temppos[ii*2+1] = g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[1];
					sol[ii*2] = 0;
					sol[ii*2+1] = 0;

				}

				int ret;
				ret = param->run(1000, temppos, sol );

				if( _finite(sol[0]) )
				{
					for(int ii = 0; ii < g.halfmesh.rings[j].size(); ii++)
					{
						g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[0] = sol[ii*2];
						g.halfmesh.vertexData[ g.halfmesh.rings[j][ii] ].uvpos[1] = sol[ii*2+1];
					}
				}

				/*param->calc_full_gradient(temppos, sol);
				for(int ii = 0; ii < 20; ii++)
				{
					printf("%.10g ", sol[ii]);
				}
				int temp = 0;
				scanf("%d",temp);*/

				//param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
				//	g.halfmesh.totalVerts,
				//	&faces_to_opF,
				//	&indcF, 
				//	//j==0 ? &g.halfmesh.col_lev[0].boundaries : &dummyboundary,
				//	//j == 0 ? &boo : &dummyboundary,
				//	&boundaryF,
				//	&g.iso_tris);

				/*for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
				{
					temppos[ii*2] = g.halfmesh.vertexData[ii].uvpos[0];
					temppos[ii*2+1] = g.halfmesh.vertexData[ii].uvpos[1];
				}*/

				//std::printf(" error %.10g", param->calc_error(temppos, 0, sol));
			/*	param->calc_full_gradient(temppos, storegrad);

				double diff = 0;
				for(int ii = 0; ii < indc[j].size(); ii++)
				{
					diff += (storegrad[indc[j][ii]*2 ] - temppos[ii*2])*(storegrad[indc[j][ii]*2 ] - temppos[ii*2]);
					diff += (storegrad[indc[j][ii]*2 + 1] - temppos[ii*2 + 1])*(storegrad[indc[j][ii]*2 + 1] - temppos[ii*2 + 1]);
				}*/

				//std::printf("diff of grads = %.10g\n", diff);

				display();
				//break;
			}
			//break;
		}
		delete[] temppos;
		delete[] sol;

		printf("Rings took = %.10g\n", get_time() - st);
	}

	else if(key == '3')
	{
		while(g.halfmesh.col_lev[0].num_col > 0)
		{
			g.halfmesh.col_lev[0].num_col--;
			
			vect2i e = g.halfmesh.col_lev[0].edge.top();
			
			//Restore to old 3D positions
			g.halfmesh.vertexData[e[0]].pos = g.halfmesh.col_lev[0].p_col.top();
			g.halfmesh.vertexData[e[1]].pos = g.halfmesh.col_lev[0].p_rem.top();

			//Delete me later... restors old uv positions
			//g.halfmesh.vertexData[e[0]].uvpos = g.halfmesh.col_lev[0].uv_p_col.top();
			//g.halfmesh.vertexData[e[1]].uvpos = g.halfmesh.col_lev[0].uv_p_rem.top();


			//Remove stack for next expansion
				//Add to expand operation later
			g.halfmesh.col_lev[0].edge.pop();
			g.halfmesh.col_lev[0].p_col.pop();
			g.halfmesh.col_lev[0].p_rem.pop();
			g.halfmesh.col_lev[0].uv_p_col.pop();
			g.halfmesh.col_lev[0].uv_p_rem.pop();


			//Calc new uvs?
			//Vertex index of collapsed
			int col = g.halfmesh.col_verts.top();
			//printf("Expand V == %d\n", col);
			g.halfmesh.vertexData[col].newindex = g.halfmesh.indchange[col];
			g.halfmesh.col_verts.pop();

			//compute dir of vertex placement
			vect2d dir;
			vect3d z;
			z[0] = 0; z[1] = 0; z[2] = 1; 

			vect3i wing = g.halfmesh.col_lev[0].wings.top();
			g.halfmesh.col_lev[0].wings.pop();
			vect3i wing2;

			bool moveotherway = false;
			if (g.halfmesh.vertexData[col].boundary)
			{ //Vertex is a boundary vertex and can only move along neighbor either boundary edge
				dir = g.halfmesh.vertexData[wing[1]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos;

				vect3d ta; ta[0] = dir[0]; ta[1] = dir[1]; ta[2] = 0;
				ta.normalize();
				vect3d norm = ta.cross(z);
				vect2d tn; tn[0] = norm[0]; tn[1] = norm[1]; //Calculates normal of edge

				wing2 = g.halfmesh.col_lev[0].wings.top();
				g.halfmesh.col_lev[0].wings.pop();

				//Check dot product of wing edge with edge norm...
				//If negative expanded vertex can cause fold
				//Move along opposit edge
				if ((g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos).dot(tn) <= 0)
				{
					moveotherway = true;
					dir = g.halfmesh.vertexData[wing2[1]].uvpos - g.halfmesh.vertexData[wing2[2]].uvpos;
				}
			}
			else
			{ //Interior collapse so find interior

				//Find perpendicular bisector of 
				vect2d a = g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos;
				vect2d b = g.halfmesh.vertexData[wing[2]].uvpos - g.halfmesh.vertexData[wing[1]].uvpos;

				vect3d ta;
				ta[0] = a[0];
				ta[1] = a[1];
				ta[2] = 0;
				ta.normalize();

				vect3d tb;
				tb[0] = b[0];
				tb[1] = b[1];
				tb[2] = 0;
				tb.normalize();

				/*vect3d c1 = z.cross(ta);
				vect3d c2 = tb.cross(z);*/

				vect3d c1 = z.cross(ta);
				vect3d c2 = z.cross(tb);

				/*vect3d c1 = ta.cross(z);
				vect3d c2 = z.cross(tb);*/

				dir[0] = (c1[0] + c2[0]) / 2.0f;
				dir[1] = (c1[1] + c2[1]) / 2.0f;
				dir.normalize();

			}

			//find position for each point
				vector<vect2i> wingtowing = g.halfmesh.col_lev[0].wing_to_wing.top();
				g.halfmesh.col_lev[0].wing_to_wing.pop();

				//loop through
				vect2d newp;
#ifdef BARY_EXPANSION
				vect2d w = g.halfmesh.col_lev[0].barys.top();
				g.halfmesh.col_lev[0].barys.pop();
				if( g.halfmesh.vertexData[col].boundary )
				{
					vector<vect2i> wingtowing2 = g.halfmesh.col_lev[0].wing_to_wing.top();
					g.halfmesh.col_lev[0].wing_to_wing.pop();

					vect2d w2 = g.halfmesh.col_lev[0].barys.top();
					g.halfmesh.col_lev[0].barys.pop();

					if (moveotherway)
					{
						w = w2;
						wing = wing2;
						wingtowing = wingtowing2;
					}
				}


				double t = 1;  //Find maximum distance along dir
				for (int j = 0; j < wingtowing.size(); j++)
				{
					double temp = find_intersection(g.halfmesh.vertexData[wing[2]].uvpos,
						g.halfmesh.vertexData[wingtowing[j][0]].uvpos,
						g.halfmesh.vertexData[wingtowing[j][1]].uvpos,
						dir);

					if (temp < t && temp > 0)
						t = temp;
				}

					
				//loop through all edges making sure not to cross
				vect2d checkp;
#ifdef ABFPP_EXPANSION

				//ABF ++ method for calculating barycentric coordinates from isometric triangles
				if( !g.halfmesh.vertexData[col].boundary )
				{
					vector<int> f;
					f.push_back( wing[0] );
					for(int j = 0; j < wingtowing.size() - 1; j++)
						f.push_back( wingtowing[j][1]);
					f.push_back(wing[1]);
					f.push_back(wing[2]);
					checkp = g.calc_new_vert(g.halfmesh.vertexData[col].pos, f, !g.halfmesh.vertexData[col].boundary);
				}
#else
				//Use from uv coordinates from befor collapsed
				vect2d checkp = g.halfmesh.vertexData[wing[0]].uvpos*w[0] +  g.halfmesh.vertexData[wing[2]].uvpos*w[1] + g.halfmesh.vertexData[wing[1]].uvpos*(1 - w[1] - w[0]);
#endif
				////////////////////////////////////////////////////////////////////////
				//Check if checkp is inside of kernel of one ring
				z[0] = 0; z[1] = 0; z[2] = 1;

				vect3d td;
				td[0] = (g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos)[0];
				td[1] = (g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos)[1];
				td[2] = 0;

				td = z.cross(td);
				td.normalize();

				vect2d offset;
				offset[0] = td[0]; offset[1] = td[1];
				bool positive = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.0000001*offset) ) >= 0;
				
				//bool positive = ((g.halfmesh.vertexData[wing[0]].uvpos+offset) - (g.halfmesh.vertexData[wing[2]].uvpos+offset)).dot( checkp - g.halfmesh.vertexData[wing[2]].uvpos ) >= 0;
				
				td[0] = (g.halfmesh.vertexData[wing[2]].uvpos - g.halfmesh.vertexData[wing[1]].uvpos)[0];
				td[1] = (g.halfmesh.vertexData[wing[2]].uvpos - g.halfmesh.vertexData[wing[1]].uvpos)[1];
				td[2] = 0;

				td = z.cross(td);
				td.normalize();
				offset[0] = td[0]; offset[1] = td[1];
				bool positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.0000001*offset) ) >= 0;
					//((g.halfmesh.vertexData[wing[2]].uvpos+offset) - (g.halfmesh.vertexData[wing[1]].uvpos+offset)).dot( checkp - g.halfmesh.vertexData[wing[2]].uvpos ) >= 0;
				
				bool keep = true;
				if( positive == positive2 )
				{ //inside wings
					//printf("ooook\n");
										
					for (int j = 0; j < wingtowing.size(); j++)
					{
						//printf("ok %d\n",j);
						td[0] = (g.halfmesh.vertexData[wingtowing[j][0]].uvpos - g.halfmesh.vertexData[wingtowing[j][1]].uvpos)[0];
						td[1] = (g.halfmesh.vertexData[wingtowing[j][0]].uvpos - g.halfmesh.vertexData[wingtowing[j][1]].uvpos)[1];
						td[2] = 0;

						td = td.cross(z);
						td.normalize();

						offset[0] = td[0]; offset[1] = td[1];

						positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wingtowing[j][0]].uvpos+.00001*offset) ) >= 0;
						//positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.000001*offset) ) >= 0;
							//((g.halfmesh.vertexData[wingtowing[j][1]].uvpos+offset) - g.halfmesh.vertexData[wingtowing[j][0]].uvpos+offset).dot( checkp - g.halfmesh.vertexData[wing[2]].uvpos ) >= 0;
						if(positive2 != positive)
						{
							//printf("shit %d\n", j);
							keep = false;
							j = wingtowing.size();
						}
					}
				}
				else
					keep = false;


				if(g.halfmesh.vertexData[col].boundary)
					keep = false;

				if(keep)
				{//point inside of offset one ring
					newp = checkp;
				}
				else if(g.halfmesh.vertexData[col].boundary)
				{

					newp = g.halfmesh.vertexData[wing[2]].uvpos + dir*(t/2.0f);
				}
				else
				{  //Check p not in kernel... project onto closest point of perpendicular bisector clamped inside kernel 
					//project
					double d = (checkp - g.halfmesh.vertexData[wing[2]].uvpos).dot(dir);
					newp = g.halfmesh.vertexData[wing[2]].uvpos + dir*min(max(t*.1, d),t*.9);
					
				}

			
				//printf("expanded p = (%f, %f)\n", newp[0], newp[1]);
				if (moveotherway)
				{
					//printf("other way\n");
					g.halfmesh.vertexData[col].uvpos = g.halfmesh.vertexData[wing[2]].uvpos;
					g.halfmesh.vertexData[wing[2]].uvpos = newp;
				}
				else
					g.halfmesh.vertexData[col].uvpos = newp;

				g.halfmesh.currentVerts++;

#ifdef PERFORM_BARY_SMOOTHING
#endif

#else

//Dont perform barycentric coordinates
//keeping out of implementation


#endif
			g.halfmesh.collapsed[col] = -1;
		}
		//else
		{
			printf("done expanding\n");
			g.fully_expanded = true;
			g.halfmesh.drawlvl = 1;
		}
	}
	else if (key == '1')
	{
		printf("finish collapsing, prepare for multi res\n");
		g.halfmesh.finish_collapsing();
		g.expand_iter = g.halfmesh.col_lev.size() - 1;
		g.halfmesh.drawlvl = g.expand_iter;
	}
	else if (key == '4')
	{
		max_iters_t =  1000000;
		printf("run optimization, full no removal\n");
		//run opt 
		
		g.debug_dirs.clear();
		g.debug_p.clear();

		if(!param->is_full_created)
		{
			param->is_full_created = true;
			param->full_mesh = new double[g.halfmesh.totalVerts*2];
			g.perform_isometric_flattening(&g.mesh, &g.iso_tris);
			param->init_lbfgs( g.halfmesh.totalVerts );
		}


		printf("num verts %d / %d",g.halfmesh.currentVerts,g.halfmesh.totalVerts );
		if(g.fully_expanded)
		{
		param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
			g.halfmesh.currentVerts,
			&g.halfmesh.col_lev[1].faces,
			&g.halfmesh.indchange, 
			&g.halfmesh.col_lev[1].boundaries,
			&g.halfmesh.col_lev[1].iso);
		}
		else
		{
			vector<vect3i> facesnew;
			for (int ii = 0; ii < g.halfmesh.totalFaces; ii++)
			{
				vect3i t;
				t = g.mesh.indices_tex[ii];
				facesnew.push_back(t);
			}
			param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/,
				g.halfmesh.currentVerts,
				&g.halfmesh.col_lev[g.expand_iter].faces,
				&g.halfmesh.indchange,
				&g.halfmesh.col_lev[g.expand_iter].boundaries,
				//&g.iso_tris);
			&g.halfmesh.col_lev[g.expand_iter].iso);
		}
		double* temppos = new double[g.halfmesh.totalVerts * 2];
		int ind = 0;
		for (int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			if (g.halfmesh.collapsed[i] == -1)
			{
				
				temppos[g.halfmesh.vertexData[i].newindex * 2] = g.halfmesh.vertexData[i].uvpos[0];
				temppos[g.halfmesh.vertexData[i].newindex * 2 + 1] = g.halfmesh.vertexData[i].uvpos[1];
				ind++;
				//printf("%f, ", temppos[ind * 2]);
			}
		}
		printf("verts in = %d\n",ind);
		double* sol = new double[g.halfmesh.totalVerts * 2];

		double start = get_time();

		int ret = -1;
		//while( ret != -998 )
		{
			ret = param->run(1000, temppos, sol );
			for(int i = 0; i < g.mesh.tex_coords.s*2; i++)
			{
				temppos[i] = sol[i];
			}
		}
		double end = get_time();

		g.timing = end - start;
		g.total_timing += g.timing;
		printf("Total time = %f, iteration = %f\n", g.total_timing, g.timing);

	//	param->print_timings();

		ind = 0;
		for (int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			//if (g.halfmesh.vertexData[i].newindex != -1)
			if (g.halfmesh.collapsed[i] == -1)
			{
				g.mesh.tex_coords[i][0] = g.halfmesh.vertexData[i].uvpos[0] = sol[g.halfmesh.vertexData[i].newindex * 2];
				g.mesh.tex_coords[i][1] = g.halfmesh.vertexData[i].uvpos[1] = sol[g.halfmesh.vertexData[i].newindex * 2 + 1];
				ind++;
			}
		}

		printf("clean up\n");
		delete[] temppos;
		delete[] sol;
		printf("clean up\n");
	}
	else if(key == 'x')
	{
		g.seamless += .0001;
		printf("seam = %f\n", g.seamless);
	}
	else if(key == 'X')
	{
		g.seamless *= 2;

		printf("seam = %f\n", g.seamless);
	}
	else if (key == '2')
	{
		max_iters_t =  10000000;
		printf("run optimization\n");
		//run opt 
		
		g.debug_dirs.clear();
		g.debug_p.clear();

		//printf("debug %d, %d, %d\n", g.halfmesh.faceData[499].vi[0],g.halfmesh.faceData[499].vi[1],g.halfmesh.faceData[499].vi[2]);
		if(!param->is_full_created)
		{
			param->is_full_created = true;
			param->full_mesh = new double[g.halfmesh.totalVerts*2];
			g.perform_isometric_flattening(&g.mesh, &g.iso_tris);
			param->init_lbfgs( g.halfmesh.totalVerts );
		}

		vector<int> indc;
		vector<vector<int>> boo; boo.push_back(g.boundaryorder);

		vector<vector<vect2d>> isos_to_opT;
		vector<vect3i> faces_to_op;
		vector<vector<vect2d>> isos_to_op;

		vector<int> remove_ind;
		vector<int> tri_inds;

		g.remove_these_tris.clear();
		printf("Remove Tris = %d\n", g.remove_these_tris.size());
		for(int i = 0; i < g.remove_these_tris.size(); i++)
		{
			for(int j = 0; j < 3; j++)
			{
				int ind = g.halfmesh.faceData[g.remove_these_tris[i]].vi[j];
			//	printf("%d\n", ind);
				bool isin = false;
				for(int ii = 0; ii < remove_ind.size(); ii++)
				{
					if(ind == remove_ind[ii])
					{
						isin = true; break;
					}
				}
				if(!isin)
				{
					//printf("%d\n",ind);
					remove_ind.push_back(ind);
				}
			}
		}
	//	printf("remove_ind = %d\n", remove_ind.size());

		int ind = 0;
		//for(int i = 0; i < g.halfmesh.totalFaces; i++)
		//{
		//	bool isin = false;
		//	for(int ii = 0; ii < g.remove_these_tris.size(); ii++)
		//	{
		//		if(i == g.remove_these_tris[ii])
		//		{
		//			isin = true; break;
		//		}
		//	}
		//	if(isin)
		//		;//tri_inds.push_back(-1);
		//	else
		//	{
		//		tri_inds.push_back(i);
		//		ind++;
		//	}
		//}

		//Get indchange
		ind = 0;
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			bool isin = false;
			for(int ii = 0; ii < remove_ind.size(); ii++)
			{
				if(i == remove_ind[ii])
				{
					isin = true; break;
				}
			}
			if(isin)
			{
				indc.push_back(-1);
				//printf("ind is -1 = %d\n", i);
			}
			else
			{
				indc.push_back(ind);
				ind++;
			}
		}

		//Get faces to opt and isos to opt
		for(int i = 0; i < g.halfmesh.totalFaces; i++)
		{
			bool isin = false;
			for(int ii = 0; ii < g.remove_these_tris.size(); ii++)
			{
				if(i == g.remove_these_tris[ii])
				{
					isin = true; break;
				}
			}
			if(!isin)
			{
				vect3i t = g.halfmesh.faceData[i].vi;
				faces_to_op.push_back(t);
				//{
				//	vector<vect2d> flat;
				//	vect2d t0, t1,t2;
				//	flat.push_back(t0);
				//	flat.push_back(t1);
				//	flat.push_back(t2);

				//	// first vertex is straight in x, with length of the edge
				//	vect3d p[3];
				//	p[0][0] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[0]].pos[0];
				//	p[0][1] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[0]].pos[1];
				//	p[0][2] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[0]].pos[2];

				//	p[1][0] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[1]].pos[0];
				//	p[1][1] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[1]].pos[1];
				//	p[1][2] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[1]].pos[2];

				//	p[2][0] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[2]].pos[0];
				//	p[2][1] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[2]].pos[1];
				//	p[2][2] = g.halfmesh.vertexData[g.halfmesh.faceData[i].vi[2]].pos[2];

				//	vect3d X, Y, Z;
				//	X = p[1] - p[0];
				//	X.normalize();
				//	Z = X % (p[2] - p[0]);
				//	Z.normalize();

				//	Y = Z % X;

				//	// store
				//	flat[0].set(0, 0);
				//	flat[1].set((p[1] - p[0]) * X, 0);
				//	flat[2].set((p[2] - p[0]) * X, (p[2] - p[0]) * Y);
				//	isos_to_opT.push_back(flat);
				//}
				//isos_to_opT.push_back(g.halfmesh.col_lev[g.expand_iter].iso[i]);

				/*vector<vect2d> tt;
				for(int ii = 0; ii < g.iso_tris[i].size(); ii++)
				{
					vect2d ttt = g.iso_tris[i][ii];
					if(!_finite(ttt[0]) || !_finite(ttt[1]))
						printf("WTF %f, %f\n",t[0],t[1] );
					tt.push_back(ttt);
				}
				isos_to_op.push_back(tt);*/
				
				//isos_to_op.push_back(g.iso_tris[i]);
			}
		}
		printf("total faces = %d, faces to opt = %d / %d\n", g.halfmesh.totalFaces, faces_to_op.size(), isos_to_opT.size());
		printf("total verts = %d, verts to opt = %d\n", g.halfmesh.totalVerts, ind);

		//set collapse level
		param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
			ind,
			&g.halfmesh.col_lev[g.expand_iter].faces,
			&indc, 
			&boo,
			&g.halfmesh.col_lev[g.expand_iter].iso);

		param->set_seamless_param(g.seamless);

		//Find matching boundaries
		param->match_boundary.clear();
		//param->match_v.clear();
		vector<int> inserted;
		inserted.clear();
		for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			MatchingBoundary t;

			bool skipme = false;
			for(int j = 0; j < inserted.size(); j++)
			{
				if(i == inserted[j])
				{
					skipme = true;
				}
			}
			if(skipme)
				continue;
			int match1 = i;
			int match2 = (i+1)%g.halfmesh.boundaries[0].size();
			for(int j = 0; j < g.halfmesh.boundaries[0].size(); j++)
			{
				//check first point
				//if(i != (j-1)%g.halfmesh.boundaries[0].size())
				{
					double len = (g.halfmesh.boundaries[0][i]->pos - g.halfmesh.boundaries[0][j]->pos).length();
					double len2 = (g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->pos - g.halfmesh.boundaries[0][(j-1)%g.halfmesh.boundaries[0].size()]->pos).length();


					if((i == 20 && len < .0001) || i == 20 && len2 < .0001 )
						printf("checker");
					if(len < 0.00000000001 && len2 < .000000000001)
					{
						if(i == j && (i+1)%g.halfmesh.boundaries[0].size() == (j-1)%g.halfmesh.boundaries[0].size())
							printf("wtf seams\n");

						if(i == (j-1)%g.halfmesh.boundaries[0].size() && (i+1)%g.halfmesh.boundaries[0].size() == j)
							printf("wtf seams\n");

						match1 = j;
						match2 = (j-1)%g.halfmesh.boundaries[0].size();
						break;
					}
					//match1 = len < .000001 ? j : match1;

					////check second point
					//if(match1 == j && (i+1)%g.halfmesh.boundaries[0].size() != (j+1)%g.halfmesh.boundaries[0].size())
					//{
					//	double len = (g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->pos - g.halfmesh.boundaries[0][(j+1)%g.halfmesh.boundaries[0].size()]->pos).length();
					//	match2 = len < .0000001 ? (j+1)%g.halfmesh.boundaries[0].size() : match2;
					//}

					////check second point
					//if(match1 == j && (i+1)%g.halfmesh.boundaries[0].size() != (j-1)%g.halfmesh.boundaries[0].size())
					//{
					//	double len = (g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->pos - g.halfmesh.boundaries[0][(j-1)%g.halfmesh.boundaries[0].size()]->pos).length();
					//	match2 = len < .0000001 ? (j-1)%g.halfmesh.boundaries[0].size() : match2;
					//}
				}
				
				
			}

			/*t.u1 = g.halfmesh.boundaries[0][i]->index;
			t.u2 = g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->index;

			t.v1 = g.halfmesh.boundaries[0][ match1 ]->index;
			t.v2 = g.halfmesh.boundaries[0][ match2 ]->index;*/

			t.u1 = g.halfmesh.boundaries[0][i]->index;
			t.v1 = g.halfmesh.boundaries[0][(i+1)%g.halfmesh.boundaries[0].size()]->index;

			t.u2 = g.halfmesh.boundaries[0][ match1 ]->index;
			t.v2 = g.halfmesh.boundaries[0][ match2 ]->index;

			/*if(t.u1 == t.u2)
				printf("wtf seam 2\n");*/

			bool insert = true;
			if(t.u1 == t.u2)
				insert = false;
			/*for(int j = 0; j < param->match_boundary.size(); j++)
			{
				if( t.u1 == param->match_boundary[j].u1 && t.u2 == param->match_boundary[j].u2 )
				{
					insert = false;
					break;
				}
				if( t.u1 == param->match_boundary[j].u2 && t.u2 == param->match_boundary[j].u1 )
				{
					insert = false;
					break;
				}

				if( t.v1 == param->match_boundary[j].u1 && t.v2 == param->match_boundary[j].u2 )
				{
					insert = false;
					break;
				}
				if( t.v1 == param->match_boundary[j].u2 && t.v2 == param->match_boundary[j].u1 )
				{
					insert = false;
					break;
				}
			}*/
			if(insert)
			{
				param->match_boundary.push_back(t);
				inserted.push_back(i);
				inserted.push_back(match2);
			}
		}

		printf("matching boundaries = %d\n", param->match_boundary.size());

		for(int ii = 0; ii < g.halfmesh.totalVerts; ii++)
		{
			param->full_mesh[ii*2] = g.halfmesh.vertexData[ii].uvpos[0];
			param->full_mesh[ii*2+1] = g.halfmesh.vertexData[ii].uvpos[1];
		}

		double* temppos = new double[ g.halfmesh.totalVerts * 2];
		double* sol = new double[ g.halfmesh.totalVerts * 2];

		ind = 0;
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			if(indc[i] != -1)
			{
				temppos[2*ind] = g.halfmesh.vertexData[i].uvpos[0];
				temppos[2*ind+1] = g.halfmesh.vertexData[i].uvpos[1];
				ind++;
			}
		}

		//if(g.fully_expanded)
		//{
		//param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
		//	g.halfmesh.currentVerts,
		//	&g.halfmesh.col_lev[1].faces,
		//	&g.halfmesh.indchange, 
		//	&g.halfmesh.col_lev[1].boundaries,
		//	&g.halfmesh.col_lev[1].iso);
		//}
		//else
		//{
		//	vector<vect3i> facesnew;
		//	for (int ii = 0; ii < g.halfmesh.totalFaces; ii++)
		//	{
		//		vect3i t;
		//		t = g.mesh.indices_tex[ii];
		//		facesnew.push_back(t);
		//	}
		//	param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/,
		//		g.halfmesh.currentVerts,
		//		&g.halfmesh.col_lev[g.expand_iter].faces,
		//		&g.halfmesh.indchange,
		//		&g.halfmesh.col_lev[g.expand_iter].boundaries,
		//		//&g.iso_tris);
		//	&g.halfmesh.col_lev[g.expand_iter].iso);
		//}
		//double* temppos = new double[g.halfmesh.totalVerts * 2];
		//ind = 0;
		//for (int i = 0; i < g.halfmesh.totalVerts; i++)
		//{
		//	if (g.halfmesh.collapsed[i] == -1)
		//	{
		//		
		//		temppos[g.halfmesh.vertexData[i].newindex * 2] = g.halfmesh.vertexData[i].uvpos[0];
		//		temppos[g.halfmesh.vertexData[i].newindex * 2 + 1] = g.halfmesh.vertexData[i].uvpos[1];
		//		ind++;
		//		//printf("%f, ", temppos[ind * 2]);
		//	}
		//}
		/*printf("verts in = %d\n",ind);
		double* sol = new double[g.halfmesh.totalVerts * 2];*/

		double start = get_time();

		int ret = -1;
		//while( ret != -998 )
		{
			ret = param->run(1000, temppos, sol );

		//	if(ret != -1001)
			for(int i = 0; i < ind*2; i++)
			{
				temppos[i] = sol[i];
			}
		}
		double end = get_time();

		g.timing = end - start;
		g.total_timing += g.timing;
		printf("Total time = %f, iteration = %f\n", g.total_timing, g.timing);

	//	param->print_timings();

		//ind = 0;
		//for (int i = 0; i < g.halfmesh.totalVerts; i++)
		//{
		//	//if (g.halfmesh.vertexData[i].newindex != -1)
		//	if (g.halfmesh.collapsed[i] == -1)
		//	{
		//		g.mesh.tex_coords[i][0] = g.halfmesh.vertexData[i].uvpos[0] = sol[g.halfmesh.vertexData[i].newindex * 2];
		//		g.mesh.tex_coords[i][1] = g.halfmesh.vertexData[i].uvpos[1] = sol[g.halfmesh.vertexData[i].newindex * 2 + 1];
		//		ind++;
		//	}
		//}


		ind = 0;
		double movement = 0;
		for(int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			if(indc[i] != -1)
			{
				movement += (g.halfmesh.vertexData[i].uvpos[0] - temppos[2*ind])*(g.halfmesh.vertexData[i].uvpos[0] - temppos[2*ind]);
				movement += (g.halfmesh.vertexData[i].uvpos[1] - temppos[2*ind+1])*(g.halfmesh.vertexData[i].uvpos[1] - temppos[2*ind+1]);
				g.mesh.tex_coords[i][0] = g.halfmesh.vertexData[i].uvpos[0] = temppos[2*ind];
				g.mesh.tex_coords[i][1] = g.halfmesh.vertexData[i].uvpos[1] = temppos[2*ind+1];
				ind++;
			}
		}

		printf("movement = %.15g\n", movement);
		if(movement < REMOVE_DIS && g.remove_tri_index != -1)
		{
			/*printf("remove triangle = %d\n", tri_inds[g.remove_tri_index] );
			g.remove_these_tris.push_back( tri_inds[g.remove_tri_index] );*/

			printf("remove triangle = %d\n", g.remove_tri_index );
			g.remove_these_tris.push_back( g.remove_tri_index );
		}


		double min = 9999999999999;
		double max = -999999;
		double tot = 0;
		//glBegin(GL_LINES);
		for(int i = 0; i < param->match_boundary.size(); i++)
		{
			double err = param->match_boundary[i].p1*param->match_boundary[i].p2*param->match_boundary[i].p3*param->match_boundary[i].p4;

			tot += err;
			min = err < min ? err : min;
			max = err > max ? err : max;
		}
		
		printf("seam error = %.10g, min = %.10g, max = %.10g\n", tot, min,max);

		//printf("clean up\n");
		delete[] temppos;
		delete[] sol;
		//printf("clean up\n");
	}
	else if (key == '#')
	{

		if (g.expand_iter > 0)
		{
			g.debug_dirs.clear();
			g.debug_p.clear();
			printf("expand mesh 1 level\n");
			//Expand operation

			//for (int i = 0; i < g.halfmesh.col_lev[g.expand_iter].num_col; i++)
			if(g.halfmesh.col_lev[0].num_col > 0)
			{
				g.halfmesh.col_lev[0].num_col--;
				//setup new index for newly made points
				int col = g.halfmesh.col_verts.top();
				g.halfmesh.vertexData[col].newindex = g.halfmesh.indchange[col];
				g.halfmesh.col_verts.pop();

				
				vect3i wing = g.halfmesh.col_lev[0].wings.top();
				g.halfmesh.col_lev[0].wings.pop();

				vect3i wing2;
			

				//compute dir of vertex placement
				vect2d dir;
				vect3d z;
				z[0] = 0; z[1] = 0; z[2] = 1; 

				bool bound = false;
				bool moveotherway = false;
				if (g.halfmesh.vertexData[col].boundary)
				{
					bound = true;
					//vect3d ta;
					//printf("%d, %d\n", wing[1], wing[2]);
					dir = g.halfmesh.vertexData[wing[1]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos;

					vect3d ta; ta[0] = dir[0]; ta[1] = dir[1]; ta[2] = 0;
					ta.normalize();
					vect3d norm = ta.cross(z);
					vect2d tn; tn[0] = norm[0]; tn[1] = norm[1];

					wing2 = g.halfmesh.col_lev[0].wings.top();
					g.halfmesh.col_lev[0].wings.pop();

					//Check dot product of wing edge with edge norm...
					//If negative expanded vertex can cause fold
					//Move along opposit edge
					if ((g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos).dot(tn) <= 0)
					{
						moveotherway = true;
						dir = g.halfmesh.vertexData[wing2[1]].uvpos - g.halfmesh.vertexData[wing2[2]].uvpos;

						//Use wing2 to calc direction
					}
				}
				else
				{
					vect2d a = g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos;
					vect2d b = g.halfmesh.vertexData[wing[2]].uvpos - g.halfmesh.vertexData[wing[1]].uvpos;

					vect3d ta;
					ta[0] = a[0];
					ta[1] = a[1];
					ta[2] = 0;
					ta.normalize();

					vect3d tb;
					tb[0] = b[0];
					tb[1] = b[1];
					tb[2] = 0;
					tb.normalize();

					/*vect3d c1 = z.cross(ta);
					vect3d c2 = tb.cross(z);*/

					vect3d c1 = z.cross(ta);
					vect3d c2 = z.cross(tb);

					/*vect3d c1 = ta.cross(z);
					vect3d c2 = z.cross(tb);*/

					dir[0] = (c1[0] + c2[0]) / 2.0f;
					dir[1] = (c1[1] + c2[1]) / 2.0f;
					dir.normalize();

				}
				
				//find position for each point
				vector<vect2i> wingtowing = g.halfmesh.col_lev[0].wing_to_wing.top();
				g.halfmesh.col_lev[0].wing_to_wing.pop();

				//loop through
				vect2d newp;
#ifdef BARY_EXPANSION

				vect2d w = g.halfmesh.col_lev[0].barys.top();
				g.halfmesh.col_lev[0].barys.pop();
				if(bound)
				{
					vector<vect2i> wingtowing2 = g.halfmesh.col_lev[0].wing_to_wing.top();
					g.halfmesh.col_lev[0].wing_to_wing.pop();

					vect2d w2 = g.halfmesh.col_lev[0].barys.top();
					g.halfmesh.col_lev[0].barys.pop();

					if (moveotherway)
					{
						w = w2;
						wing = wing2;
						wingtowing = wingtowing2;
					}
				}


				double t = 1;
				for (int j = 0; j < wingtowing.size(); j++)
				{
					double temp = find_intersection(g.halfmesh.vertexData[wing[2]].uvpos,
						g.halfmesh.vertexData[wingtowing[j][0]].uvpos,
						g.halfmesh.vertexData[wingtowing[j][1]].uvpos,
						dir);

					if (temp < t && temp > 0)
						t = temp;
				}

					
				//loop through all edges making sure not to cross
				vect2d checkp = g.halfmesh.vertexData[wing[0]].uvpos*w[0] +  g.halfmesh.vertexData[wing[2]].uvpos*w[1] + g.halfmesh.vertexData[wing[1]].uvpos*(1 - w[1] - w[0]);

#ifdef ABFPP_EXPANSION
				if(!bound)
				{
					vector<int> f;
					f.push_back( wing[0] );
					for(int j = 0; j < wingtowing.size() - 1; j++)
						f.push_back( wingtowing[j][1]);
					f.push_back(wing[1]);
					f.push_back(wing[2]);
					checkp = g.calc_new_vert(g.halfmesh.vertexData[col].pos, f, !bound);

					vect2d tempp = g.halfmesh.vertexData[wing[0]].uvpos*w[0] +  g.halfmesh.vertexData[wing[2]].uvpos*w[1] + g.halfmesh.vertexData[wing[1]].uvpos*(1 - w[1] - w[0]);
					printf("");


				}
#endif

				//is in one ring if so keep it...
				
				z[0] = 0; z[1] = 0; z[2] = 1;

				vect3d td;
				td[0] = (g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos)[0];
				td[1] = (g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos)[1];
				td[2] = 0;

				td = z.cross(td);
				td.normalize();

				vect2d offset;
				offset[0] = td[0]; offset[1] = td[1];
				bool positive = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.0000001*offset) ) >= 0;
				
				//bool positive = ((g.halfmesh.vertexData[wing[0]].uvpos+offset) - (g.halfmesh.vertexData[wing[2]].uvpos+offset)).dot( checkp - g.halfmesh.vertexData[wing[2]].uvpos ) >= 0;
				
				td[0] = (g.halfmesh.vertexData[wing[2]].uvpos - g.halfmesh.vertexData[wing[1]].uvpos)[0];
				td[1] = (g.halfmesh.vertexData[wing[2]].uvpos - g.halfmesh.vertexData[wing[1]].uvpos)[1];
				td[2] = 0;

				td = z.cross(td);
				td.normalize();
				offset[0] = td[0]; offset[1] = td[1];
				bool positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.0000001*offset) ) >= 0;
					//((g.halfmesh.vertexData[wing[2]].uvpos+offset) - (g.halfmesh.vertexData[wing[1]].uvpos+offset)).dot( checkp - g.halfmesh.vertexData[wing[2]].uvpos ) >= 0;
				
				bool keep = true;
				if( positive == positive2 )
				{ //inside wings
					//printf("ooook\n");
										
					for (int j = 0; j < wingtowing.size(); j++)
					{
						//printf("ok %d\n",j);
						td[0] = (g.halfmesh.vertexData[wingtowing[j][0]].uvpos - g.halfmesh.vertexData[wingtowing[j][1]].uvpos)[0];
						td[1] = (g.halfmesh.vertexData[wingtowing[j][0]].uvpos - g.halfmesh.vertexData[wingtowing[j][1]].uvpos)[1];
						td[2] = 0;

						td = td.cross(z);
						td.normalize();

						offset[0] = td[0]; offset[1] = td[1];

						positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wingtowing[j][0]].uvpos+.00001*offset) ) >= 0;
						//positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.000001*offset) ) >= 0;
							//((g.halfmesh.vertexData[wingtowing[j][1]].uvpos+offset) - g.halfmesh.vertexData[wingtowing[j][0]].uvpos+offset).dot( checkp - g.halfmesh.vertexData[wing[2]].uvpos ) >= 0;
						if(positive2 != positive)
						{
							//printf("shit %d\n", j);
							keep = false;
							j = wingtowing.size();
						}
					}
				}
				else
					keep = false;

				if(bound)
					keep = false;
				if(keep)
				{//point inside of offset one ring
					newp = checkp;
					g.debug_p.push_back(checkp);

					vect2d tt;
					tt[0] = (g.halfmesh.vertexData[wing[2]].uvpos)[0];
					tt[1] = (g.halfmesh.vertexData[wing[2]].uvpos)[1];
					vector<vect2d> pushme;
					pushme.push_back(tt);
					vect2d tt2;
					tt2[0] = (g.halfmesh.vertexData[wing[0]].uvpos)[0];
					tt2[1] = (g.halfmesh.vertexData[wing[0]].uvpos)[1];
					pushme.push_back(tt2);
					g.debug_dirs.push_back(pushme);

					vect2d att;
					att[0] = (g.halfmesh.vertexData[wing[2]].uvpos)[0];
					att[1] = (g.halfmesh.vertexData[wing[2]].uvpos)[1];
					vector<vect2d> apushme;
					apushme.push_back(att);
					vect2d att2;
					att2[0] = (g.halfmesh.vertexData[wing[1]].uvpos)[0];
					att2[1] = (g.halfmesh.vertexData[wing[1]].uvpos)[1];
					apushme.push_back(att2);
					g.debug_dirs.push_back(apushme);

					//if(g.debug_p.size() == 2)
					//{
					//	printf("");
					//	for (int j = 0; j < wingtowing.size(); j++)
					//	{
					//		//printf("ok %d\n",j);
					//		td[0] = (g.halfmesh.vertexData[wingtowing[j][0]].uvpos - g.halfmesh.vertexData[wingtowing[j][1]].uvpos)[0];
					//		td[1] = (g.halfmesh.vertexData[wingtowing[j][0]].uvpos - g.halfmesh.vertexData[wingtowing[j][1]].uvpos)[1];
					//		td[2] = 0;

					//		td = td.cross(z);
					//		td.normalize();

					//		offset[0] = td[0]; offset[1] = td[1];

					//	//	positive2 = offset.dot( checkp - (g.halfmesh.vertexData[wing[2]].uvpos+.00001*offset) ) >= 0;
					//		g.debug_p.push_back(g.halfmesh.vertexData[wingtowing[j][0]].uvpos+.001*offset  );
					//		g.debug_p.push_back(g.halfmesh.vertexData[wingtowing[j][1]].uvpos+.001*offset  );
					//	}
					//}

				//	printf("hurray\n");
				}
				else if(bound)
				{

					newp = g.halfmesh.vertexData[wing[2]].uvpos + dir*(t/2.0f);
				}
				else
				{
					//project
					double d = (checkp - g.halfmesh.vertexData[wing[2]].uvpos).dot(dir);
					newp = g.halfmesh.vertexData[wing[2]].uvpos + dir*min(max(t*.1, d),t*.9);
					
				}

			
				if (moveotherway)
				{
					//printf("other way\n");
					g.halfmesh.vertexData[col].uvpos = g.halfmesh.vertexData[wing[2]].uvpos;
					g.halfmesh.vertexData[wing[2]].uvpos = newp;
				}
				else
					g.halfmesh.vertexData[col].uvpos = newp;
				g.halfmesh.collapsed[col] = -1;
				
				g.halfmesh.currentVerts++;
#else

//#ifdef ABFPP_EXPANSION
//
//	// Compute new position for uv
//		//bring back 3d position
//		//compute weights in 3D
//		//compute new position
//	
//	//Check if position inside of kernel
//		//check halfspaces
//
//#ifdef ONE_RING_SMOOTHING
//	//Optional:
//	//Do one ring smoothing
//		//basically step one
//#endif

//
//#else
//#endif

				if (bound)
				{
					vector<vect2i> wingtowing2 = g.halfmesh.col_lev[g.expand_iter].wing_to_wing.top();
					g.halfmesh.col_lev[g.expand_iter].wing_to_wing.pop();
					//loop through all edges making sure not to cross

					if (moveotherway)
						wingtowing = wingtowing2;
					double t = 1;
					for (int j = 0; j < wingtowing.size(); j++)
					{
						//g.debug_p.push_back(g.halfmesh.vertexData[wingtowing[j][0]].uvpos);
						//g.debug_p.push_back(g.halfmesh.vertexData[wingtowing[j][1]].uvpos);

						/*vector<vect2d> temp1;
						temp1.push_back(g.halfmesh.vertexData[wingtowing[j][0]].uvpos);
						temp1.push_back(g.halfmesh.vertexData[wingtowing[j][1]].uvpos);
						g.debug_dirs.push_back(temp1);*/

						double temp = find_intersection(g.halfmesh.vertexData[wing[2]].uvpos,
							g.halfmesh.vertexData[wingtowing[j][0]].uvpos,
							g.halfmesh.vertexData[wingtowing[j][1]].uvpos,
							dir);

						if (temp < t && temp > 0)
							t = temp;
					}

				//	printf("%d checks -> t for bound = %f\n", wingtowing.size(), t);
					newp = (t*.5)*dir + g.halfmesh.vertexData[wing[2]].uvpos;
					//g.debug_p.push_back(newp);
					/*g.debug_p.push_back(g.halfmesh.vertexData[wing[0]].uvpos);
					g.debug_p.push_back(g.halfmesh.vertexData[wing[1]].uvpos);
					g.debug_p.push_back(g.halfmesh.vertexData[wing[2]].uvpos);*/
					//newp = g.halfmesh.vertexData[wing[2]].uvpos;
					/*vector<vect2d> temp;
					temp.push_back(g.halfmesh.vertexData[wing[2]].uvpos);
					temp.push_back(newp);
					g.debug_dirs.push_back(temp);*/
				}
				else
				{
					if (wingtowing.size() == 0)
					{

						//use length to wing edge
						newp = (g.halfmesh.vertexData[wing[0]].uvpos - g.halfmesh.vertexData[wing[2]].uvpos).length() * dir + g.halfmesh.vertexData[wing[2]].uvpos;
					}
					else
					{
						//loop through all edges making sure not to cross
						double t = 1;
						for (int j = 0; j < wingtowing.size(); j++)
						{
							double temp = find_intersection(g.halfmesh.vertexData[wing[2]].uvpos,
								g.halfmesh.vertexData[wingtowing[j][0]].uvpos,
								g.halfmesh.vertexData[wingtowing[j][1]].uvpos,
								dir);

							if (temp < t && temp > 0)
								t = temp;
						}
						newp = (t / 2.0f) * dir + g.halfmesh.vertexData[wing[2]].uvpos;
					}

					/*vector<vect2d> temp;
					temp.push_back(g.halfmesh.vertexData[wing[2]].uvpos);
					temp.push_back(dir + g.halfmesh.vertexData[wing[2]].uvpos);

					vector<vect2d> temp1;
					temp1.push_back(g.halfmesh.vertexData[wing[0]].uvpos);
					temp1.push_back(g.halfmesh.vertexData[wing[2]].uvpos);

					g.debug_dirs.push_back(temp);
					g.debug_dirs.push_back(temp1);

					if (wing[1] != -1)
					{

						vector<vect2d> temp2;
						temp2.push_back(g.halfmesh.vertexData[wing[1]].uvpos);
						temp2.push_back(g.halfmesh.vertexData[wing[2]].uvpos);
						g.debug_dirs.push_back(temp2);
					}*/

				}
				
				//g.halfmesh.vertexData[col].newindex = g.halfmesh.currentVerts;
				//g.halfmesh.indchange[col] = g.halfmesh.currentVerts;

				if (moveotherway)
				{
					//printf("other way\n");
					g.halfmesh.vertexData[col].uvpos = g.halfmesh.vertexData[wing[2]].uvpos;
					g.halfmesh.vertexData[wing[2]].uvpos = newp;
				}
				else
					g.halfmesh.vertexData[col].uvpos = newp;
				g.halfmesh.collapsed[col] = -1;
				
				g.halfmesh.currentVerts++;
				//printf("%d, ", col);
#endif
			}
			else
			{
				printf("done expanding\n");
				g.halfmesh.drawlvl = 1;
			}

			//increase index

			printf("num verts = %d\n", g.halfmesh.currentVerts);
			//g.expand_iter--;
			//g.halfmesh.drawlvl--;
		}
		else
			printf("Done expanding\n");
	}
	else if(key == 'f')
	{
		for(int i = 0; i < 100; i++)
			g.halfmesh.collapseMinimal();
	}
	else if(key == 'd')
	{
		if(g.display_error)
		{
			g.display_error = false;
		}
		else
		{
			printf("Display error\n");
			g.display_error = true;
			calc_display_error();
			
		}
	}
	else if ( key == 'c' )
	{
		if (!g.draw_halfmesh)
		{
			g.build_half_mesh();
			g.draw_halfmesh = true;
		}
		else
		{
			int count = 100;
			//g.halfmesh.collapse_one_iteration(g.mesh.indices_tex);
			if(first_collapse)
			{
				first_collapse = false;
				
				CollapseLevel lvl;
				lvl.num_col = 0;
				g.halfmesh.col_lev.push_back(lvl);

				g.halfmesh.gen_collapse_lvl();
			}
			while(count > 0)
			{
			g.halfmesh.collapseMinimal();
			count--;
			}
		}
		g.fully_expanded = false;
		printf("\r collapse to  = %d", g.halfmesh.currentVerts);
			
	}
	else if( key == 'C')
	{
		printf("Finish collapse and recompute floaters\n");

		g.halfmesh.gen_collapse_lvl();
#ifdef PERFORM_FLOATERS
		for(int ii = 0; ii < g.halfmesh.boundaries.size(); ii++)
		{
			g.halfmesh.tuttes_embedding_chart(ii);
		}

#endif
		g.halfmesh.scale_halmesh(100);

		for (int i = 0; i < g.halfmesh.totalVerts; i++)
		{
			g.mesh.positions[i] = g.halfmesh.vertexData[i].pos;
			g.mesh.tex_coords[i] = g.halfmesh.vertexData[i].uvpos;
		}
		printf("num levels debug %d\n", g.halfmesh.col_lev.size());
		g.halfmesh.col_lev.clear();
		g.halfmesh.gen_collapse_lvl();

		g.halfmesh.analyze_boundary(); //Boundaries
		
		g.halfmesh.err.clear();
		
		double tot = 0;
		double maxi = -99999;
		for(int i = 0; i < g.halfmesh.totalFaces;i++)
		{
			double t = g.tri_error( g.halfmesh.col_lev[0].iso[i][1][0],
				g.halfmesh.col_lev[0].iso[i][2][0],
				g.halfmesh.col_lev[0].iso[i][2][1],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos[0],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[0] ].uvpos[1],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos[0],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[1] ].uvpos[1],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos[0],
				g.halfmesh.vertexData[ g.halfmesh.faceData[i].vi[2] ].uvpos[1]
				) ;
				
				t = t > 10000 && false ? 10000 : t;
			g.halfmesh.err.push_back( .2*-1*log( 1.0f / (t+1)));
			tot += t;
			if(t > maxi)
				maxi = t;
		}

		printf("max err = %f.. avg = %f\n", maxi, tot/g.halfmesh.totalFaces);
		for(int i = 0; i < g.halfmesh.totalFaces;i++)
		{

		}

		g.halfmesh.bin_errors.clear();
		g.halfmesh.bins.clear();
		for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			g.halfmesh.bin_errors.push_back(0);
			g.halfmesh.bins.push_back(0);
		}
		vect2d s;
		s[0] = 1;s[1] = 0;
		for(int i = 0; i < g.halfmesh.boundaries[0].size(); i++)
		{
			g.halfmesh.bins[i] = g.find_ang(s, g.halfmesh.boundaries[0][i]->uvpos);
		}
	}
	else if(key == 'g')
	{
		if(!g.draw_grad)
		{
			//printf("compute gradients, not working right now\n");

			//param->set_collapse_level(g.halfmesh.totalVerts/* - g.halfmesh.col_lev[g.halfmesh.col_lev.size() - 1].num_col*/, 
			//g.halfmesh.currentVerts,
			//&g.halfmesh.col_lev[g.expand_iter].faces,
			//&g.halfmesh.indchange, 
			//&g.halfmesh.col_lev[g.expand_iter].boundaries,
			//&g.halfmesh.col_lev[g.expand_iter].iso);

			g.draw_grad = true;
			double* temppos = new double[g.mesh.tex_coords.s * 2];
			for(int i = 0; i < g.mesh.tex_coords.s; i++)
			{
				temppos[i*2] = g.halfmesh.vertexData[i].uvpos[0];
				temppos[i*2+1] = g.halfmesh.vertexData[i].uvpos[1];
			}
			//param->ind_tex = g.mesh.indices_tex;
			double* sol = new double[g.mesh.tex_coords.s * 2];
			double* nsol = new double[g.mesh.tex_coords.s * 2];

			param->calc_full_gradient(temppos, sol);

			for(int i = 0; i < g.mesh.tex_coords.s*2; i++)
			{
				nsol[i] = -1*sol[i];
				g.gradDraw[i] = 100*nsol[i];

			}
//				g.gradDraw[i] = nsol[i] = -1*sol[i];

			double f;
			//g.max_param = param->max_parameter(&f, temppos, sol, nsol, g.mesh.tex_coords.s * 2 );
			//printf("Max param = %.10g\n", g.max_param);
	//		g.max_param = 0;

			delete[] temppos;
			delete[] sol;
			delete[] nsol;


		//	g.draw_grad = false;
		}
		else
		{
			g.draw_grad = false;
		}
	}
	printf("before out\n");
	glutPostRedisplay();
}

void motion(int x, int y)
{
	y = SCREEN_HEIGHT - y;

	if(g.mleft && false)
	{
		float dx = (x - g.mprevX);
			float dy = (y - g.mprevY);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glRotated(dx, 0, 1, 0);
			glRotated(dy, -1, 0, 0);
			glMultMatrixd(g.rotation);
			glGetDoublev(GL_MODELVIEW_MATRIX, g.rotation);
	}
	else if (g.mmiddle || g.mleft)
	{
		double dx = (x - g.mprevX);
		double dy = (y - g.mprevY);

		vect3d X, Y;
		X[0] = g.rotation[0];
		X[1] = g.rotation[4];
		X[2] = g.rotation[8];
		Y[0] = g.rotation[1];
		Y[1] = g.rotation[5];
		Y[2] = g.rotation[9];

		g.focus -= X*(5e-4*g.zoom*dx) + Y*(5e-4*g.zoom*dy);

	}
	else if (g.mright)
	{
		double dy = (y - g.mprevY);
		g.zoom -= g.zoom*5e-3*dy;

		if(g.zoom < .001)
			g.zoom = .001;

		if(g.zoom > MAX_ZOOM)
			g.zoom = MAX_ZOOM;
	}
	
	g.mprevX = x;
	g.mprevY = y;

	g.mhomeX = x;
	g.mhomeY = y;

	glutPostRedisplay();
}

void click(int button, int state, int x, int y)
{
	y = SCREEN_HEIGHT - y;

	// Mouse state that should always be stored on pressing
	if (state == GLUT_DOWN)
	{
		g.mhomeX = x;
		g.mhomeY = y;
		g.mprevX = x;
		g.mprevY = y;
	}

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		g.mleft = true;
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		g.mleft = false;
	}

	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		g.mright = true;
	}
	if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
	{
		g.mright = false;
	}
	
	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
	{
		g.mmiddle = true;
	}
	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_UP)
	{
		g.mmiddle = false;
	}

	glutPostRedisplay();
}

void load_mesh(string name)
{
	printf("Loading mesh %s\n", name.c_str());

#ifdef LOAD_OBJ_JOETYPE
	load_obj2(name + ".obj", g.mesh);
#else
	load_obj(name + ".obj", g.mesh);
#endif
	g.charts.merge(g.mesh);
	g.get_boundary_order(&g.mesh,&g.charts,  &g.boundaryorder);
	printf("boundary size = %d\n", g.boundaryorder.size());
}


void init_gl()
{

}


int main(int argc, char **argv)
{
#ifdef _WIN64
	printf("x64\n");
#else
#ifdef _WIN32
	printf("x86\n");
#endif
#endif

	if(argc == 4)
	{
		string meshname = argv[1]; 
		load_mesh(meshname);
	
		printf("Vertices Loaded = %d\n", g.mesh.positions.s);
		printf("Texture V Loaded = %d\n", g.mesh.tex_coords.s);
		printf("Boundary Edges = %d\n", g.boundaryorder.size());

		g.init();

		g.charts.merge(g.mesh);

		g.build_half_mesh();
#ifdef FLIP_UVS
		g.halfmesh.flip_tex();
#endif
		g.draw_halfmesh = true;
		vector<vector<int>> t;
		t.push_back(g.boundaryorder);
		param = new JasonFull("jason method\n", g.mesh.tex_coords.s, t, g.iso_tris);

		switch( atoi(argv[2]) ){
		
		case 1:
			g.opt_epsilon = 1e-1;
			break;
			case 2:
			g.opt_epsilon = 1e-2;
			break;
			case 3:
			g.opt_epsilon = 1e-3;
			break;
			case 4:
			g.opt_epsilon = 1e-4;
			break;
			case 5:
			g.opt_epsilon = 1e-5;
			break;
			case 6:
			g.opt_epsilon = 1e-6;
			break;
			case 7:
			g.opt_epsilon = 1e-7;
			break;
			case 8:
			g.opt_epsilon = 1e-8;
			break;
			case 9:
			g.opt_epsilon = 1e-9;
			break;
			case 10:
			g.opt_epsilon = 1e-10;
			break;
			case 11:
			g.opt_epsilon = 1e-11;
			break;
			case 12:
			g.opt_epsilon = 1e-12;
			break;
			case 13:
			g.opt_epsilon = 1e-13;
			break;
			case 14:
			g.opt_epsilon = 1e-14;
			break;
			case 15:
			g.opt_epsilon = 1e-15;
			break;
			case 16:
			g.opt_epsilon = 1e-16;
			break;
			case 17:
			g.opt_epsilon = 1e-17;
			break;
			case 18:
			g.opt_epsilon = 1e-18;
			break;
			case 19:
			g.opt_epsilon = 1e-19;
			break;
			case 20:
			g.opt_epsilon = 1e-20;
			break;			
		}
		
		g.do_boundary = atoi( argv[3] );

		g.takefrom = 1.0;

		run_full_optimization();
	}
	else
	{
		string meshname = MESHNAME; 
		load_mesh(SAVETO + meshname + "/" + meshname);
	
		printf("Vertices Loaded = %d\n", g.mesh.positions.s);
		printf("Texture V Loaded = %d\n", g.mesh.tex_coords.s);
		printf("Boundary Edges = %d\n", g.boundaryorder.size());

		g.init();

		g.charts.merge(g.mesh);
		g.build_half_mesh();
	#ifdef FLIP_UVS
		g.halfmesh.flip_tex();
	#endif
	
		g.draw_halfmesh = true;

		vector<vector<int>> t;
		t.push_back(g.boundaryorder);
		param = new JasonFull("jason method\n", g.mesh.tex_coords.s, t, g.iso_tris);

		g.opt_epsilon = 1e-4;
		g.do_boundary = true;

		g.takefrom = 1.0;
	
	}

	g.do_boundary = true;

	printf("Running main loop\n");
	////Open Gl
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(SCREEN_WDITH, SCREEN_HEIGHT);
	glutCreateWindow("No Fold Parameterization");

	//// display
	display();

	//// Glut Setup
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(input);
	glutMouseFunc(click);
	glutMotionFunc(motion);

	glutMainLoop();

	if (true)
		system("pause");

	return 0;
	
}