#include "global.h"
#include "timer.h"

#include "JasonFull.h"

ParamMethod *par;
void Global::perform_multi_resolution()
{
	//Init data structure to full size of mesh
	double *pos = new double[halfmesh.totalVerts];
	double *sol = new double[halfmesh.totalVerts];

	//For number of collapse levels
	for (int i = 0; i < halfmesh.col_lev.size(); i++)
	{
		//Init iteration
		//need to set connectivity
		//param->ind_tex = g.mesh.indices_tex;

		//Run optimization
		double start = get_time();

		int ret = -1;
		//while( ret != -998 )
		{
			ret = par->run(1000, pos, sol );
			for(int i = 0; i < g.mesh.tex_coords.s*2; i++)
			{
				pos[i] = sol[i];
			}
		}
		double end = get_time();

		g.timing = end - start;
		g.total_timing += g.timing;

		//Expand back out
	}

	//Set solution

	delete[] pos;
	delete[] sol;
}