#include "mesh_charts.h"

MeshChartVertex* MeshChart::find(MeshChartVertex* e)
{
	if (e->parent != e)
		e->parent = find(e->parent);
	return e->parent;
}

void MeshChart::merge(MeshChartVertex* x, MeshChartVertex* y)
{
	MeshChartVertex* x_root = find(x);
	MeshChartVertex* y_root = find(y);

	if (x_root == y_root)
		return;

	if (x_root->rank < y_root->rank)
		x_root->parent = y_root;
	else if (y_root->rank < x_root->rank)
		y_root->parent = x_root;
	else
	{
		y_root->parent = x_root;
		x_root->rank++;
	}
}

void MeshChart::make_set(MeshChartVertex* e)
{
	e->parent = e;
}

void MeshChart::merge(MeshData &mesh)
{
	verts.clear();
	verts.resize(mesh.tex_coords.s);
	
	// create all of the combined vertices from the triangle data.
	for (int it = 0; it < mesh.indices_pos.s; it++)
	{
		vect3i tt = mesh.indices_tex[it];
		vect3i pp = mesh.indices_pos[it];

		for (int iv = 0; iv < 3; iv++)
		{
			MeshChartVertex mv(tt[iv], pp[iv]);
			verts[tt[iv]] = mv;
		}
	}

	// set pointers to self (make sets for all verts)
	for (int i = 0; i < verts.s; i++)
	{
		make_set(&verts[i]);
	}

	// combine vertices from triangle connectivity
	for (int it = 0; it < mesh.indices_pos.s; it++)
	{
		vect3i tt = mesh.indices_tex[it];

		for (int iv = 0; iv < 3; iv++)
		{
			int ivn = (iv + 1) % 3;

			MeshChartVertex &mv0 = verts[tt[iv]];
			MeshChartVertex &mv1 = verts[tt[ivn]];
			merge(&mv0, &mv1);
		}
	}

	// build list of vertices in charts
	for (int i = 0; i < verts.s; i++)
	{
		int r = find(&verts[i])->index_mesh;
		verts[i].index_chart = chart_data[r].verts.s;
		chart_data[r].verts.push_back(i);
	}

	// build list of triangles in chart
	for (int  i = 0; i < mesh.indices_tex.s; i++)
	{
		int r = find(&verts[mesh.indices_tex[i][0]])->index_mesh; // get root of any verts to find chart
		chart_data[r].tris.push_back(i);
	}

	// find lowest manhatan distance pinned vertices for each chart
	{
		// construct pinned map
		map<int, vect2i> pinned;
		for (int i = 0; i < verts.s; i++)
		{
			int r = find(&verts[i])->index_mesh;

			if (pinned.count(r) == 0)
			{
				// always best if the first
				pinned[r] = i;
			}
			else
			{
				// find lowest and highest
				vect2i &v = pinned[r];
				vect2d p = mesh.tex_coords[i];
				vect2d l = mesh.tex_coords[v[0]];
				vect2d h = mesh.tex_coords[v[1]];

				float pt = p[0] + p[1];
				float lt = l[0] + l[1];
				float ht = h[0] + h[1];

				if (pt < lt)
					v[0] = i;
				if (pt > ht)
					v[1] = i;
			}
		}

		// copy into chart_data map, which also initializes the chart_data structure for the mesh!!!!
		for (map<int, vect2i>::iterator it = pinned.begin(); it != pinned.end(); ++it)
		{
			root_index.push_back(it->first);
			chart_data[it->first].pinned = it->second;
#ifdef PIN_ONE_VERTEX
		// if pinning one vertex, choose something at random (such as the root) that is more likely to be in the middle
			chart_data[it->first].pinned[0] = it->first;
#endif
		}

	}

	// count edges
	for (int it = 0; it < mesh.indices_tex.s; it++)
	{
		for (int ie = 0; ie < 3; ie++)
		{
			int ien = (ie + 1) % 3;
			vect2i edge;
			edge[0] = mesh.indices_tex[it][ie];
			edge[1] = mesh.indices_tex[it][ien];
			if (edge[1] < edge[0])
				swap(edge[0], edge[1]);
			edge_count[edge]++;
		}
	}

	// set boundary flags
	for (map<vect2i, int>::iterator it = edge_count.begin(); it != edge_count.end(); ++it)
	{
		if (it->second == 1)
		{
			verts[it->first.v[0]].is_boundary = true;
			verts[it->first.v[1]].is_boundary = true;
		}
	}
}

float det(vect2d x, vect2d y)
{
	return (x[0]*y[1] - x[1]*y[0]) * .5f;
}

void MeshChart::calc_areas(MeshData &mesh, int idx)
{
	// init to zero
	for (map<int, MeshChartData>::iterator it = chart_data.begin(); it != chart_data.end(); ++it)
	{
		it->second.area[idx] = 0;
	}

	for (int it = 0; it < mesh.indices_tex.s; it++)
	{
		vect3i tt = mesh.indices_tex[it];

		vect2d &v0 = mesh.tex_coords[tt[0]];
		vect2d &v1 = mesh.tex_coords[tt[1]];
		vect2d &v2 = mesh.tex_coords[tt[2]];

		double a = abs(det(v1 - v0, v2 - v0));

		int r = find(&verts[tt[0]])->index_mesh;

		chart_data[r].area[idx] += a;
	}
}
