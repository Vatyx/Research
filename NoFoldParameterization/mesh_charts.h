#pragma once

// written by Josiah Manson
// this code uses union-find to determine what chart a vertex belongs to.
// all of the vertices in a chart are in a tree, with one tree per chart.
// each chart is represented by a single vertex that is the root of the tree for the chart.

#include "vect.h"
#include "array.h"
#include "obj_data.h"
#include <map>


using namespace std;

struct MeshChartVertex
{
	MeshChartVertex *parent;
	int index_mesh; // index in the mesh texture
	int index_mesh_pos; // index in the mesh positions
	int index_chart; // index in the chart
	int rank;
	bool is_boundary;

	MeshChartVertex()
	{
		index_mesh = -1;
		index_mesh_pos = -1;
		index_chart = -1;
		parent = 0;
		rank = 0;
		is_boundary = false;
	}
	MeshChartVertex(int vi, int vip)
	{
		index_mesh = vi;
		index_mesh_pos = vip;
		index_chart = -1;
		parent = 0;
		rank = 0;
		is_boundary = false;
	}
};

struct MeshChartData
{
	Array<int> verts; // index of texture vert in mesh
	Array<int> tris; // triangle indexing chart
	vect2i pinned;
	double area[2]; // area [1] before and [2] after reparameterization
};

struct MeshChart
{
	Array<MeshChartVertex> verts;

	map<int, MeshChartData> chart_data; // root of chart -> data
	Array<int> root_index; // index of the root of the ith chart
	map<vect2i, int> edge_count; // count the number of trinagles adjecent to an edge

	bool is_pinned_arap(int idx)
	{
		for (map<int, MeshChartData>::iterator it = chart_data.begin(); it != chart_data.end(); ++it)
		{
#ifdef PIN_ONE_VERTEX
			if (it->second.pinned[0] == idx)
#else
			if (it->second.pinned[0] == idx || it->second.pinned[1] == idx)
#endif
				return true;
		}
		return false;
	}
	
	bool is_pinned_lscm(int idx)
	{
		for (map<int, MeshChartData>::iterator it = chart_data.begin(); it != chart_data.end(); ++it)
		{
			if (it->second.pinned[0] == idx || it->second.pinned[1] == idx)
				return true;
		}
		return false;
	}

	int chart_num(int i)
	{
		MeshChartVertex *root = find(&verts[i]);
		int c = 0;
		for (map<int, MeshChartData>::iterator it = chart_data.begin(); it != chart_data.end(); ++it, c++)
		{
			if (it->first == root->index_mesh)
				return c;
		}
		assert(false); // better have found it
		return 0;
	}

	MeshChartVertex* find(MeshChartVertex* e);
	void merge(MeshChartVertex* x, MeshChartVertex* y);
	void make_set(MeshChartVertex* e);

	void merge(MeshData &mesh);

	void calc_areas(MeshData &mesh, int idx);
};
