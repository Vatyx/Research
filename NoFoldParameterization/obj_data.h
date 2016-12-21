#pragma once

#include "array.h"
#include "vect.h"
#include <string>

using namespace std;

struct MeshData
{
	// data
	Array<vect3d> positions;
	Array<vect2d> tex_coords;
	Array<vect3d> normals;

	// triangles' indices into data
	Array<vect3i> indices_pos;
	Array<vect3i> indices_tex;
	Array<vect3i> indices_norm;
};

void load_obj(string filename, MeshData &mesh);
void load_obj2(string filename, MeshData &mesh);
void load_off(string filename, MeshData &mesh);