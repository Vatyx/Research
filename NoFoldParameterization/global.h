#include <windows.h>
#include <stdio.h>

//#include "ParamMethod.h"
#include <vector>

#include "vect.h"
#include "obj_data.h"
#include "mesh_charts.h"

#include "PointGrid.h"
#include "ResizeArray.h"

#include "HalfEdge.h"

#include <GL/glut.h>

struct flapOpSave
{
	vect3i verts;
	double alpha;
	double L;
};

struct Global
{
	Global();
	void init();


	/*********************************************/
	//GUI vars
	/*********************************************/
	int ringI,ringJ;
	//Display Code
	vect3d bg_color; // background color

	double opt_epsilon;
	bool do_boundary;

	double rotation[16];
	vect3d focus;
	double zoom;
	bool printInfo;

	GLint viewport[4]; // View orientation
	GLdouble modelview[16];	
	GLdouble proj[16];

	//Mouse input
	bool mleft, mright, mmiddle;
	double mhomeX, mhomeY;
	double mprevX, mprevY;
	bool draw_full_mesh;
	void init_draw();
	
	bool draw_grad, draw_both_grads, draw_uv01;
	bool display_error;
	bool draw_halfmesh;
	double *err;

	bool multiple_charts;
	int init_type;
	int remove_tri_index;
	vector<int> remove_these_tris;
	/*********************************************/
	//prog vars
	/*********************************************/
	bool debug_seamless;
	int which_seamless;


	//Input mesh and tex coordinates
	MeshData mesh;
	MeshData smesh;
	//Helper class for mesh charts
	MeshChart charts, scharts;

	HalfEdgeMesh halfmesh;

	//Isometric triangles
	vector<vector<vect2d>> iso_tris;
	vector<vector<vect2d>> siso_tris;

	//Ordering of boundary vertices
	vector<int> boundaryorder;
	vector<int> sbo;

	bool fully_expanded;
	bool draw_simplified;
	//Gradients drawn on screen after computing them
	double *gradDraw;
	double *globalboundarygrad;
	double *interiorgrad;

	double max_param;
	int max_param_ind;
	double seamless;

	int iteration_count;
	double timing;
	double total_timing;

	int expand_iter;

	double takefrom;

	vector<vector<vect2d>> debug_dirs;
	vector<vect2d> debug_p;

	void print_debug_info(int ind);

	void perform_multi_resolution();
	/*Store result in iso_tris*/
	void perform_isometric_flattening(MeshData *m, vector<vector<vect2d>> *isot);

	/*Store result in boundary order*/
	void get_boundary_order(MeshData *m, MeshChart *mc, vector<int> *bo);

	void build_half_mesh();

	void perform_vertex_smoothing(int v_ind);

	/*Normalize mesh size and center*/
	void center_mesh();
	void normalize_mesh();

	void flip_uv( MeshData *m );

	//Output Mesh
	void save_mesh(string name);

	//outputjoes
	void save_mesh2(string name);

	//Floater mean value parameterization
	void perform_floater_param(MeshData *m, MeshChart *c, vector<int> *bo);

	void convert_halfmesh();

	double find_ang(vect2d v, vect2d w);

	double tri_error(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2);
	vect2d calc_new_vert( vect3d p, vector<int> f, bool interior );

	vector<flapOpSave> flapLoopup;
};


extern Global g;