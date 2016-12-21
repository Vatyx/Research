#ifndef HALF_EDGE_MESH_H
#define HALF_EDGE_MESH_H

//#include <vector>
#include <stack>
//#include "vect.h"
#include "QEF.h"

#include <map>
#include <set>

#include "TrackingHeap.h"
#include "obj_data.h"
#include "mesh_charts.h"


#define BOUNDARY_THICKNESS 6.0f

#define MIN_BOUNDARY_SIZE 5
//#define SCOTT_COLLAPSE
#define PI 3.14159265359
#define NUM_TYPE double
//#define QEF_TYPE QEFQR
#define QEF_TYPE QEFCenter
//#define QEF_TYPE QEFNormal

#define QEF_VERTEX

class Mesh;
class ExpandOp;

using namespace std;

class HE_Vertex;
class HE_Edge;
class HE_Face;



class HE_HalfEdge
{
public:
	HE_Vertex *v; // points to the next vertex
	HE_Face *f;
	HE_HalfEdge *flip, *next, *prev;
	HE_Edge *e;
	
};

class HE_Vertex 
{
public:
	HE_HalfEdge *he; // the edge that comes into the vertex v->he->v == v

	vect3d pos;
	vect2d uvpos;

	QEF_TYPE<NUM_TYPE, 3> qef;

	bool collapsable;
	bool boundary;
	bool recompute_qef;

	int recalc;

	int which_boundary;

	int newindex;
	int index;

	int chart_ind;
};

class HE_Edge : public TrackingHeap::TrackingHeapData
{
public:
	HE_HalfEdge *he;
	// could remove this info, it's just for speeding up algorithm
	bool dirty;
	int index;
	//VEC_TYPE pos;
	vect3d pos;
	HE_Edge ( void ) /*: TrackingHeap::TrackingHeapData ( )*/ {}
};

class HE_Face
{
public:
	int removed;
	bool hole; // assume no holes in mesh
	int index;
	HE_HalfEdge *he;
	bool inside_rings;
	vect3i vi;
	
};

struct HE_ExpandOp
{
	// all of the objects that are removed in an edge collapse
	HE_Edge *e[2], *ec; // the edges from the collapsed wings, and the edge that was collapsed
	HE_HalfEdge *hec[2][3]; // the half edges from the collapsed wings
	HE_Vertex *v; // the removed vertex
	HE_Face *f[2]; // the faces that were deleted

	HE_HalfEdge *he[2]; // kept half edges that determine collapse
	
	vect3d pos[2]; // positions of the collapsed vertices
	vect2d uvpos[2];

};

class CollapseLevel
{
public:
	int num_col;
	//where to put current vertices into next level for opt
	//vector<int> vert_ind_change; 

	//just face indicies
	vector<vect3i> faces;
	vector<vector<vect2d>> iso;
	vector<vector<int>> boundaries;
	vector<vector<int>> onerings;

	stack<vect3d> p_col;
	stack<vect3d> p_rem;
	stack<vect2d> uv_p_col;
	stack<vect2d> uv_p_rem;
	stack<vect2i> edge;

	stack<vect2d> barys;
	stack<vect3i> wings;
	stack<vector<vect2i>> wing_to_wing;

	vector<HE_ExpandOp> expands;
};


class PeakInterval
{
public:
	int start, end, peak;
	float vs, ve, vp;
	bool circ;

	PeakInterval()
	{
		start = end = peak = 0;
		vs =ve = vp = 0;
		circ = false;
	}

	void join( PeakInterval j )
	{
		if(j.start <  start)
			start = j.start;

		if(j.end > end)
			end = j.end;

		if(j.vp > vp)
		{
			vp = j.vp;
			peak = j.peak;
		}
	}
};

class HalfEdgeMesh
{
friend class MultiresMesh;
public:

	double time_col, time_init, time_floater;
	int collapse_iter;
	int drawlvl;

	stack<int> col_verts;
	vector<int> newindex;
	vector<int> indchange;
	vector<double> err;
	vector<int> peaks;

	HE_ExpandOp *expandOps;
	vector<HE_Vertex> vertexData;
	vector<HE_Face> faceData;
	vector<HE_Edge> edgeData;
	vector<HE_HalfEdge> heData;
	vector<int> collapsed;
	vector<double> bin_errors;
	vector<double> bins;
	int totalFaces;
	int totalVerts;

	vector<PeakInterval> pi;

	vect2d zero_param;
	vect2d pie2_param;
	vector<vect2d> error_view;
	vector<int> error_ind;
	vector<vector<vect2d>> iso_tris;
	vector<bool> inlist;
	vector<bool> finlist;
	vector<vector<int>> rings;
	vector<vector<vect3i>> rings_faces;
	vector<vector<int>> rings_find;
	vector<vector<vector<vect2d>>> rings_isos;

	TrackingHeap errQueue;

	vector<vector<HE_Vertex*>> boundaries;
	vector<int> current_boundary_size;
	int chart_boundary_ind;

	vector<vect2d> debugedge;

	vector<CollapseLevel> col_lev;
	void HalfEdgeMesh::analyze_boundary();
	void perform_floaters();

	int currentVerts;
	void HalfEdgeMesh::print_debug_info(int ind);
	void build_rings();
	void flip_tex();
	void finish_collapsing();
	void collapseMinimal ( void );
	bool isSafeCollapse ( HE_Edge *e );

	void tuttes_embedding_chart(int c);
	/*float minimizeEdge ( HE_Edge *e );
	void buildQEFs ( void );
	void combineQEFs ( HE_Vertex *v1, HE_Vertex *v2 );*/
	ExpandOp *collapse ( HE_Edge *e, int* &head, int &f1, int &f2, int &v );
	
	void draw_3Dfull();
	void scale_halmesh(double rad);
	void gen_collapse_lvl();
	void set_to_recalc_qef( HE_Vertex *v);
	void ring_and_spoke( HE_Vertex *v );
	void collapse_one_iteration ( Array<vect3i> indices_pos);

	bool isboundary(HE_HalfEdge *he );
	bool isboundary(HE_Vertex *he );

#ifdef SCOTT_COLLAPSE
	HE_HalfEdge* find_collapsable_edge_boundary(HE_Vertex v);
	bool collapse ( HE_HalfEdge *e );
#else
	bool collapse ( HE_Edge *e, HE_Vertex *v );
	HE_Edge* find_collapsable_edge_boundary(HE_Vertex v);
#endif

	#ifdef SCOTT_COLLAPSE
	void set_one_ring_uncollapsable(HE_Vertex *v);
#else
	void sanity();
	void set_one_ring_uncollapsable(HE_Vertex v);
	bool is_3_valence( HE_Vertex *v );
#endif

	void make_collapsable();


	float minimizeEdge ( HE_Edge *e );
	void buildQEFs ( void );
	void print_buildQEFs ( HE_Vertex *v );
	void combineQEFs ( HE_Vertex *v1, HE_Vertex *v2 );

	void save_mesh();
	void save_mesh2();
	void setup_boundaries();
	void setup_first_col_lvl();
public:
	HalfEdgeMesh ( void )
	{
		/*totalFaces = totalVerts = currentVerts = 0;
		vertexData = NULL;
		faceData = NULL;
		edgeData = NULL;
		heData = NULL;*/
		expandOps = 0;
		mesh_bookmark = 0;
	}

	~HalfEdgeMesh ( void )
	{
		if ( totalFaces + totalVerts > 0 )
		{
			/*delete[] vertexData;
			delete[] faceData;
			delete[] edgeData;
			delete[] heData;*/
			//delete[] expandOps;
		}
		//delete mesh_bookmark;
	}

	void init ( MeshData *mesh, MeshChart*charts );

	void draw(void);
	void draw_full(void);
	void draw_col_lvl(void);

	//ExpandOp *collapseMinimal ( int* &head, int &f1, int &f2, int &v );
	//void collapseMinimal ( void );
	void expand();

	int getNumVerts ( void ) { return currentVerts; }

	
	// bookmarking
	HalfEdgeMesh *mesh_bookmark;
	void addBookmark()
	{
		//printf("currentVerts = %d\n", currentVerts);
		delete mesh_bookmark;
		mesh_bookmark = new HalfEdgeMesh();
		copyTo(mesh_bookmark);
	}

	void restoreBookmark()
	{
		mesh_bookmark->copyTo(this);
	}
	
	void copyTo(HalfEdgeMesh *m)
	{
		m->currentVerts = currentVerts;
		m->edgeData = edgeData;
		m->faceData = faceData;
		m->heData = heData;
		m->totalFaces = totalFaces;
		m->totalVerts = totalVerts;
		m->vertexData = vertexData;
		//m->expandOps = expandOps;
	}
};


#endif