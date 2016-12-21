#include "HalfEdge.h"

#include <set>
#include <map>
#include <algorithm>
#include <assert.h>
#include <float.h>
#include "vect.h"
#include <queue>

#include "SparseMatrix.h"

#include <GL/glut.h>

void HalfEdgeMesh::analyze_boundary()
{
	double tot_length = 0;
	double max_length = -9999;
	double min_length = 9999;
	double root = 0;
	for(int i = 0; i < this->boundaries[0].size(); i++)
	{
		
		double t = ((this->boundaries[0][i]->pos) - (this->boundaries[0][(i+1)%this->boundaries[0].size()]->pos)).length();
			
			tot_length += t;
			root += ((this->boundaries[0][i]->pos) - (this->boundaries[0][(i+1)%this->boundaries[0].size()]->pos)).length2();

			max_length = t > max_length ? t : max_length;
			min_length = t < min_length ? t : min_length;

	}

	printf("Boundaries = %d\n", boundaries[0].size());
	printf("tot = %f,   avg = %f\n", tot_length, tot_length /  boundaries[0].size());
	printf("tot = %f,   avg = %f\n", root, sqrt(root /  boundaries[0].size()));
	printf("max = %f,   min = %f\n", max_length, min_length);
}
void HalfEdgeMesh::scale_halmesh(double rad)
{
	double currad = 0;
	double maxx, minx;

	minx = vertexData[0].uvpos[0];
	maxx = vertexData[0].uvpos[0];
	for (int i = 1; i < totalVerts; i++)
	{
		minx = min(minx, vertexData[i].uvpos[0]);
		maxx = max(maxx, vertexData[i].uvpos[0]);
	}

	double tarea2 = 0;
	for (int it = 0; it < totalFaces; it++)
	{
		//if (!faceData[it].hole)
		{
			//QUESTION.... SHOULD I ADD HOLE POLYS
			vect<3, vect2d> flat;

			// first vertex is straight in x, with length of the edge
			vect3d p[3];
			p[0][0] = faceData[it].he->v->uvpos[0];
			p[0][1] = faceData[it].he->v->uvpos[1];
			p[0][2] = 0;

			p[1][0] = faceData[it].he->next->v->uvpos[0];
			p[1][1] = faceData[it].he->next->v->uvpos[1];
			p[1][2] = 0;

			p[2][0] = faceData[it].he->next->next->v->uvpos[0];
			p[2][1] = faceData[it].he->next->next->v->uvpos[1];
			p[2][2] = 0;

			vect3d X, Y, Z;
			X = p[1] - p[0];
			X.normalize();
			Z = X % (p[2] - p[0]);
			Z.normalize();

			Y = Z % X;

			// store
			flat[0].set(0, 0);
			flat[1].set((p[1] - p[0]) * X, 0);
			flat[2].set((p[2] - p[0]) * X, (p[2] - p[0]) * Y);

			tarea2 += flat[1][0] * flat[2][1];
		}
	}

	double uvarea = 0;
	for(int i = 0; i < totalFaces; i++)
	{
		vect3d a;
		a[0] = vertexData[faceData[i].vi[2]].uvpos[0] - vertexData[faceData[i].vi[0]].uvpos[0];
		a[1] = vertexData[faceData[i].vi[2]].uvpos[1] - vertexData[faceData[i].vi[0]].uvpos[1];
		a[2] = 0;

		vect3d b;
		b[0] =  vertexData[faceData[i].vi[1]].uvpos[0] - vertexData[faceData[i].vi[0]].uvpos[0];
		b[1] = vertexData[faceData[i].vi[1]].uvpos[1] - vertexData[faceData[i].vi[0]].uvpos[1];
		b[2] = 0;

		uvarea += abs(a.cross(b)[2]);
	}

	printf("3d a = %.10g, uv a  = %.10g\n", tarea2, uvarea);
	//#ifdef PERFORM_SCALE
	currad = abs(maxx - minx) / 2.0f;

	/*for (int i = 0; i < totalVerts; i++)
	{
		vertexData[i].uvpos *= rad / currad;
	}*/

	double curs = 3.14159265359 * currad*currad;
	double prevs = 3.14159265359 * rad*rad;

	printf("%.10g , %.10g\n", curs, prevs);
	printf("%.10g\n", sqrt(prevs / curs));

	for (int i = 0; i < totalVerts; i++)
	{
		//vertexData[i].pos *= (rad / currad);
		 vertexData[i].pos *= sqrt(prevs / curs);
		 vertexData[i].uvpos *= sqrt(prevs / tarea2) ;
	}

	printf("cur = %.10g\n", currad);
	minx = vertexData[0].uvpos[0];
	maxx = vertexData[0].uvpos[0];
	for (int i = 1; i < totalVerts; i++)
	{
		minx = min(minx, vertexData[i].uvpos[0]);
		maxx = max(maxx, vertexData[i].uvpos[0]);
	}
	

	currad = abs(maxx - minx) / 2.0f;
	printf("cur = %.10g\n", currad);



	double tarea = 0;
	for (int it = 0; it < totalFaces; it++)
	{
		//if (!faceData[it].hole)
		{
			//QUESTION.... SHOULD I ADD HOLE POLYS
			vect<3, vect2d> flat;

			// first vertex is straight in x, with length of the edge
			vect3d p[3];
			p[0] = faceData[it].he->v->pos;
			p[1] = faceData[it].he->next->v->pos;
			p[2] = faceData[it].he->next->next->v->pos;

			vect3d X, Y, Z;
			X = p[1] - p[0];
			X.normalize();
			Z = X % (p[2] - p[0]);
			Z.normalize();

			Y = Z % X;

			// store
			flat[0].set(0, 0);
			flat[1].set((p[1] - p[0]) * X, 0);
			flat[2].set((p[2] - p[0]) * X, (p[2] - p[0]) * Y);

			tarea += flat[1][0] * flat[2][1];
		}
	}
	printf("3D = %.10g\n", tarea);

	tarea2 = 0;
	for (int it = 0; it < totalFaces; it++)
	{
		//if (!faceData[it].hole)
		{
			//QUESTION.... SHOULD I ADD HOLE POLYS
			vect<3, vect2d> flat;

			// first vertex is straight in x, with length of the edge
			vect3d p[3];
			p[0][0] = faceData[it].he->v->uvpos[0];
			p[0][1] = faceData[it].he->v->uvpos[1];
			p[0][2] = 0;

			p[1][0] = faceData[it].he->next->v->uvpos[0];
			p[1][1] = faceData[it].he->next->v->uvpos[1];
			p[1][2] = 0;

			p[2][0] = faceData[it].he->next->next->v->uvpos[0];
			p[2][1] = faceData[it].he->next->next->v->uvpos[1];
			p[2][2] = 0;

			vect3d X, Y, Z;
			X = p[1] - p[0];
			X.normalize();
			Z = X % (p[2] - p[0]);
			Z.normalize();

			Y = Z % X;

			// store
			flat[0].set(0, 0);
			flat[1].set((p[1] - p[0]) * X, 0);
			flat[2].set((p[2] - p[0]) * X, (p[2] - p[0]) * Y);

			tarea2 += flat[1][0] * flat[2][1];
		}
	}
	printf("UV = %.10g\n", tarea2);

#ifdef PERFORM_DIFF_SCALE
	for (int i = 0; i < totalVerts; i++)
	{
		//vertexData[i].pos *= (rad / currad);
		 vertexData[i].pos /= sqrt(tarea / tarea2);
		// vertexData[i].uvpos *= sqrt(prevs / tarea2) ;
	}
	tarea = 0;
	for (int it = 0; it < totalFaces; it++)
	{
		//if (!faceData[it].hole)
		{
			//QUESTION.... SHOULD I ADD HOLE POLYS
			vect<3, vect2d> flat;

			// first vertex is straight in x, with length of the edge
			vect3d p[3];
			p[0] = faceData[it].he->v->pos;
			p[1] = faceData[it].he->next->v->pos;
			p[2] = faceData[it].he->next->next->v->pos;

			vect3d X, Y, Z;
			X = p[1] - p[0];
			X.normalize();
			Z = X % (p[2] - p[0]);
			Z.normalize();

			Y = Z % X;

			// store
			flat[0].set(0, 0);
			flat[1].set((p[1] - p[0]) * X, 0);
			flat[2].set((p[2] - p[0]) * X, (p[2] - p[0]) * Y);

			tarea += flat[1][0] * flat[2][1];
		}
	}
	printf("3D = %.10g\n", tarea);
#endif
	//tarea / tarea2
}
bool is_3_val(HE_HalfEdge *he)
{
	HE_Vertex *wv = he->next->v;

	HE_HalfEdge *start = wv->he;
	HE_HalfEdge *now = wv->he;

	int count = 0;
	do
	{
		if (now->flip == NULL)
		{
			return false;
		}

		now = now->flip->next->next;
		count++;
	} while (start != now && count <1000);
	if (count == 1000)
		printf("SHIT\n");
	if (count == 3)
		return true;
	return false;
}

#define AREA_TOL 0.00000001
bool is_0_area(HE_HalfEdge * he)
{
	vect2d p = he->v->uvpos;
	vect3d pp = he->v->pos;
	 
	vect2d a,b;
	vect3d aa,bb;
	if(he->next->flip != NULL)
	{
		//check top
		a = p - he->next->flip->next->v->uvpos;
		b = he->next->v->uvpos - he->next->flip->next->v->uvpos;
		if(abs( -1*a[1]*b[0] + a[0]*b[1]) < AREA_TOL)
		{
			return true;
		}

		aa = he->next->flip->next->v->pos;
		bb = he->next->v->pos;
		if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
			return true;
		if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
			return true;
	}

	if(he->next->next->flip != NULL)
	{
		//check bot
		a = p - he->next->next->flip->next->v->uvpos;
		b = he->next->v->uvpos - he->next->next->flip->next->v->uvpos;
		if(abs( -1*a[1]*b[0] + a[0]*b[1]) < AREA_TOL)
		{
			return true;
		}

		aa = he->next->next->flip->next->v->pos;
		bb = he->next->v->pos;
		if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
			return true;
		if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
			return true;
	}

	
	if(he->flip != NULL)
	{//check other wing
		HE_HalfEdge *hef = he->flip;
		if(hef->next->flip != NULL)
		{
			//check top
			a = p - hef->next->flip->next->v->uvpos;
			b = hef->next->v->uvpos - hef->next->flip->next->v->uvpos;
			if(abs( -1*a[1]*b[0] + a[0]*b[1])  < AREA_TOL)
			{
				return true;
			}

			aa = hef->next->flip->next->v->pos;
			bb = hef->next->v->pos;
			if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
				return true;
			if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
				return true;
		}

		if(hef->next->next->flip != NULL)
		{
			//check bot
			a = p - hef->next->next->flip->next->v->uvpos;
			b = hef->next->v->uvpos - hef->next->next->flip->next->v->uvpos;
			if(abs( -1*a[1]*b[0] + a[0]*b[1]) < AREA_TOL)
			{
				return true;
			}

			aa = hef->next->next->flip->next->v->pos;
			bb = hef->next->v->pos;
			if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
				return true;
			if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
				return true;
		}
	}
	return false;
}

bool is_0_areaF(HE_HalfEdge * he)
{
	vect2d p = he->v->uvpos;
	vect3d pp = he->v->pos;
	 
	vect2d a,b;
	vect3d aa,bb;
	if(he->next->flip != NULL)
	{
		//check top
		a = p - he->next->flip->next->v->uvpos;
		b = he->next->v->uvpos - he->next->flip->next->v->uvpos;
		if(abs( -1*a[1]*b[0] + a[0]*b[1]) < AREA_TOL)
		{
			return true;
		}

		aa = he->next->flip->next->v->pos;
		bb = he->next->v->pos;
		if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
			return true;
		if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
			return true;
	}

	if(he->next->next->flip != NULL)
	{
		//check bot
		a = p - he->next->next->flip->next->v->uvpos;
		b = he->next->v->uvpos - he->next->next->flip->next->v->uvpos;
		if(abs( -1*a[1]*b[0] + a[0]*b[1]) < AREA_TOL)
		{
			return true;
		}

		aa = he->next->next->flip->next->v->pos;
		bb = he->next->v->pos;
		if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
			return true;
		if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
			return true;
	}

	
	if(he->flip != NULL)
	{//check other wing
		HE_HalfEdge *hef = he->flip;
		if(hef->next->flip != NULL)
		{
			//check top
			a = p - hef->next->flip->next->v->uvpos;
			b = hef->next->v->uvpos - hef->next->flip->next->v->uvpos;
			if(abs( -1*a[1]*b[0] + a[0]*b[1])  < AREA_TOL)
			{
				return true;
			}

			aa = hef->next->flip->next->v->pos;
			bb = hef->next->v->pos;
			if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
				return true;
			if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
				return true;
		}

		if(hef->next->next->flip != NULL)
		{
			//check bot
			a = p - hef->next->next->flip->next->v->uvpos;
			b = hef->next->v->uvpos - hef->next->next->flip->next->v->uvpos;
			if(abs( -1*a[1]*b[0] + a[0]*b[1]) < AREA_TOL)
			{
				return true;
			}

			aa = hef->next->next->flip->next->v->pos;
			bb = hef->next->v->pos;
			if(aa[0] == pp[0] && aa[1] == pp[1] && aa[2] == pp[2])
				return true;
			if(bb[0] == pp[0] && bb[1] == pp[1] && bb[2] == pp[2])
				return true;
		}
	}
	return false;
}

void HalfEdgeMesh::flip_tex()
{
	for (int i = 0; i < totalVerts; i++)
	{
		vertexData[i].uvpos[0] *= -1;
	}
}
// debug code
bool insideV ( HE_Vertex *v, HE_HalfEdge *he )
{
	HE_HalfEdge *curr = v->he;

	do
	{
		if ( he == curr )
		{
			return true;
		}
		curr = curr->flip->next;
	} while ( curr != v->he );
	return false;
}

void HalfEdgeMesh::make_collapsable()
{
	for(int i = 0; i < vertexData.size(); i++)
	{
		vertexData[i].collapsable = true;
	}
}

void HalfEdgeMesh::setup_boundaries()
{
	//Build Boundary
	for (int i = 0; i < totalFaces; i++)
	{
		HE_HalfEdge *he = (faceData[i].he);
		if (he->flip == NULL)
		{
			he->v->boundary = true;
			he->next->next->v->boundary = true;
		}

		he = he->next;
		if (he->flip == NULL)
		{
			he->v->boundary = true;
			he->next->next->v->boundary = true;
		}

		he = he->next;
		if (he->flip == NULL)
		{
			he->v->boundary = true;
			he->next->next->v->boundary = true;
		}
	}

	//Setup multi boundaries
	boundaries.clear();
	current_boundary_size.clear();
	int bound_ind = 0;
	for (int i = 0; i < totalVerts; i++)
	{
		
		if ( vertexData[i].which_boundary < 0 && vertexData[i].boundary )
		{
			vector<HE_Vertex*> temp;
			//start wrapping around for boundary
			
			//find he on boundary pointing away from v
			vertexData[i].which_boundary = bound_ind;
			temp.push_back(&vertexData[i]);
			HE_HalfEdge *he = vertexData[i].he;
			HE_HalfEdge *now = he->next;

			while (now->flip != NULL)
			{
				now = now->flip->next;
			}

			//now = proper half edge on boundary
			he = now;
			while (he->v != &vertexData[i])
			{
				if (!he->v->boundary)
					printf("shit\n");
				he->v->which_boundary = bound_ind;
				temp.push_back(he->v);
				now = he->next;

				while (now->flip != NULL)
				{
					now = now->flip->next;
				}
				he = now;
			}
			boundaries.push_back(temp);
			bound_ind++;
		}
		
	}

	printf("Number of boundaries found = %d\n", boundaries.size());
	int maxbound = 0;
	for (int i = 0; i < boundaries.size(); i++)
	{
		if (boundaries[i].size() > maxbound)
		{
			maxbound = boundaries[i].size();
			chart_boundary_ind = i;
		}
		current_boundary_size.push_back(boundaries[i].size());
		printf("\t num boundary = %d\n", boundaries[i].size());
	}
	printf("Main boundary considered = %d\n", chart_boundary_ind);
}

void HalfEdgeMesh::setup_first_col_lvl()
{
	if (col_lev.size() == 0)
	{
		CollapseLevel lvl;
		lvl.num_col = 0;

		vect3d X, Y, Z;
		for (int i = 0; i < totalFaces; i++)
		{
			if (!faceData[i].hole)
			{
				vector<vect2d> t;
				vect2d t1, t2, t3;
				t.push_back(t1);
				t.push_back(t2);
				t.push_back(t3);

				vect3d a = faceData[i].he->v->pos;
				vect3d b = faceData[i].he->next->v->pos;
				vect3d c = faceData[i].he->next->next->v->pos;

				X = b - a;
				X.normalize();
				Z = X % (c - a);
				Z.normalize();

				Y = Z % X;

				//// scale triangle by its initial scale (to get approximately the same size mesh as the input)
				double amb = 1;
				X *= amb;
				Y *= amb;

				// store
				t[0].set(0, 0);
				t[1].set((b - a) * X, 0);
				t[2].set((c - a) * X, (c - a) * Y);

				lvl.iso.push_back(t);
			}
		}
		printf("\tnum iso tris = %d\n", lvl.iso.size());
		//boundary
		for (int i = 0; i < boundaries.size(); i++)
		{
			vector<int> t;
			for (int j = 0; j < boundaries[i].size(); j++)
			{
				if (collapsed[boundaries[i][j]->index] == -1)
				{
					t.push_back(boundaries[i][j]->index);
				}
			}
			lvl.boundaries.push_back(t);
		}
		printf("\tnum boundaries = %d\n", lvl.boundaries.size());
		for (int i = 0; i < lvl.boundaries.size(); i++)
		{
			printf("\t\t %d\n", lvl.boundaries[i].size());
		}

		//faces
		for (int i = 0; i < totalFaces; i++)
		{
			if (!faceData[i].hole)
			{
				vect3i t;
				t[0] = faceData[i].he->v->index;
				t[1] = faceData[i].he->next->v->index;
				t[2] = faceData[i].he->next->next->v->index;

				lvl.faces.push_back(t);
			}
		}

		//one ring
		for (int i = 0; i < totalVerts; i++)
		{
			if (collapsed[i] == -1)
			{
				HE_HalfEdge *start;
				start = vertexData[i].he;
				vector<int> or;
				if (vertexData[i].boundary)
				{
					//find start
					while (start->next->flip != NULL)
					{
						start = start->next->flip;
					}
					start = start->next;
					or.push_back(start->v->index);

					//wrap to end
					while (start->prev->flip != NULL)
					{
						start = start->prev->flip;
						or.push_back(start->v->index);
					}
					or.push_back(start->next->v->index);
				}
				else
				{
					//wrap around to start
					HE_HalfEdge *now = start;
					do
					{
						now = now->flip;
						or.push_back(now->v->index);
						now = now->next->next;
					} while (now != start);

				}
				lvl.onerings.push_back(or);
			}
		}
		printf("\tonerings = %d\n", lvl.onerings.size());
		//push lvl
		col_lev.push_back(lvl);

	}
}

void HalfEdgeMesh::build_rings()
{
	

	inlist.resize( vertexData.size() );
	for(int i = 0; i < vertexData.size(); i++)
	{
		inlist[i] = false;
	}

	finlist.resize( faceData.size() );
	for(int i = 0; i < faceData.size(); i++)
	{
		finlist[i] = false;
	}

	vector<int> first;
	for(int i = 0; i < boundaries[this->chart_boundary_ind].size(); i++)
	{
		first.push_back( boundaries[this->chart_boundary_ind][i]->index );
		inlist[ boundaries[this->chart_boundary_ind][i]->index ] = true;
	}

	vector<vect3i> f_first;
	vector<int> fi_first;
	for(int  i = 0; i < first.size(); i++)
	{
		HE_HalfEdge *start = vertexData [ first[i] ].he, *curr;

		bool isbound = vertexData[ first[i] ].boundary;

		//if boundary vertex move start to left triangle
		if(isbound)
		{
			while(start->next->flip != NULL)
				start = start->next->flip;
		}

		curr = start;
		bool cont = true;
		do
		{
			//Looop through faces
			//outer ring
			if(!curr->f->inside_rings)
			{
				curr->f->inside_rings = true;
				vect3i t = curr->f->vi;
			/*	t[0] = curr->f->he->v->index;
				t[1] = curr->f->he->next->v->index;
				t[2] = curr->f->he->next->next->v->index;*/
				f_first.push_back(t);
				fi_first.push_back(curr->f->index);
			}

			if( curr->flip != NULL && !curr->flip->f->inside_rings)
			{
				curr->flip->f->inside_rings = true;
				vect3i t = curr->flip->f->vi;
				//t[0] = curr->flip->f->he->v->index;
				//t[1] = curr->flip->f->he->next->v->index;
				//t[2] = curr->flip->f->he->next->next->v->index;
				f_first.push_back(t);
				fi_first.push_back(curr->flip->f->index);
			}


			if(isbound)
			{
				if(curr->flip != NULL)
				{
					curr = curr->flip->next->next;
					cont = true;
				}
				else
					cont = false;
			}
			else
			{
				curr = curr->flip->next->next;
				cont = curr != start;
			}
		} while ( cont );
	}

	rings.push_back(first);
	rings_faces.push_back(f_first);
	rings_find.push_back(fi_first);

	int count = 0;
	int ringind = 0;
	do
	{
		count = 0;
		vector<int> new_ring;
		//loop through faces of prev level
		for(int i = 0; i < rings_faces[ringind].size(); i++)
		{
			//add vertices that are not in the current ring
			for(int ii = 0; ii < 3; ii++)
			{
				//bool addme = true;
				if(!inlist[rings_faces[ringind][i][ii]])
				{
					new_ring.push_back( rings_faces[ringind][i][ii] );
					inlist[rings_faces[ringind][i][ii]] = true;
					count++;
				}
			}
		}

		vector<vect3i> new_faces;
		vector<int> ni;

		for(int  i = 0; i < new_ring.size(); i++)
		{
			HE_HalfEdge *start = vertexData [ new_ring[i] ].he, *curr;

			bool isbound = vertexData[ new_ring[i] ].boundary;

			//if boundary vertex move start to left triangle
			if(isbound)
			{
				while(start->next->flip != NULL)
					start = start->next->flip;
			}

			curr = start;
			bool cont = true;
			do
			{
				//Looop through faces
				//outer ring
				if(!curr->f->inside_rings)
				{
					curr->f->inside_rings = true;
					vect3i t = curr->f->vi;
	//				t[0] = curr->f->he->v->index;
		//			t[1] = curr->f->he->next->v->index;
			//		t[2] = curr->f->he->next->next->v->index;
					new_faces.push_back(t);
					ni.push_back(curr->f->index);
				}

				if( curr->flip != NULL && !curr->flip->f->inside_rings)
				{
					curr->flip->f->inside_rings = true;
					vect3i t = curr->flip->f->vi;
					//t[0] = curr->flip->f->he->v->index;
					//t[1] = curr->flip->f->he->next->v->index;
					//t[2] = curr->flip->f->he->next->next->v->index;
					new_faces.push_back(t);
					ni.push_back(curr->flip->f->index);
				}


				if(isbound)
				{
					if(curr->flip != NULL)
					{
						curr = curr->flip->next->next;
						cont = true;
					}
					else
						cont = false;
				}
				else
				{
					curr = curr->flip->next->next;
					cont = curr != start;
				}
			} while ( cont );
		}
		

		ringind++;
		printf("count ring = %d\n", count);
		if(count > 0)
		{
		rings.push_back(new_ring);
		rings_faces.push_back(new_faces);
		rings_find.push_back(ni);
		//rings_isos.push_back(new_isos);
		}
	}
	while(count > 0);

}
void HalfEdgeMesh::init ( MeshData *mesh, MeshChart*charts )
{
	drawlvl = -1;
	collapse_iter = 0;

#ifdef TRACE_DEBUG
	trace_file = trace_file1;
#endif

	printf ( "%d verts, %d faces\n", mesh->positions.s, mesh->indices_pos.s );
	printf ( "building half edge structure\n" );
	//Timer::singleton ( )->resetTime ( 0 );
	map<pair<int, int>, HE_Edge*> edgeMap;
	int i, j;

	vertexData.resize( mesh->tex_coords.s );
	faceData.resize( mesh->indices_tex.s );
	heData.resize( mesh->indices_tex.s * 3   );
	// a closed triangle mesh always has an even number of faces
	edgeData.resize( ( mesh->indices_tex.s * 3 ) );
	collapsed.resize( mesh->tex_coords.s );

	totalFaces = mesh->indices_tex.s;
	totalVerts = mesh->tex_coords.s;
	int currentHE = 0;
	int currentEdge = 0;

	printf("CHARTS = %d\n", charts->chart_data.size());
	vector<int> chart_id;
	for(map<int, MeshChartData>::iterator it = charts->chart_data.begin(); it != charts->chart_data.end(); it++)
	{
		chart_id.push_back( it->first);
	}

	for (i = 0; i < mesh->tex_coords.s; i++)
	{
		vertexData[i].he = NULL;
		vertexData[i].pos = mesh->positions[i];
		vertexData[i].uvpos = mesh->tex_coords[i];
		vertexData[i].collapsable = true;
		vertexData[i].boundary = false;
		vertexData[i].which_boundary = -1;
		vertexData[i].index = i;
		vertexData[i].newindex = -1;
		vertexData[i].recompute_qef = true;

		for(j = 0; j < chart_id.size(); j++)
		{
			if(charts->find(&charts->verts[i])->index_mesh == chart_id[j])
				vertexData[i].chart_ind = j;
		}

		vertexData[i].qef.zero();
#ifdef LOAD_OBJ_JOETYPE
#else
		// WTF? - Sahil
		//		vertexData[i].pos[0] *= -1;
#endif
		collapsed[i] = -1;
	}

	for ( i = 0; i < mesh->indices_tex.s; i++ )
	{
		faceData [ i ].hole = false;
		faceData[i].inside_rings = false;
		faceData[i].removed = -1;
		faceData[i].index = i;
		faceData[i].vi[0] = mesh->indices_tex[i][0];
		faceData[i].vi[1] = mesh->indices_tex[i][1];
		faceData[i].vi[2] = mesh->indices_tex[i][2];

		HE_HalfEdge *prev = NULL, *start;
		for ( j = 0; j < 3; j++ )
		{
			int next = mesh->indices_tex[i] [( j + 1 ) % 3];
			int curr = mesh->indices_tex[i] [j];
			pair<int, int> key;
			HE_HalfEdge *he = &(heData [ currentHE++ ]);

			if ( j != 0 )
			{
				prev->next = he;
				he->prev = prev;
			}
			else
			{
				faceData [ i ].he = he;
				start = he;
			}
			prev = he;

			he->f = &(faceData [ i ]);
			he->v = &(vertexData [ next ]);
			he->flip = NULL;
			he->next = NULL;

#ifdef SCOTT_COLLAPSE
			vertexData [ curr ].he = he;
#else	
			vertexData [ next ].he = he;
#endif

			if ( curr < next )
			{
				key = pair<int, int> ( curr, next );
			}
			else
			{
				key = pair<int, int> ( next, curr );
			}

			map<pair<int, int>, HE_Edge*>::iterator it = edgeMap.find ( key );
			if ( it == edgeMap.end ( ) )
			{
				he->e = &(edgeData [ currentEdge++ ]);
				he->e->dirty = false;
				he->e->he = he;
				edgeMap [ key ] = he->e;
				he->flip = NULL;
			}
			else
			{
				he->e = it->second;
				he->flip = he->e->he;
				he->e->he->flip = he;
			}
		}
		prev->next = start;
		start->prev = prev;
	}

	currentVerts = totalVerts;
	//printf("totalVerts = %d\n", totalVerts);
	setup_boundaries();

	// allocate memory for the expansions
	//expandOps = new HE_ExpandOp [ totalVerts ];


//	printf ( "building qefs\n" );
//	//Timer::singleton ( )->resetTime ( 0 );
//	buildQEFs ( );
//	//printf ( "Took %f to build qefs\n", Timer::singleton ( )->resetTime ( 0 ) );
//
//	printf ( "inserting into queue\n" );
//	//Timer::singleton ( )->resetTime ( 0 );
//	
//	/*#ifdef TRACE_DEBUG
//	for ( map<pair<int, int>, HE_Edge *>::iterator it = edgeMap.begin ( ); it != edgeMap.end ( ); it++ )
//		fprintf(trace_file1, "edge %d %d\n", it->first.first, it->first.second);
//#endif
//	for ( map<pair<int, int>, HE_Edge *>::iterator it = edgeMap.begin ( ); it != edgeMap.end ( ); it++ )
//		errQueue.insert ( it->second, minimizeEdge ( it->second ) );
//*/
//	for ( map<pair<int, int>, HE_Edge *>::iterator it = edgeMap.begin ( ); it != edgeMap.end ( ); it++ )
//		errQueue.insert ( it->second, minimizeEdge ( it->second ) );
//	
//	//build_rings();
}

void HalfEdgeMesh::set_to_recalc_qef( HE_Vertex *v )
{
	v->recompute_qef = true;
	HE_HalfEdge *start = vertexData [ v->index ].he, *curr;

	bool isbound = vertexData[v->index].boundary;

	//if boundary vertex move start to left triangle
	if(isbound)
	{
		while(start->next->flip != NULL)
			start = start->next->flip;
	}

	curr = start;
	bool cont = true;
	do
	{
		curr->v->recompute_qef = true;
		curr->next->v->recompute_qef = true;
		curr->next->next->v->recompute_qef = true;

		if(isbound)
		{
			if(curr->flip != NULL)
			{
				curr = curr->flip->next->next;
				cont = true;
			}
			else
				cont = false;
		}
		else
		{
			curr = curr->flip->next->next;
			cont = curr != start;
		}
	} while ( cont );
}
void HalfEdgeMesh::buildQEFs ( void )
{
	int i;
	vect3d v1, v2, v3, norm;
	NUM_TYPE eq [ 4 ];

	for ( i = 0; i < totalVerts; i++ )
	{
		if(vertexData[i].recompute_qef)
		{
			vertexData [ i ].qef.zero ( );
			HE_HalfEdge *start = vertexData [ i ].he, *curr;

			bool isbound = vertexData[i].boundary;

			//if boundary vertex move start to left triangle
			if(isbound)
			{
				while(start->next->flip != NULL)
					start = start->next->flip;
			}

			curr = start;
			bool cont = true;
			do
			{
	//			if ( curr->f->hole == false )
				{
					v1 = curr->v->pos;
					v2 = curr->next->v->pos;
					v3 = curr->prev->v->pos;

#ifdef QEF_VERTEX
					eq[0] = v2[0];
					eq[1] = v2[1];
					eq[2] = v2[2];

					vertexData[ i ].qef.combineSelf(eq);
#else
	//				norm = (v2 - v1).cross(v3 - v1);
					norm = (v1 - v3).cross(v2 - v3);
	//				norm.normalize ( );
					eq [ 0 ] = norm[0];
					eq [ 1 ] = norm[1];
					eq [ 2 ] = norm[2];
					//eq [ 3 ] = -1*( norm.dot(v1) );
					eq [ 3 ] = -1*( norm.dot(v3) );
					vertexData [ i ].qef.combineSelf ( eq );

					// use centroid in cases where matrix is degenerate
					float w = 0.01 * norm.length ( );
					vect3d centroid = ( v1 + v2 + v3 ) / 3;
					eq [ 0 ] = w;
					eq [ 1 ] = 0;
					eq [ 2 ] = 0;
					eq [ 3 ] = -centroid[0] * w;
					vertexData [ i ].qef.combineSelf ( eq );
					eq [ 0 ] = 0;
					eq [ 1 ] = w;
					eq [ 2 ] = 0;
					eq [ 3 ] = -centroid[1] * w;
					vertexData [ i ].qef.combineSelf ( eq );
					eq [ 0 ] = 0;
					eq [ 1 ] = 0;
					eq [ 2 ] = w;
					eq [ 3 ] = -centroid[2] * w;
					vertexData [ i ].qef.combineSelf ( eq );
#endif
				}
			
				if(isbound)
				{
					if(curr->flip != NULL)
					{
						curr = curr->flip->next->next;
						cont = true;
					}
					else
						cont = false;
				}
				else
				{
					curr = curr->flip->next->next;
					cont = curr != start;
				}
			} while ( cont );
			vertexData[i].recompute_qef = false;

#ifdef QEF_VERTEX
			if(vertexData[i].boundary)
			{
					eq[0] = curr->next->next->v->pos[0];
					eq[1] = curr->next->next->v->pos[1];
					eq[2] = curr->next->next->v->pos[2];

					vertexData[ i ].qef.combineSelf(eq);
			}
#endif
		}
		//printf("qef %f %f %f %f %f %f %f %f %f %f\n", vertexData[i].qef.data[0], vertexData[i].qef.data[1], vertexData[i].qef.data[2], vertexData[i].qef.data[3], vertexData[i].qef.data[4], vertexData[i].qef.data[5], vertexData[i].qef.data[6], vertexData[i].qef.data[7], vertexData[i].qef.data[8], vertexData[i].qef.data[9]);
	}
}

void HalfEdgeMesh::print_buildQEFs ( HE_Vertex *v )
{
	int i;
	vect3d v1, v2, v3, norm;
	NUM_TYPE eq [ 4 ];

	//for ( i = 0; i < totalVerts; i++ )
	{
		//if(vertexData[i].recompute_qef)
		{
			//v->qef.zero ( );
			HE_HalfEdge *start = v->he, *curr;

			bool isbound = v->boundary;

			//if boundary vertex move start to left triangle
			if(isbound)
			{
				while(start->next->flip != NULL)
					start = start->next->flip;
			}

			curr = start;
			bool cont = true;
			do
			{
	//			if ( curr->f->hole == false )
				{
					v1 = curr->v->pos;
					v2 = curr->next->v->pos;
					v3 = curr->prev->v->pos;

					printf("{");
					printf("{%f, %f, %f},\n", v1[0], v1[1], v1[2]);
					printf("{%f, %f, %f},\n", v2[0], v2[1], v2[2]);
					printf("{%f, %f, %f} }\n", v3[0], v3[1], v3[2]);
				}
			
				if(isbound)
				{
					if(curr->flip != NULL)
					{
						curr = curr->flip->next->next;
						cont = true;
					}
					else
						cont = false;
				}
				else
				{
					curr = curr->flip->next->next;
					cont = curr != start;
				}
			} while ( cont );
			//vertexData[i].recompute_qef = false;
		}
		//printf("qef %f %f %f %f %f %f %f %f %f %f\n", vertexData[i].qef.data[0], vertexData[i].qef.data[1], vertexData[i].qef.data[2], vertexData[i].qef.data[3], vertexData[i].qef.data[4], vertexData[i].qef.data[5], vertexData[i].qef.data[6], vertexData[i].qef.data[7], vertexData[i].qef.data[8], vertexData[i].qef.data[9]);
	}
}


void HalfEdgeMesh::draw_3Dfull()
{
	float color_surf[4] = {.8,.8,.8, 1};
	int i, j;

	glEnable(GL_LIGHT0);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_CULL_FACE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_surf);

	//glColor3f(0,0,1);
	
	for(i = 0; i < this->faceData.size(); i++)
	{
		glBegin( GL_POLYGON );

		vect3f t;
		t= (( vertexData[faceData[i].vi[2]].pos - vertexData[faceData[i].vi[0]].pos ).cross( vertexData[faceData[i].vi[1]].pos - vertexData[faceData[i].vi[0]].pos ));
		t.normalize();
		glNormal3fv(
			t.v
			);

		glVertex3dv( vertexData[faceData[i].vi[0]].pos.v );
		glVertex3dv( vertexData[faceData[i].vi[2]].pos.v );
		glVertex3dv( vertexData[faceData[i].vi[1]].pos.v );

		glEnd();
	}

	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	//glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHT0);

	glLineWidth(15);
	glColor3f(0,0,1);
	glBegin(GL_LINES);
	for(int i = 0; i < this->boundaries[0].size(); i++)
	{
		glVertex3dv(this->boundaries[0][i]->pos.v);
		glVertex3dv(this->boundaries[0][(i+1)%this->boundaries[0].size()]->pos.v);
		/*double t = (g.smesh.tex_coords[g.sbo[i]] - g.smesh.tex_coords[g.sbo[(i+1)%g.sbo.size()]]).length();
			
			tot_length += t;

			max_length = t > max_length ? t : max_length;
			min_length = t < min_length ? t : min_length;*/

	}
	glEnd();

	glPointSize(BOUNDARY_THICKNESS);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(15);
	//glBegin(GL_POINTS);
	////glPointSize(20);
	//for(int i = 0; i < this->boundaries[0].size(); i++)
	//{
	//	glVertex3dv(this->boundaries[0][i]->pos.v);
	//}
	//glEnd();
}
void HalfEdgeMesh::draw_col_lvl(void)
{
	int i;
	glColor3f(0, 1, 1);
	glBegin(GL_LINES);
	for (i = 0; i < col_lev[drawlvl].faces.size(); i++)
	{
		int i0, i1, i2;
		/*i0 = indchange[col_lev[drawlvl].faces[i][0]];
		i1 = indchange[col_lev[drawlvl].faces[i][1]];
		i2 = indchange[col_lev[drawlvl].faces[i][2]];*/

		i0 = col_lev[drawlvl].faces[i][0];
		i1 = col_lev[drawlvl].faces[i][1];
		i2 = col_lev[drawlvl].faces[i][2];

		glVertex2dv((vertexData[i0].uvpos.v));
		glVertex2dv((vertexData[i1].uvpos.v));

		glVertex2dv((vertexData[i1].uvpos.v));
		glVertex2dv((vertexData[i2].uvpos.v));

		glVertex2dv((vertexData[i0].uvpos.v));
		glVertex2dv((vertexData[i2].uvpos.v));
	}
	glEnd();

	

	

	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	for(int i = 0; i < col_lev[drawlvl].boundaries.size(); i++)
	{
		for(int j = 0; j < col_lev[drawlvl].boundaries[i].size(); j++)
		{
			glVertex2dv( vertexData[col_lev[drawlvl].boundaries[i][j]].uvpos.v);
			glVertex2dv( vertexData[(col_lev[drawlvl].boundaries[i][(j+1)% col_lev[drawlvl].boundaries[i].size()]) ].uvpos.v);
		}
	}
	glEnd();
	
}

double find_a(vect2d v, vect2d w)
{
	vect3d vv;
	vv[0] = v[0]; vv[1] = v[1]; vv[2] = 0;

	vect3d ww;
	ww[0] = w[0]; ww[1] = w[1]; ww[2] = 0;

	return vv.cross(ww)[2] >= 0 ?  acos( v.dot(w) / (v.length()*w.length() ) ) : 2*PI - acos( v.dot(w) / (v.length()*w.length() ) );
}

void HalfEdgeMesh::draw(void)
{
	unsigned int i;

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_LINE_SMOOTH);
	for (i = 0; i < totalFaces; i++ )
	{
		if ( !faceData [ i ].hole )
		{
			HE_HalfEdge *he = faceData [ i ].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;

#ifdef DEBUG_COLORING
			if(this->err.size() > 0)
			{
				glColor3f(err[i],0,0);
			}
			else
#endif
				glColor3f(.95,.95,.95);
			glBegin(GL_POLYGON);
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));
			glEnd();


			glLineWidth(.01);
			glColor3f(.4, .4, .4);
			glBegin ( GL_LINES );
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));

			glVertex2dv((v3.v));
			glVertex2dv((v1.v));
			glEnd ( );
		}
	}

	for(i = 0; i < error_ind.size(); i++)
	{
		HE_HalfEdge *he = faceData [ error_ind[i] ].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;
		glColor3f(1,1,0);
			glBegin(GL_POLYGON);
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));
			glEnd();
	}

	

	for ( i = 0; i < totalFaces; i++ )
	{
		if ( !faceData [ i ].hole )
		{
			HE_HalfEdge *he = faceData [ i ].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;

			/*glColor3f(.9,.9,.9);
			glBegin(GL_POLYGON);
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));
			glEnd();*/


			glLineWidth(.5);
			/*if(this->err.size() > 0)
			{
				glColor3f(err[i],0,0);
			}
			else*/
				glColor3f(.4, .4, .4);
			glBegin ( GL_LINES );
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));

			glVertex2dv((v3.v));
			glVertex2dv((v1.v));
			glEnd ( );
		}
	}
	

	
	glLineWidth(BOUNDARY_THICKNESS);
	
	glBegin(GL_LINES);
	//

    glColor3f(0,0,1);
	
	for ( i = 0; i < totalFaces; i++ )
	{
		if ( !faceData [ i ].hole )
		{
			HE_HalfEdge *he = faceData [ i ].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;

			/*vect2d norm = MathVec3Cross(v2 - v1, v3 - v1);
			norm.normalize ( );
			GL_NORMAL ( &(norm.x) );*/
			if(he->next->flip == NULL)
			{
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			}
			if(he->next->next->flip == NULL)
			{
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));
			}
			if(he->next->next->next->flip == NULL)
			{
			glVertex2dv((v3.v));
			glVertex2dv((v1.v));
			}
		}

	}
	glEnd();

	glPointSize(BOUNDARY_THICKNESS);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);
	for ( i = 0; i < totalFaces; i++ )
	{
		if ( !faceData [ i ].hole )
		{
			HE_HalfEdge *he = faceData [ i ].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;

			/*vect2d norm = MathVec3Cross(v2 - v1, v3 - v1);
			norm.normalize ( );
			GL_NORMAL ( &(norm.x) );*/
			if(he->next->flip == NULL)
			{
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			}
			if(he->next->next->flip == NULL)
			{
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));
			}
			if(he->next->next->next->flip == NULL)
			{
			glVertex2dv((v3.v));
			glVertex2dv((v1.v));
			}
		}

	}
	glEnd();

	/*glBegin(GL_POINTS);

	glColor3f(0,1,0);
	glVertex2dv(zero_param.v);

	glColor3f(1,1,0);
	glVertex2dv(pie2_param.v);

	glEnd();*/
#ifdef DEBUG_COLORING
	if(error_view.size() > 0)
	{
	glBegin(GL_LINES);
	glColor3f(0,1,0);
	glVertex3f(0,0,0);
	glVertex3f(PI * 2 * 10,0,0);

	glVertex3f(0,10,0);
	glVertex3f(PI * 2 * 10,10,0);

	glVertex3f(PI * 2 * 10,0,0);
	glVertex3f(PI * 2 * 10,10,0);

	glVertex3f(PI  * 10,0,0);
	glVertex3f(PI  * 10,10,0);

	glVertex3f(PI *.5 * 10,0,0);
	glVertex3f(PI *.5 * 10,10,0);

	glVertex3f(PI *1.5 * 10,0,0);
	glVertex3f(PI *1.5 * 10,10,0);

	glVertex3f(0,0,0);
	glVertex3f(0,10,0);
	glEnd();
	}
	glBegin(GL_POINTS);
	for(int i = 0; i < this->error_view.size(); i++)
	{
		glColor3f(0,0,1);
		glVertex2dv(error_view[i].v);
	}
	glEnd();


	vect2d s;
	s[0] = 1; s[1] = 0;
	
	if(bins.size() > 0)
	{
	glBegin(GL_LINES);
	glColor3f(0,0,1);

	for(i = 0; i < this->boundaries[0].size(); i++)
	{
		vect2d t;
		t[0] = bins[i]*10;//10*find_a(s,boundaries[0][i]->uvpos);

		t[1] = 0;

		glColor3f(10*bin_errors[i],0,0);

		//glVertex2dv(t.v);
		glVertex2dv(boundaries[0][i]->uvpos.v);
		glVertex2dv(boundaries[0][(i+1)%boundaries[0].size()]->uvpos.v);
	}
	glEnd();

	glBegin(GL_LINES);
	
	glColor3f(0,1,1);

	for(i = 0; i < pi.size(); i++)
	{
		glVertex2dv( boundaries[0][pi[i].start]->uvpos.v );
		glVertex2dv( boundaries[0][pi[i].end]->uvpos.v );
	}
	glEnd();

	glBegin(GL_POINTS);
	glColor3f(0,1,0);
	for(i = 0; i < peaks.size(); i++)
	{
		glVertex2dv(boundaries[0][peaks[i]]->uvpos.v);
	}
	glEnd();
	}
#endif
	
}
	



vect2d j_func_grad(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2, int ind)
{
	vect2d ret;
	//printf("grad");
	if(ind == 0)
	{
		ret[0] = (2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*(((2.0)*(u0))+(-(u1))+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u1)+(-(u2))))+(((v0)+(-(v2)))*((v1)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+(((L)*(L)*(L))*((-((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*((v0)+(-(v2)))*((v1)+(-(v2))))))+(((u1)*(u1))*((v0)+(-(v2))))+(-(u0)*((u1)+(-(u2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u1)*(u2)*((-(v1))+(v2))))*(xx)*((y)*(y)))+(-((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((-((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((xx)*(xx))*((y)*(y)))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((y)*(y)*(y)*(y)))))));
		ret[1] = (1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((2.0)*((u1)+(-(u2)))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((u1)+(-(u2)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((2.0)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((v0)+(-(v2))))+((L)*(((-2.0)*(v0))+(v1)+(v2))*(xx))+(((v0)+(-(v1)))*(((xx)*(xx))+((y)*(y))))))));
		printf("NumberForm[{%.10g, %.10g},10]\n",
			(2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*(((2.0)*(u0))+(-(u1))+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u1)+(-(u2))))+(((v0)+(-(v2)))*((v1)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+(((L)*(L)*(L))*((-((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*((v0)+(-(v2)))*((v1)+(-(v2))))))+(((u1)*(u1))*((v0)+(-(v2))))+(-(u0)*((u1)+(-(u2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u1)*(u2)*((-(v1))+(v2))))*(xx)*((y)*(y)))+(-((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((-((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((xx)*(xx))*((y)*(y)))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((y)*(y)*(y)*(y)))))))
			,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((2.0)*((u1)+(-(u2)))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((u1)+(-(u2)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((2.0)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((v0)+(-(v2))))+((L)*(((-2.0)*(v0))+(v1)+(v2))*(xx))+(((v0)+(-(v1)))*(((xx)*(xx))+((y)*(y)))))))));
	}
	else if(ind == 1)
	{
		ret[0] = (1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+((-2.0)*((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+((2.0)*((L)*(L)*(L))*((((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*(((v0)+(-(v2)))*((v0)+(-(v2)))))))+((u1)*(u2)*((v0)+(-(v2))))+(((u0)*(u0))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u0)*(((u1)*((-(v0))+(v2)))+((u2)*(((-3.0)*(v0))+((2.0)*(v1))+(v2))))))*(xx)*((y)*(y)))+((2.0)*((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+((-2.0)*((L)*(L))*((v0)+(-(v1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))));
		ret[1] = (2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*((v0)+(-(v2)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((u0)+(-(u2)))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L)*(L))*(((-2.0)*((u0)*(u0)*(u0)))+((2.0)*((u0)*(u0))*((u1)+((2.0)*(u2))))+((u1)*(((2.0)*((u2)*(u2)))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+(-(u0)*(((2.0)*(u2)*(((2.0)*(u1))+(u2)))+(((v0)+(-(v2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+((u2)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+(((v0)+(-(v1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((u0)+(-(u1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))));
		printf("NumberForm[{%.10g, %.10g},10]\n",
		(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+((-2.0)*((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+((2.0)*((L)*(L)*(L))*((((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*(((v0)+(-(v2)))*((v0)+(-(v2)))))))+((u1)*(u2)*((v0)+(-(v2))))+(((u0)*(u0))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u0)*(((u1)*((-(v0))+(v2)))+((u2)*(((-3.0)*(v0))+((2.0)*(v1))+(v2))))))*(xx)*((y)*(y)))+((2.0)*((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+((-2.0)*((L)*(L))*((v0)+(-(v1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))))
		,(2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*((v0)+(-(v2)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((u0)+(-(u2)))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L)*(L))*(((-2.0)*((u0)*(u0)*(u0)))+((2.0)*((u0)*(u0))*((u1)+((2.0)*(u2))))+((u1)*(((2.0)*((u2)*(u2)))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+(-(u0)*(((2.0)*(u2)*(((2.0)*(u1))+(u2)))+(((v0)+(-(v2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+((u2)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+(((v0)+(-(v1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((u0)+(-(u1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y)))))));
	}
	else
	{
		ret[0] = (1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((L)*((u0)+(-(u2))))+(((-(u0))+(u1))*(xx)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y)))))+((2.0)*((-(v0))+(v1))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((-(v0))+(v1))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))));
		ret[1] = (2.0/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*((((v0)+(-(v1)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(-((L)*(L)*(L))*((u0)+(-(u2)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L))*(((2.0)*((u0)*(u0)*(u0)))+((-2.0)*((u1)*(u1))*(u2))+((-2.0)*((u0)*(u0))*(((2.0)*(u1))+(u2)))+(-(u2)*(((v0)+(-(v1)))*((v0)+(-(v1)))))+((u0)*(((2.0)*(u1)*((u1)+((2.0)*(u2))))+(((v0)+(-(v1)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+(-(u1)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+((L)*((((v0)+(-(v2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((xx)*(xx))*((y)*(y)))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((y)*(y)*(y)*(y)))))));
		printf("NumberForm[{%.10g, %.10g},10]\n",(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((L)*((u0)+(-(u2))))+(((-(u0))+(u1))*(xx)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y)))))+((2.0)*((-(v0))+(v1))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((-(v0))+(v1))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))))
	, (2.0/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*((((v0)+(-(v1)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(-((L)*(L)*(L))*((u0)+(-(u2)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L))*(((2.0)*((u0)*(u0)*(u0)))+((-2.0)*((u1)*(u1))*(u2))+((-2.0)*((u0)*(u0))*(((2.0)*(u1))+(u2)))+(-(u2)*(((v0)+(-(v1)))*((v0)+(-(v1)))))+((u0)*(((2.0)*(u1)*((u1)+((2.0)*(u2))))+(((v0)+(-(v1)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+(-(u1)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+((L)*((((v0)+(-(v2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((xx)*(xx))*((y)*(y)))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((y)*(y)*(y)*(y))))))));
	}
	return ret;
	/*printf("grad = {%.10g, %.10g, %.10g, %.10g,%.10g, %.10g}\n", (2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*(((2.0)*(u0))+(-(u1))+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u1)+(-(u2))))+(((v0)+(-(v2)))*((v1)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+(((L)*(L)*(L))*((-((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*((v0)+(-(v2)))*((v1)+(-(v2))))))+(((u1)*(u1))*((v0)+(-(v2))))+(-(u0)*((u1)+(-(u2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u1)*(u2)*((-(v1))+(v2))))*(xx)*((y)*(y)))+(-((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((-((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((xx)*(xx))*((y)*(y)))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((y)*(y)*(y)*(y)))))))
	,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((2.0)*((u1)+(-(u2)))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((u1)+(-(u2)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((2.0)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((v0)+(-(v2))))+((L)*(((-2.0)*(v0))+(v1)+(v2))*(xx))+(((v0)+(-(v1)))*(((xx)*(xx))+((y)*(y))))))))
	,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+((-2.0)*((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+((2.0)*((L)*(L)*(L))*((((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*(((v0)+(-(v2)))*((v0)+(-(v2)))))))+((u1)*(u2)*((v0)+(-(v2))))+(((u0)*(u0))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u0)*(((u1)*((-(v0))+(v2)))+((u2)*(((-3.0)*(v0))+((2.0)*(v1))+(v2))))))*(xx)*((y)*(y)))+((2.0)*((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+((-2.0)*((L)*(L))*((v0)+(-(v1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))))
	,(2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*((v0)+(-(v2)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((u0)+(-(u2)))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L)*(L))*(((-2.0)*((u0)*(u0)*(u0)))+((2.0)*((u0)*(u0))*((u1)+((2.0)*(u2))))+((u1)*(((2.0)*((u2)*(u2)))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+(-(u0)*(((2.0)*(u2)*(((2.0)*(u1))+(u2)))+(((v0)+(-(v2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+((u2)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+(((v0)+(-(v1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((u0)+(-(u1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))))
	,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((L)*((u0)+(-(u2))))+(((-(u0))+(u1))*(xx)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y)))))+((2.0)*((-(v0))+(v1))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((-(v0))+(v1))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))))
	, (2.0/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*((((v0)+(-(v1)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(-((L)*(L)*(L))*((u0)+(-(u2)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L))*(((2.0)*((u0)*(u0)*(u0)))+((-2.0)*((u1)*(u1))*(u2))+((-2.0)*((u0)*(u0))*(((2.0)*(u1))+(u2)))+(-(u2)*(((v0)+(-(v1)))*((v0)+(-(v1)))))+((u0)*(((2.0)*(u1)*((u1)+((2.0)*(u2))))+(((v0)+(-(v1)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+(-(u1)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+((L)*((((v0)+(-(v2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((xx)*(xx))*((y)*(y)))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((y)*(y)*(y)*(y)))))))
	);*/
}
void HalfEdgeMesh::print_debug_info(int ind)
{
	printf("vert %d\n", ind);
	HE_HalfEdge *start =vertexData [ ind ].he, *curr;
		bool isbound = vertexData[ ind ].boundary;

		//if boundary vertex move start to left triangle
		if(isbound)
		{
			while(start->next->flip != NULL)
				start = start->next->flip;
		}

		for(int j = 0; j < totalFaces; j++)
		{
			faceData[j].inside_rings = false;
		}

		curr = start;
		vect2d grad;
		grad[0] = 0; grad[1] = 0;
		bool cont = true;
		do
		{
			//Looop through faces
			//outer ring
			if(!curr->f->inside_rings)
			{
				curr->f->inside_rings = true;
				vect3i t;
				t[0] = curr->f->he->v->index;
				t[1] = curr->f->he->next->v->index;
				t[2] = curr->f->he->next->next->v->index;
				int vert_ind;
				if(t[0] == ind)
					vert_ind = 0;
				if(t[1] == ind)
					vert_ind = 1;
				if(t[2] == ind)
					vert_ind = 2;
				/*printf("%.10g, %.10g, %.10g\n", g.iso_tris[593][1][0],
						g.iso_tris[593][2][0],
						g.iso_tris[593][2][1]);*/
				int faceind = curr->f->index;
			/*	printf("NumberForm[GradJ[{{%.10g, %.10g}\n, {%.10g, %.10g}\n, {%.10g, %.10g}},\n %.10g, %.10g, %.10g],10] \n ",
					vertexData[t[0]].uvpos[0],vertexData[t[0]].uvpos[1],
					vertexData[t[1]].uvpos[0],vertexData[t[1]].uvpos[1],
					vertexData[t[2]].uvpos[0],vertexData[t[2]].uvpos[1],
					iso_tris[faceind][1][0],iso_tris[faceind][2][0],iso_tris[faceind][2][1]
					
					);*/

				grad+= j_func_grad(iso_tris[faceind][1][0],iso_tris[faceind][2][0],iso_tris[faceind][2][1],
					vertexData[t[0]].uvpos[0],vertexData[t[0]].uvpos[1],
					vertexData[t[1]].uvpos[0],vertexData[t[1]].uvpos[1],
					vertexData[t[2]].uvpos[0],vertexData[t[2]].uvpos[1],
					vert_ind
					);

			}

			if( curr->flip != NULL && !curr->flip->f->inside_rings)
			{
				curr->flip->f->inside_rings = true;
				vect3i t;
				t[0] = curr->flip->f->he->v->index;
				t[1] = curr->flip->f->he->next->v->index;
				t[2] = curr->flip->f->he->next->next->v->index;
				int vert_ind;
				if(t[0] == ind)
					vert_ind = 0;
				if(t[1] == ind)
					vert_ind = 1;
				if(t[2] == ind)
					vert_ind = 2;
				int faceind = curr->flip->f->index;


				iso_tris[faceind][1][0];

				/*printf("NumberForm[GradJ[{{%.10g, %.10g}\n, {%.10g, %.10g}\n, {%.10g, %.10g}},\n %.10g, %.10g, %.10g],10] \n",
					vertexData[t[0]].uvpos[0],vertexData[t[0]].uvpos[1],
					vertexData[t[1]].uvpos[0],vertexData[t[1]].uvpos[1],
					vertexData[t[2]].uvpos[0],vertexData[t[2]].uvpos[1],
					iso_tris[faceind][1][0],iso_tris[faceind][2][0],iso_tris[faceind][2][1]
					
					);*/

					grad+= j_func_grad(iso_tris[faceind][1][0],iso_tris[faceind][2][0],iso_tris[faceind][2][1],
					vertexData[t[0]].uvpos[0],vertexData[t[0]].uvpos[1],
					vertexData[t[1]].uvpos[0],vertexData[t[1]].uvpos[1],
					vertexData[t[2]].uvpos[0],vertexData[t[2]].uvpos[1],
					vert_ind
					);
			}

			if(isbound)
			{
				if(curr->flip != NULL)
				{
					curr = curr->flip->next->next;
					cont = true;
				}
				else
					cont = false;
			}
			else
			{
				curr = curr->flip->next->next;
				cont = curr != start;
			}
		} while ( cont );

		//printf("NumberForm[{%.10g, %.10g},10];\n", grad[0], grad[1]);

}

void HalfEdgeMesh::draw_full(void)
{
	int i;

	glColor3f(0, 1, 1);
	glBegin(GL_LINES);
	int count = 0;
	for (i = 0; i < totalFaces; i++)
	{
		if (!faceData[i].hole)
		{
			HE_HalfEdge *he = faceData[i].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;

			/*vect2d norm = MathVec3Cross(v2 - v1, v3 - v1);
			norm.normalize ( );
			GL_NORMAL ( &(norm.x) );*/
			/*if(he != he->prev->prev->prev)
			glColor3f(1, 0, 0);
			else
			glColor3f(0, 1, 1);*/

			glVertex2dv((v1.v));
			glVertex2dv((v2.v));

			glVertex2dv((v2.v));
			glVertex2dv((v3.v));

			glVertex2dv((v3.v));
			glVertex2dv((v1.v));
			//count++;
		}
	}
	glEnd();



	


	//glColor3f(0, 1, 1);
	//glBegin ( GL_LINES );
	//for ( i = 0; i < totalFaces; i++ )
	//{
	//	if ( !faceData [ i ].hole )
	//	{
	//		HE_HalfEdge *he = faceData [ i ].he;
	//		vect2d v1 = he->v->uvpos;
	//		vect2d v2 = he->next->v->uvpos;
	//		vect2d v3 = he->prev->v->uvpos;

	//		/*vect2d norm = MathVec3Cross(v2 - v1, v3 - v1);
	//		norm.normalize ( );
	//		GL_NORMAL ( &(norm.x) );*/
	//		if((he->flip == NULL  || he->flip->flip != he ) && (he->next->flip == NULL || he->next->flip->flip != he->next) && 
	//			( he->next->next->flip == NULL ||  he->next->next->flip->flip != he->next->next))
	//		{
	//			glColor3f(1, 0, 0);
	//			glVertex2dv((v1.v));
	//			glVertex2dv((v2.v));
	//		
	//			glVertex2dv((v2.v));
	//			glVertex2dv((v3.v));

	//			glVertex2dv((v3.v));
	//			glVertex2dv((v1.v));
	//		}
	//		else
	//			
	//		{
	//		}

	//		
	//	}
	//}
	//glEnd ( );

	//glColor3f(1,0,0);
	//glPointSize(10);

	//glBegin(GL_POINTS);
	//
	////show collapseable
	////for(int i = 0; i < totalVerts; i++)
	////{
	////	//if(!vertexData[i].collapsable && vertexData[i].when_col == -1)
	////	
	////	if (collapsed[i] == -1)
	////		//if( vertexData[i].boundary )
	////			if( !vertexData[i].collapsable )
	////		glVertex2dv( vertexData[i].uvpos.v );
	////}

	//for (int i = 0; i < totalVerts; i++)
	//{
	//	if (collapsed[i] == -1)
	//	{
	//		if (vertexData[i].which_boundary > -1)
	//		{
	//			if (vertexData[i].which_boundary == 0)
	//			{
	//				//count++;
	//				glColor3f(0, 1, 0);
	//			}
	//			else if (vertexData[i].which_boundary == 1)
	//				glColor3f(1,0, 0);
	//			else 
	//				glColor3f(0, 0, 1);
	//			glVertex2dv(vertexData[i].uvpos.v);
	//		}
	//		else
	//		{
	//			glColor3f(0, 1, 1);
	//			glVertex2dv(vertexData[i].uvpos.v);
	//		}
	//	}
	//}
	////printf("...%d\n", count);
	//glEnd();

	

	glColor3f(1, 0, 1);
	glBegin(GL_LINES);
	for (i = 0; i < totalFaces; i++)
	{
		if (!faceData[i].hole)
		{

			float u0 = faceData[i].he->v->uvpos[0];
			float v0 = faceData[i].he->v->uvpos[1];

			float u1 = faceData[i].he->next->v->uvpos[0];
			float v1 = faceData[i].he->next->v->uvpos[1];

			float u2 = faceData[i].he->next->next->v->uvpos[0];
			float v2 = faceData[i].he->next->next->v->uvpos[1];

			float A = (u2*(-v0 + v1) + u1*(v0 - v2) + u0*(-v1 + v2));
			//if (area >= 0)
			if (A < 0)
			{
				printf("\n flipped = %d\n", i);
			/*	glVertex2d(u0, v0);
				glVertex2d(u1, v1);

				glVertex2d(u1, v1);
				glVertex2d(u2, v2);

				glVertex2d(u2, v2);
				glVertex2d(u0, v0);

				 u0 = faceData[i].he->v->pos[0];
			 v0 = faceData[i].he->v->pos[1];

			 u1 = faceData[i].he->next->v->pos[0];
			 v1 = faceData[i].he->next->v->pos[1];

			 u2 = faceData[i].he->next->next->v->pos[0];
			 v2 = faceData[i].he->next->next->v->pos[1];

			 printf("%.10g, %.10g, %.10g\n %.10g, %.10g, %.10g\n %.10g, %.10g, %.10g", 
				 u0,v0,faceData[i].he->v->pos[2],
				 u0,v0,faceData[i].he->next->v->pos[2],
				 u0,v0,faceData[i].he->next->next->v->pos[2]);

			 printf("%d, %d, %d\n", faceData[i].he->v->index, faceData[i].he->next->v->index, faceData[i].he->next->next->v->index);
			 printf("{{%.10g, %.10g},{%.10g, %.10g},{%.10g, %.10g}}", vertexData[faceData[i].he->v->index].uvpos[0],
				 vertexData[faceData[i].he->v->index].uvpos[1], 
				 vertexData[faceData[i].he->next->v->index].uvpos[0], 
				 vertexData[faceData[i].he->next->v->index].uvpos[1], 
				 vertexData[faceData[i].he->next->next->v->index].uvpos[0], 
				 vertexData[faceData[i].he->next->next->v->index].uvpos[1]);
			 print_debug_info(faceData[i].he->v->index);
			 print_debug_info(faceData[i].he->next->v->index);
			 print_debug_info(faceData[i].he->next->next->v->index);

			 int s;
			 scanf("%d",s);*/
			}
		}
	}
	


	glBegin(GL_LINES);
	//

    glColor3f(1,0,0);
	for ( i = 0; i < totalFaces; i++ )
	{
		if ( !faceData [ i ].hole )
		{
			HE_HalfEdge *he = faceData [ i ].he;
			vect2d v1 = he->v->uvpos;
			vect2d v2 = he->next->v->uvpos;
			vect2d v3 = he->prev->v->uvpos;

			/*vect2d norm = MathVec3Cross(v2 - v1, v3 - v1);
			norm.normalize ( );
			GL_NORMAL ( &(norm.x) );*/
			if(he->next->flip == NULL)
			{
			glVertex2dv((v1.v));
			glVertex2dv((v2.v));
			}
			if(he->next->next->flip == NULL)
			{
			glVertex2dv((v2.v));
			glVertex2dv((v3.v));
			}
			if(he->next->next->next->flip == NULL)
			{
			glVertex2dv((v3.v));
			glVertex2dv((v1.v));
			}
		}

		/*if(debugedge.size() > 0)
		{
			glVertex2dv(debugedge[0].v);
			glVertex2dv(debugedge[1].v);
		}*/
	}
	glEnd();



	//glColor3f(1, 0, 0);
	//glBegin(GL_POINTS);
	///*for (int i = 0; i < totalVerts; i++)
	//{
	//	if (collapsed[i] == -1)
	//	{
	//		if (!this->vertexData[i].collapsable)
	//			glColor3f(1, 0, 0);
	//		else
	//			glColor3f(0, 1, 1);
	//		glVertex2dv(vertexData[i].uvpos.v);
	//	}
	//}*/
	//glVertex2dv(vertexData[482].uvpos.v);
	//glEnd();
	glPointSize(1);
	glBegin(GL_POINTS);
	for(int i = 0; i < rings.size(); i++)
	{
		glColor3f ( (i%8) >> 2, (((i%8)>>1)&1), ((i%8)&1) );
		if(i%3 == 0)
			glColor3f(1,0,0);
		else if(i%3==1)
			glColor3f(0,1,0);
		else
			glColor3f(0,0,1);
		for(int j = 0; j < rings[i].size(); j++)
		{
			glVertex2dv(vertexData[rings[i][j]].uvpos.v);
		}
	}
	glEnd();

	glPointSize(5);
	glColor3f(0,1,0);
	glBegin(GL_POINTS);
	/*glVertex2dv( vertexData[137078].uvpos.v);
	glVertex2dv( vertexData[137119].uvpos.v);
	glVertex2dv( vertexData[137090].uvpos.v);*/

	glEnd();

}

bool HalfEdgeMesh::isSafeCollapse ( HE_Edge *e )
{

	if(e->he->v->boundary && e->he->prev->v->boundary)
	{//if two boundary vertices

		if(e->he->flip != NULL)
			return false;//boundary to boundary not flipped

		if(e->he->next->flip == NULL && e->he->prev->flip == NULL)
			return false;//Would create hangline line

		//This was wing test... was a bit to restrictive
		//if (e->he->next->v->boundary)
		//	return false; //Would expand to a single line... will collapse more appropriately with later collapse

		if(e->he->next->flip == NULL || e->he->next->next->flip == NULL) //Checks if more than one edge in triangle is a boundary
			return false; //Would expand to single line
	}


	HE_HalfEdge *he = e->he;
	if (e->he->flip == NULL)
	{
		if(is_0_areaF(e->he))
		{
			printf("wtf this collapse would cause 0 area\n");
			return false;
		}

		//need to check if collapse makes 0 iso...
		vect3d tind = e->he->prev->v->pos;
		HE_HalfEdge *the = e->he;
		while(the->next->flip != NULL)
			the = the->next->flip;

		vect3d tind2 = the->next->v->pos;
		if((tind - tind2).length2() < .000000001)
			return false;
	}
	else if( isboundary(he) )
	{

		if(is_0_area(he))
		{
			printf("wtf this collapse would cause 0 area\n");
			return false;
		}
		//return false;

		/*if (is_3_val(he))
			return false;

		if (is_3_val(he->flip))
			return false;*/

		if( !(( he->next->flip != NULL || he->prev->flip != NULL ) && ( he->flip->next->flip != NULL || he->flip->prev->flip != NULL ) ) )
		{
			//Wing Neighbors don't exist... would create haninging edge
			return false;
		}
	}
	else
	{
		if(is_0_area(he))
		{
			printf("wtf this collapse would cause 0 area\n");
			return false;
		}
		//return false;


		/*if (is_3_val(he))
			return false;

		if (is_3_val(he->flip))
			return false;*/
	}

	//Check valence 3 neighbors shared with both vertices to be collapsed
	set<HE_Vertex *> set1, set2, result;
	bool flip = false;

	HE_HalfEdge *start = e->he;
	HE_HalfEdge *curr = start;

	set1.insert(start->next->next->v);
	set1.insert(start->v);

	do
	{
		if (curr->flip == NULL)
		{
			flip = true;
			break;
		}
		curr = curr->flip->next;
		set1.insert(curr->v);
	} while (curr != start);

	if (flip)
	{
		curr = start->next->next;
		set1.insert(curr->prev->v);

		while (curr->flip != NULL)
		{
			curr = curr->flip;
			curr = curr->next->next;
			set1.insert(curr->prev->v);
		}
	}

	flip = false;
	start = e->he->next;
	set2.insert(start->next->next->v);
	set2.insert(start->v);
	curr = start;

	do
	{
		if (curr->flip == NULL)
		{
			flip = true;
			break;
		}
		curr = curr->flip->next;
		set2.insert(curr->v);
	} while (curr != start);

	if (flip)
	{
		curr = start->next->next;
		set2.insert(curr->prev->v);

		while (curr->flip != NULL)
		{
			curr = curr->flip;
			curr = curr->next->next;
			set2.insert(curr->prev->v);
		}
	}

	flip = false;

	
	insert_iterator<set<HE_Vertex *, less<HE_Vertex *> > > resultIt ( result, result.begin ( ) );
	set_intersection ( set1.begin ( ), set1.end ( ), set2.begin ( ), set2.end ( ), resultIt );

	if(e->he->flip != NULL)
		return result.size ( ) == 4;
	else
		return result.size ( ) == 3;
}

#ifdef SCOTT_COLLAPSE
void HalfEdgeMesh::set_one_ring_uncollapsable(HE_Vertex *v)
{
	v->collapsable = false;
	HE_HalfEdge *now = v->he;

	bool goback = false;
	now->v->collapsable = false;

	int count = 0;
	do
	{
		now->next->v->collapsable=false;
		now = now->next->next->flip;
		if(now==NULL)
		{
			goback = true;
			break;
		}
		count++;
	}while( now != v->he && count < 100);
	
	if(goback)
	{
		now = v->he;
		do
		{
			now->v->collapsable = false;
			now = now->flip;
			if(now==NULL)
			{
				break;
			}
			now = now->next;
			count++;
		}while( now != v->he && count < 100 );
	}
	if(count == 100)
		printf("WTF\n");
	
}

HE_HalfEdge* HalfEdgeMesh::find_collapsable_edge_boundary(HE_Vertex v)
{
	HE_HalfEdge *now = v.he;

	now = now->next->next;
	if(now->flip == NULL)
		return now;

	int count = 0;
	do
	{
		if(now->flip == NULL)
			return now;

		now = now->flip->next->next;
		count++;
	}while(/*now->flip != NULL*/ count < 100);

//	bool goback = false;
//	if (now == NULL)
//	{
//		goback = true;
////		return v.he;
//	}
//
//	int count = 0;
//	if(goback)
//	{
//	}
//	else
//	{
//
//	do
//	{
//		if(now->next->flip == NULL)
//		{
//			/*if(v.boundary && now->next->v->boundary)
//			{
//				printf("good\n");
//			}*/
//			if(&v != )
//			return now->next;
//		}
//		now = now->next->flip;
//		count++;
//	}while(  count < 100 );
//	}


	printf("shit find boundaryedge = %d\n", count);
	return v.he;
}


void HalfEdgeMesh::collapse_one_iteration ( Array<vect3i> indices_pos )
{
	/*for (int i = 0; i < totalVerts; i++)
	{
		if (collapsed[i] == -1)
			vertexData[i].collapsable = true;
	}*/

	for(int i = 0; i < vertexData.size(); i++)
	{
		if(collapsed[i] == -1 && vertexData[i].collapsable)
		{
			//find valid edge
			if(vertexData[i].boundary)
			{
				HE_HalfEdge *collapse_me = find_collapsable_edge_boundary( vertexData[i] );
				bool didc = collapse(collapse_me);

				if (didc)
				{
					HE_Vertex *v1 = collapse_me->v;
					HE_Vertex *v2 = collapse_me->next->next->v;

					vertexData[i].collapsable = false;
				//	v1->collapsable = false;
				//	v2->collapsable = false;
				//	set_one_ring_uncollapsable(vertexData[i]);
				//	set_one_ring_uncollapsable(*vertexData[i].he->v);

					collapsed[i] = 1;// this->collapse_iter;
				}
				break;
			}
			else
			{
				//Collapse vert i int one_ring[i][0]
			//	printf("interior\n");
				if(vertexData[i].he->flip->v == &vertexData[i])
				{
					printf("check");
				}
				HE_HalfEdge *he = vertexData[i].he;
				bool didc = collapse(vertexData[i].he->flip);
				if (didc)
				{
				//	HE_Vertex *v1 = vertexData[i].he->e->he->v;
					//HE_Vertex *v2 = vertexData[i].he->e->he->next->next->v;
					//	if(vertexData[i].he->flip == NULL)

					vertexData[i].collapsable = false;
					he->v->collapsable = false;
					set_one_ring_uncollapsable(he->v);
				//	set_one_ring_uncollapsable(&vertexData[i]);
					//v1->collapsable = false;
					//v2->collapsable = false;
					//set_one_ring_uncollapsable(vertexData[i]);
					//set_one_ring_uncollapsable(*vertexData[i].he->v);

					//vertexData[i].when_col = this->collapse_iter;

					//break;
					collapsed[i] = 1;//this->collapse_iter;
				}
				break;
			}
		}
	}

	collapse_iter++;
}

#else
void HalfEdgeMesh::set_one_ring_uncollapsable(HE_Vertex v)
{
	v.collapsable = false;
	HE_HalfEdge *now = v.he;
	now = now->flip;

	
	HE_Vertex *start;

	int count = 0;
	bool goback = false;

	if (now == NULL)
	{
		now = v.he->next;
		start = now->v;
		start->collapsable = false;
	}
	else
	{
		now->v->collapsable = false;
		now->flip->v->collapsable = false;
		
		start = now->v;
		start->collapsable = false;
		
	}

	count = 0;
	goback = false;
	do
	{
		if (now->next->next->flip == NULL)
		{
			goback = true;
			now->next->v->collapsable = false;
			break;
		}
		else
		{
			now = now->next->next->flip;
			now->v->collapsable = false;
		}
		count++;

	} while (now->v != start && count < 100);
	
	count = 0;
	if(goback)
	{
		now = v.he->next;
		do
		{
			if(now->flip == NULL)
			{
				now->v->collapsable = false;
				break;
			}
			else
			{
				now->v->collapsable = false;
				now = now->flip->next;
			}
			count++;

		}while( count < 100 );
	}
	if(count >90)
	printf("shit set un = %d\n", count);
}

HE_Edge* HalfEdgeMesh::find_collapsable_edge_boundary(HE_Vertex v)
{
	HE_HalfEdge *now = v.he;

	now = now->flip;
	if (now == NULL)
		return v.he->e;

	int count = 0;
	do
	{
		if(now->prev->flip == NULL)
		{
			return now->prev->e;
		}
		now = now->prev->flip;
		count++;
	}while(  count < 100 );

	//int count = 0;
	//do
	//{
	//	if( now->next->next->flip == NULL )
	//	{
	//		return now->next->next->e;
	//	}
	//	now = now->next->next->flip;

	///*	if( now->v->boundary )
	//	{
	//		return now->e;
	//	}*/
	//	count++;

	//}while( now->v != start && count < 100 );

	printf("shit find boundaryedge = %d\n", count);
	return v.he->e;
}

void HalfEdgeMesh::finish_collapsing()
{
	int ind = 0;
	newindex.clear();
	for (int i = 0; i < totalVerts; i++)
	{
		if (collapsed[i] == -1)
		{
			vertexData[i].newindex = ind;
			newindex.push_back(ind);
			ind++;
		}
	}

	stack<int> restack;
	int si = col_verts.size();
	for (int i = 0; i < si; i++)
	{
		newindex.push_back(ind);
		vertexData[col_verts.top()].newindex = ind;
		restack.push(col_verts.top());
		col_verts.pop();
		ind++;
	}

	indchange.clear();
	for (int i = 0; i < totalVerts; i++)
	{
		indchange.push_back(vertexData[i].newindex);
	}

	si = restack.size();
	for (int i = 0; i < si; i++)
	{
		col_verts.push( restack.top() );
		restack.pop(); 
	}
	/*for (int i = 0; i < totalVerts; i++)
	{
	printf("%d,  %d -> %d\n", i,vertexData[i].index, vertexData[i].newindex);
	}*/

	

}

vect2d calc_barys( vect2d a, vect2d b, vect2d c, vect2d p)
{
	vect2d ret;
	ret[0] = ((b[1]-c[1])* (-c[0]+p[0])+(-b[0]+c[0])* (-c[1]+p[1]))/((-b[0]+c[0])* (a[1]-c[1])+(a[0]-c[0])* (b[1]-c[1]));
	ret[1] = ((-a[1]+c[1])* (-c[0]+p[0])+(a[0]-c[0])* (-c[1]+p[1]))/((-b[0]+c[0])* (a[1]-c[1])+(a[0]-c[0])*(b[1]-c[1]));
	return ret;
}

vect2d calc_barys( /*wing0*/vect3d b, /*remaining*/vect3d a, /*wing1*/vect3d c, /*collapsing*/vect3d p)
//vect2d calc_barys( /*wing0*/vect3d a, /*remaining*/vect3d b, /*wing1*/vect3d c, /*collapsing*/vect3d p)
{
	vect3d V1 = b-a;
	vect3d V2 = c-a;
	vect3d V3 = p-a;

	double A1 = (V1.cross(V3)).length() * .5;
	double A2 = (V2.cross(V3)).length() * .5;
	
	double A = (1.0f / (V3.length2()) ) * ( A1*(V2.dot(V3)) + A2*(V1.dot(V3)) );

	vect2d ret;
	ret[0] = A2 / A;
	ret[1] = ( A - A1 - A2 ) / A;
	return ret;
}

void HalfEdgeMesh::gen_collapse_lvl()
{
	CollapseLevel lvl;;
	vect3d X, Y, Z;
	for (int i = 0; i < totalFaces; i++)
	{
		if (!faceData[i].hole)
		{
			vector<vect2d> t;
			vect2d t1, t2, t3;
			t.push_back(t1);
			t.push_back(t2);
			t.push_back(t3);

			/*vect3d a = vertexData[faceData[i].vi[0]].pos;
			vect3d b = vertexData[faceData[i].vi[1]].pos;
			vect3d c = vertexData[faceData[i].vi[2]].pos;*/

			vect3d a = faceData[i].he->v->pos;
			vect3d b = faceData[i].he->next->v->pos;
			vect3d c = faceData[i].he->next->next->v->pos;

			X = b - a;
			X.normalize();
			Z = X % (c - a);
			Z.normalize();

			Y = Z % X;

			//// scale triangle by its initial scale (to get approximately the same size mesh as the input)
			double amb = 1;
			X *= amb;
			Y *= amb;

			// store
			t[0].set(0, 0);
			t[1].set((b - a) * X, 0);
			t[2].set((c - a) * X, (c - a) * Y);

			lvl.iso.push_back(t);
		}
	}
	//printf("\tnum iso tris = %d\n", lvl.iso.size());
	//boundary
	for (int i = 0; i < boundaries.size(); i++)
	{
		vector<int> t;
		for (int j = 0; j < boundaries[i].size(); j++)
		{
			if (collapsed[boundaries[i][j]->index] == -1)
			{
				t.push_back(boundaries[i][j]->index);
			}
		}
		lvl.boundaries.push_back(t);
	}
	printf("\tnum boundaries = %d\n", lvl.boundaries.size());
	for (int i = 0; i < lvl.boundaries.size(); i++)
	{
		printf("\t\t %d\n", lvl.boundaries[i].size());
	}

	//faces
	for (int i = 0; i < totalFaces; i++)
	{
		if (!faceData[i].hole)
		{
			vect3i t;
			t[0] = faceData[i].he->v->index;
			t[1] = faceData[i].he->next->v->index;
			t[2] = faceData[i].he->next->next->v->index;

			lvl.faces.push_back(t);
		}
	}

	//one ring
	for (int i = 0; i < totalVerts; i++)
	{
		if (collapsed[i] == -1)
		{
			HE_HalfEdge *start;
			start = vertexData[i].he;
			vector<int> or;
			if (vertexData[i].boundary)
			{
				//find start
				while (start->next->flip != NULL)
				{
					start = start->next->flip;
				}			
				start = start->next;				
				or.push_back( start->v->index );

				//wrap to end
				while (start->prev->flip != NULL)
				{
					start = start->prev->flip;
					or.push_back(start->v->index);
				}
				or.push_back(start->next->v->index);
			}
			else
			{
				//wrap around to start
				HE_HalfEdge *now = start;
				do
				{
					now = now->flip;
					or.push_back( now->v->index );
					now = now->next->next;
				} while (now != start);

			}
			lvl.onerings.push_back(or);
		}
	}
	//printf("\tonerings = %d\n", lvl.onerings.size());
	//push lvl
	col_lev.push_back(lvl);

	collapse_iter++;
}
void HalfEdgeMesh::collapse_one_iteration ( Array<vect3i> indices_pos )
{
	for (int i = 0; i < totalVerts; i++)
	{
		if (collapsed[i] == -1)
			vertexData[i].collapsable = true;
	}

	CollapseLevel lvl;
	lvl.num_col = 0;

	//sanity();
	for(int i = 0; i < vertexData.size(); i++)
	{
		if (collapsed[i] == -1 && vertexData[i].collapsable)
		{
			vect2d p = vertexData[i].uvpos;
			vect3d P = vertexData[i].pos;
			//find valid edge
			if (vertexData[i].boundary)
			{
				if (current_boundary_size[vertexData[i].which_boundary] > MIN_BOUNDARY_SIZE)
					//if(vertexData[i].when_col == -1 && isboundary(&vertexData[i]))
				{
					HE_Edge *collapse_me = find_collapsable_edge_boundary(vertexData[i]);
					vect3i wings;
					wings[0] = collapse_me->he->next->v->index;
				
					//Find next boundary
					HE_HalfEdge *now1 = collapse_me->he;
					while (now1->next->flip != NULL)
					{
						now1 = now1->next->flip;
					}
					wings[1] = now1->next->v->index;

					//wings[1] = -1;
					wings[2] = collapse_me->he->prev->v->index;
				
					vector<vect2i> wTow;
					int prev = wings[0];

					//HE_HalfEdge *now = vertexData[i].he->next;
					HE_HalfEdge *now = collapse_me->he->next;
					while (now->flip != NULL)
					{
						vect2i addme;
						addme[0] = prev;
						now = now->flip->next;
						addme[1] = now->v->index;
						wTow.push_back(addme);
						prev = addme[1];

					}

					vect3i wings2;
					vector<vect2i> wTow2;

					wings2[0] = wings[0];
					wings2[2] = wings[2];

					now1 = collapse_me->he;
					while (now1->prev->flip != NULL)
					{
						now1 = now1->prev->flip;
					}
					wings2[1] = now1->next->v->index;

					prev = wings[0];
					//HE_HalfEdge *now = vertexData[i].he->next;
					now = collapse_me->he->prev;
					while (now->flip != NULL)
					{
						vect2i addme;
						addme[0] = prev;
						now = now->flip->prev;
						addme[1] = now->prev->v->index;
						wTow2.push_back(addme);
						prev = addme[1];
					}



					bool didc = collapse(collapse_me, &vertexData[i]);


					if (didc)
					{

						HE_Vertex *v1 = collapse_me->he->v;
						HE_Vertex *v2 = collapse_me->he->next->next->v;

						vertexData[i].collapsable = false;
						v1->collapsable = false;
						v2->collapsable = false;
						set_one_ring_uncollapsable(*v2);
						set_one_ring_uncollapsable(*v1);

						current_boundary_size[vertexData[i].which_boundary]--;
						collapsed[i] = 1;// this->collapse_iter;
						col_verts.push(i);
						//printf("b%d, ", i);
						lvl.num_col++;

						lvl.wings.push(wings2);
						lvl.wings.push(wings);

						lvl.barys.push( calc_barys(vertexData[wings2[0]].pos, vertexData[wings2[2]].pos, vertexData[wings2[1]].pos ,P) );
						lvl.barys.push( calc_barys(vertexData[wings[0]].pos, vertexData[wings[2]].pos, vertexData[wings[1]].pos ,P) );

						//lvl.barys.push( calc_barys(vertexData[wings2[0]].uvpos, vertexData[wings2[2]].uvpos, vertexData[wings2[1]].uvpos ,p) );
						//lvl.barys.push( calc_barys(vertexData[wings[0]].uvpos, vertexData[wings[2]].uvpos, vertexData[wings[1]].uvpos ,p) );
						
						lvl.wing_to_wing.push(wTow2);
						lvl.wing_to_wing.push(wTow);
						
					}

				}
			}
			else
			{
				vect3i wings;
				wings[0] = vertexData[i].he->next->v->index;
				wings[1] = vertexData[i].he->flip->next->v->index;
				wings[2] = vertexData[i].he->prev->v->index;

				vector<vect2i> wTow;
				int prev = wings[0];

				HE_HalfEdge *now = vertexData[i].he->next;
				while (now->v->index != wings[1])
				{
					if (now->flip == NULL)
					{
						printf("SHIT assumption wrong\n");
						break;
					}
					vect2i addme;
					addme[0] = prev;
					now = now->flip->next;
					addme[1] = now->v->index;
					wTow.push_back(addme);
					prev = addme[1];

				}

				bool didc = collapse(vertexData[i].he->e, &vertexData[i]);
				if (didc)
				{
					HE_Vertex *v1 = vertexData[i].he->e->he->v;
					HE_Vertex *v2 = vertexData[i].he->e->he->next->next->v;

					vertexData[i].collapsable = false;
					v1->collapsable = false;
					v2->collapsable = false;
					set_one_ring_uncollapsable(*v2);
					set_one_ring_uncollapsable(*v1);

					collapsed[i] = 1;//this->collapse_iter;
					col_verts.push(i);
					lvl.num_col++;
					lvl.wings.push(wings);
					lvl.wing_to_wing.push(wTow);
					lvl.barys.push( calc_barys(vertexData[wings[0]].pos, vertexData[wings[2]].pos, vertexData[wings[1]].pos ,P) );

				//	lvl.barys.push( calc_barys(vertexData[wings[0]].uvpos, vertexData[wings[2]].uvpos, vertexData[wings[1]].uvpos ,p) );

				//	printf("%d, ", i);
				}
			}
		}
	}
//	sanity();

	printf("num collapses = %d and wings = %d\n", lvl.num_col, lvl.wings.size());
	//Calc iso
	vect3d X, Y, Z;
	for (int i = 0; i < totalFaces; i++)
	{
		if (!faceData[i].hole)
		{
			vector<vect2d> t;
			vect2d t1, t2, t3;
			t.push_back(t1);
			t.push_back(t2);
			t.push_back(t3);

			vect3d a = faceData[i].he->v->pos;
			vect3d b = faceData[i].he->next->v->pos;
			vect3d c = faceData[i].he->next->next->v->pos;

			X = b - a;
			X.normalize();
			Z = X % (c - a);
			Z.normalize();

			Y = Z % X;

			//// scale triangle by its initial scale (to get approximately the same size mesh as the input)
			double amb = 1;
			X *= amb;
			Y *= amb;

			// store
			t[0].set(0, 0);
			t[1].set((b - a) * X, 0);
			t[2].set((c - a) * X, (c - a) * Y);

			lvl.iso.push_back(t);
		}
	}
	//printf("\tnum iso tris = %d\n", lvl.iso.size());
	//boundary
	for (int i = 0; i < boundaries.size(); i++)
	{
		vector<int> t;
		for (int j = 0; j < boundaries[i].size(); j++)
		{
			if (collapsed[boundaries[i][j]->index] == -1)
			{
				t.push_back(boundaries[i][j]->index);
			}
		}
		lvl.boundaries.push_back(t);
	}
	printf("\tnum boundaries = %d\n", lvl.boundaries.size());
	for (int i = 0; i < lvl.boundaries.size(); i++)
	{
		printf("\t\t %d\n", lvl.boundaries[i].size());
	}

	//faces
	for (int i = 0; i < totalFaces; i++)
	{
		if (!faceData[i].hole)
		{
			vect3i t;
			t[0] = faceData[i].he->v->index;
			t[1] = faceData[i].he->next->v->index;
			t[2] = faceData[i].he->next->next->v->index;

			lvl.faces.push_back(t);
		}
	}

	//one ring
	for (int i = 0; i < totalVerts; i++)
	{
		if (collapsed[i] == -1)
		{
			HE_HalfEdge *start;
			start = vertexData[i].he;
			vector<int> or;
			if (vertexData[i].boundary)
			{
				//find start
				while (start->next->flip != NULL)
				{
					start = start->next->flip;
				}			
				start = start->next;				
				or.push_back( start->v->index );

				//wrap to end
				while (start->prev->flip != NULL)
				{
					start = start->prev->flip;
					or.push_back(start->v->index);
				}
				or.push_back(start->next->v->index);
			}
			else
			{
				//wrap around to start
				HE_HalfEdge *now = start;
				do
				{
					now = now->flip;
					or.push_back( now->v->index );
					now = now->next->next;
				} while (now != start);

			}
			lvl.onerings.push_back(or);
		}
	}
	//printf("\tonerings = %d\n", lvl.onerings.size());
	//push lvl
	col_lev.push_back(lvl);

	collapse_iter++;
	printf("%d,  verts = %d\n", collapse_iter, currentVerts);
}
#endif
bool HalfEdgeMesh::isboundary(HE_Vertex *v )
{
	HE_HalfEdge *he = v->he;
	if(he->flip == NULL)
		return true;

	HE_HalfEdge *start = he;
	he = he->next->flip;
	
	int count = 0;
	while( start != he &&  count < 100)
	{
		if(he == NULL)
			return true;

		he = he->next->flip;
		count++;
	}
	;
	if(count == 100)
		printf("wtf\n");
	return false;
}

bool HalfEdgeMesh::isboundary(HE_HalfEdge *he)
{
	HE_HalfEdge *now = he;
	HE_HalfEdge *start = now;

	int count = 0;
	do
	{
		if (now->next->flip == NULL)
		{
			return true;
		}
		else
		{
			now = now->next->flip;
		}
		count++;

	} while (now != start && count < 100);
	if(count > 90)
		printf("shit %d\n", count);
	now = he;
	now = now->flip;
	start = now;

	count = 0;
	do
	{
		if (now->next->flip == NULL)
		{
			return true;
		}
		else
		{
			now = now->next->flip;
		}
		count++;

	} while (now != start && count < 100);

	return false;
}
#define FLIP_VERT_POINT
// pointer to e no longer valid after call to this function

void HalfEdgeMesh::save_mesh()
{
	printf("saving %s\n", "halfmesh_out.obj");

	FILE *f = fopen("halfmesh_out.obj", "wb");
	vector<int> id;

	int vind = 0;
	for (int i = 0; i < totalVerts; i++)
	{
		if(collapsed[i] == -1)
		{
			id.push_back(vind);
			vind++;
			fprintf(f, "v %f %f %f\n", vertexData[i].pos[0], vertexData[i].pos[1], vertexData[i].pos[2]);
		}
		else
			id.push_back(-1);
	}
	
	for (int i = 0; i < totalVerts; i++)
	{
		if(collapsed[i] == -1)
			fprintf(f, "vt %f %f\n", vertexData[i].uvpos[0], vertexData[i].uvpos[1]);
	}
	
	
	for (int i = 0; i < col_lev[col_lev.size()-1].faces.size(); i++)
	{
		fprintf(f, "f %d/%d %d/%d %d/%d\n", id[col_lev[col_lev.size()-1].faces[i][0]]+1, id[col_lev[col_lev.size()-1].faces[i][0]]+1,
											id[col_lev[col_lev.size()-1].faces[i][1]]+1, id[col_lev[col_lev.size()-1].faces[i][1]]+1,
											id[col_lev[col_lev.size()-1].faces[i][2]]+1, id[col_lev[col_lev.size()-1].faces[i][2]]+1);
	}

	fclose(f);
}

void HalfEdgeMesh::save_mesh2()
{
	printf("saving %s\n", "halfmesh_out.obj");

	FILE *f = fopen("halfmesh_out.obj", "wb");
	vector<int> id;

	int vind = 0;
	for (int i = 0; i < totalVerts; i++)
	{
		if(collapsed[i] == -1 || true)
		{
			id.push_back(vind);
			vind++;
			fprintf(f, "v %f %f %f\n", vertexData[i].pos[0], vertexData[i].pos[1], vertexData[i].pos[2]);
		}
		else
			id.push_back(-1);
	}
	
	for (int i = 0; i < totalVerts; i++)
	{
		//if(collapsed[i] == -1)
			fprintf(f, "vt %f %f\n", vertexData[i].uvpos[0], vertexData[i].uvpos[1]);
	}
	
	
	for (int i = 0; i < col_lev[1].faces.size(); i++)
	{
		fprintf(f, "f %d/%d %d/%d %d/%d\n", id[col_lev[1].faces[i][0]]+1, id[col_lev[1].faces[i][0]]+1,
											id[col_lev[1].faces[i][1]]+1, id[col_lev[1].faces[i][1]]+1,
											id[col_lev[1].faces[i][2]]+1, id[col_lev[1].faces[i][2]]+1);
	}

	fclose(f);
}
#ifdef SCOTT_COLLAPSE
bool HalfEdgeMesh::collapse ( HE_HalfEdge *he )
{
	if(he->v->boundary && he->prev->v->boundary)
	{//if two boundary vertices

		if(he->flip != NULL)
			return false;//boundary to boundary not flipped

		if(he->next->flip == NULL && he->prev->flip == NULL)
			return false;//Would create dangline line

	}

	if (he->flip == NULL)
	{
		//Boundary to boundary
		he->v->boundary = true;
		he->next->next->v->boundary = true;
	}
	else if( he->v->boundary || he->next->next->v->boundary)
	{
		//Interior to Boundary
		he->v->boundary = true;
		he->next->next->v->boundary = true;
	}



// update vertex position
    he->v->pos = he->next->next->v->pos;
 
    // common to all cases, we must make edges point from old vertex to new one
    HE_HalfEdge *curr = he;
    if ( curr->flip != NULL )
    {
            curr = curr->flip->next->next; // advance beyond wing polygon, not absolutely necessary
    }
	int count = 0;
    while ( curr->flip != NULL && curr->flip->next->next != he )
    {
            curr = curr->flip->next->next;
            curr->v = he->next->next->v;
			count++;
			if(count > 100)
			{
				printf("SHIT\n");
			}
    }
    if ( curr->flip == NULL ) // hit a boundary, go the other way until we hit another boundary
    {
            curr = he;
            while ( curr->next->flip != NULL )
            {
                    curr = curr->next->flip;
                    curr->v = he->next->next->v;
            }
    }
 
    // fix new vertex
    if ( he->next->flip != NULL )
    {
            he->v->he = he->next->flip->next;
    }
    else
    { // he->next->next->flip must exist
            he->v->he = he->next->next->flip;
    }
 
    // fix surviving half-edges
    if ( he->next->flip != NULL )
    {
            he->next->v->he = he->next->flip; // fix wing vertex
            he->next->flip->flip = he->next->next->flip;
    }
    if ( he->next->next->flip != NULL )
    {
            he->next->v->he = he->next->next->flip->next; // fix wing vertex
            he->next->next->flip = he->next->flip;
    }
    he->f->hole = true;
    if ( he->flip != NULL )
    {
            if ( he->flip->next->flip != NULL )
            {
                    he->flip->next->v->he = he->flip->next->flip; // fix wing vertex
                    he->flip->next->flip->flip = he->flip->next->next->flip;
            }
            if ( he->flip->next->next->flip != NULL )
            {
                    he->flip->next->v->he = he->flip->next->next->flip->next; // fix wing vertex
                    he->flip->next->next->flip->flip = he->flip->next->flip;
            }
            he->flip->f->hole = true;
    }
 
    currentVerts--;

}
#else

void HalfEdgeMesh::ring_and_spoke( HE_Vertex *v )
{
	HE_HalfEdge *start = vertexData [ v->index ].he, *curr;

	bool isbound = vertexData[v->index].boundary;

	//if boundary vertex move start to left triangle
	if(isbound)
	{
		while(start->next->flip != NULL)
			start = start->next->flip;
	}

	curr = start;
	bool cont = true;
	do
	{
		//outer ring
		errQueue.remove(curr->next->next->e);
		errQueue.insert ( curr->next->next->e, minimizeEdge ( curr->next->next->e ) );

		//spoke
		errQueue.remove(curr->next->e);
		errQueue.insert ( curr->next->e, minimizeEdge ( curr->next->e ) );

		if(isbound)
		{
			if(curr->flip != NULL)
			{
				curr = curr->flip->next->next;
				cont = true;
			}
			else
				cont = false;
		}
		else
		{
			curr = curr->flip->next->next;
			cont = curr != start;
		}
	} while ( cont );

	if(isbound)
	{
		errQueue.remove(curr->e);
		errQueue.insert ( curr->e, minimizeEdge ( curr->e ) );
	}
}
void HalfEdgeMesh::collapseMinimal ( void )
{
	static int counter = 0;
	counter++;
	HE_Edge *e;
	float err = FLT_MAX;
	bool reverse = false;
	do
	{
		err = errQueue.peekNextPriority ( );
		for ( e = (HE_Edge *)errQueue.serveNext ( ); e->dirty; e = (HE_Edge *)errQueue.serveNext ( ) )
		{ 
			errQueue.insert ( e, minimizeEdge ( e ) );
			err = errQueue.peekNextPriority ( );
		}
		if ( isSafeCollapse ( e ) )
		{
				HE_Vertex *vv = e->he->v;
				HE_Vertex *keepme = e->he->next->next->v;
				
				if(vv->boundary && !e->he->next->next->v->boundary)
				{
					keepme = vv;
					vv = e->he->next->next->v;

					e->he = e->he->flip;
				}

				//printf("Collapse V = %d\n", vv->index);
				

				//Build wings and wing2wing
				vect3i wings;
				vector<vect2i> wTow;

				vect3i wings2;
				vector<vect2i> wTow2;

				if(keepme->boundary && vv->boundary)
				{
					//boundary edge
					wings[0] = e->he->next->v->index;

					HE_HalfEdge *now1 = e->he;
					while (now1->next->flip != NULL)
					{
						now1 = now1->next->flip;
					}
					wings[1] = now1->next->v->index;

					wings[2] = e->he->prev->v->index;

					
					int prev = wings[0];
					//HE_HalfEdge *now = vertexData[i].he->next;
					HE_HalfEdge *now = e->he->next;
					while (now->flip != NULL)
					{
						vect2i addme;
						addme[0] = prev;
						now = now->flip->next;
						addme[1] = now->v->index;
						wTow.push_back(addme);
						prev = addme[1];

					}

					wings2[0] = wings[0];
					wings2[2] = wings[2];

					now1 = e->he;
					while (now1->prev->flip != NULL)
					{
						now1 = now1->prev->flip;
					}
					wings2[1] = now1->next->v->index;

					prev = wings[0];
					//HE_HalfEdge *now = vertexData[i].he->next;
					now = e->he->prev;
					while (now->flip != NULL)
					{
						vect2i addme;
						addme[0] = prev;
						now = now->flip->prev;
						addme[1] = now->prev->v->index;
						wTow2.push_back(addme);
						prev = addme[1];
					}
				}
				else
				{
					wings[0] = e->he->next->v->index;
					wings[1] = e->he->flip->next->v->index;
					wings[2] = e->he->prev->v->index;

					int prev = wings[0];

					HE_HalfEdge *now = e->he->next;
					while (now->v->index != wings[1])
					{
						if (now->flip == NULL)
						{
							printf("SHIT assumption wrong\n");
							break;
						}
						vect2i addme;
						addme[0] = prev;
						now = now->flip->next;
						addme[1] = now->v->index;
						wTow.push_back(addme);
						prev = addme[1];

					}
				}

				col_lev[0].p_col.push(vv->pos);
				col_lev[0].p_rem.push(keepme->pos);
				col_lev[0].uv_p_col.push( vv->uvpos );;
				col_lev[0].uv_p_rem.push( keepme->uvpos);
				vect2i tt;
				tt[0] = vv->index;
				tt[1] = keepme->index;
				col_lev[0].edge.push(tt);


				collapse( e, vv );
				//remove vertex

				col_verts.push( vv->index );
				col_lev[0].num_col++;

				vect3d P = vv->pos;
				if(keepme->boundary && vv->boundary)
				{
					col_lev[0].wings.push(wings2);
					col_lev[0].wing_to_wing.push(wTow2);
					col_lev[0].barys.push( calc_barys(vertexData[wings2[0]].pos, vertexData[wings2[2]].pos, vertexData[wings2[1]].pos ,P) );

					current_boundary_size[vv->which_boundary]--;
				}

				col_lev[0].barys.push( calc_barys(vertexData[wings[0]].pos, vertexData[wings[2]].pos, vertexData[wings[1]].pos ,P) );
				col_lev[0].wings.push(wings);
				col_lev[0].wing_to_wing.push(wTow);
				//col_lev[0].barys.push( calc_barys(vertexData[wings[0]].pos, vertexData[wings[2]].pos, vertexData[wings[1]].pos ,P) );


				//if(vv->boundary && keepme->boundary)
				//{
				//	collapsed[keepme->index] = 1;
				//	//vertexData[vv->index].collapsable = false;
				//	//walk around keepme and set to rebuild qefs
				//	set_to_recalc_qef(vv);

				//	//rebuild for vertices
				//	buildQEFs();

				//	//go through edges... one ring and spokes and remove, recalc, insert
				//	ring_and_spoke(vv);
				//}
				//else
				{
					
					collapsed[vv->index] = 1;
					//vertexData[vv->index].collapsable = false;
					//walk around keepme and set to rebuild qefs
					set_to_recalc_qef(keepme);

					//rebuild for vertices
					buildQEFs();

					//go through edges... one ring and spokes and remove, recalc, insert
					ring_and_spoke(keepme);
				}


			/*}
			else
				reverse = true;*/
		//	printf("curr = %d\n", currentVerts);
			break;
			
		}
		e->index = -3; 
	} while ( true ); // keep removing edges until one is safe

	if ( counter % 1000 == 0 )
		printf ( "error = %g\n", err );
	
	

	/*bool collapsed = collapse ( e, e->he->v );
	if(!collapsed)
	{
		collapsed = collapse ( e, e->he->next->next->v );
	}

	if(collapsed)
	{
	}*/
}

void HalfEdgeMesh::sanity()
{
	int i;
	for (i = 0; i < totalVerts; i++)
	{
		if (!vertexData[i].boundary && collapsed[i] != 1)
		{
			HE_HalfEdge *he = vertexData[i].he;

			HE_HalfEdge *now = he;
			int count = 0;
			do
			{
				count++;
				if (now != now->flip->flip)
					printf("shitttttyy");

				now = now->flip->prev;
			} while (he != now && count < 1000);
			if (count == 1000)
				printf("shitty\n");
		}
	}


	//Check valence 2
	/*for(i = 0; i < totalVerts; i++)
	{
		if(!vertexData[i].boundary && collapsed[i] != 1)
		{
			int val = 0;
			vertexData[i].he
		}
	}*/
}

void HalfEdgeMesh::combineQEFs ( HE_Vertex *v1, HE_Vertex *v2 )
{
	#ifdef TRACE_DEBUG
	fprintf(trace_file1, "combine verts %d, %d\n", v1 - &vertexData[0], v2 - &vertexData[0]);
#endif
	v1->qef.combineSelf ( v2->qef );
}

float HalfEdgeMesh::minimizeEdge ( HE_Edge *e )
{
	double error = 0;
	int i;
	QEF_TYPE<NUM_TYPE, 3> qef;

	e->dirty = false;

	e->pos = ( e->he->v->pos + e->he->next->next->v->pos ) / 2.0f;

	qef = e->he->v->qef;
	qef.combineSelf ( e->he->next->next->v->qef );
	double guess [ 3 ], point [ 3 ];
	guess [ 0 ] = e->pos[0];
	guess [ 1 ] = e->pos[1];
	guess [ 2 ] = e->pos[2];
	#ifdef TRACE_DEBUG
	fprintf(trace_file1, "guess = %I32x %I32x %I32x\n", *(int*)&guess[0], *(int*)&guess[1], *(int*)&guess[2]);
	//fprintf(trace_file1, "guess = %g, %g, %g\n", guess[0], guess[1], guess[2]);
#endif

	qef.calcPoint ( guess, point );
	e->pos[0] = point [ 0 ];
	e->pos[1] = point [ 1 ];
	e->pos[2] = point [ 2 ];
	#ifdef TRACE_DEBUG
	fprintf(trace_file1, "point = %I32x %I32x %I32x\n", *(int*)&point[0], *(int*)&point[1], *(int*)&point[2]);
	//fprintf(trace_file1, "point = %g, %g, %g\n", point[0], point[1], point[2]);
#endif
	error = qef.minimizerError ( point );

	#ifdef TRACE_DEBUG
	fprintf(trace_file1, "error = %I64x\n", *(__int64*)&error);
#endif

	//printf("%f, %f, %f and %f, %f, %f -> %f\n", e->he->v->pos[0], e->he->v->pos[1], e->he->v->pos[2], e->he->flip->v->pos[0], e->he->flip->v->pos[1], e->he->flip->v->pos[2], error);
	//printf("%f, %f, %f\n", e->pos[0], e->pos[1], e->pos[2]);
	return error;
}


bool HalfEdgeMesh::is_3_valence( HE_Vertex *v )
{
	HE_HalfEdge *he = v->he;

	int count = 0;

	do
	{
		if(he->flip == NULL)
			return false;
		he = he->flip;
		count++;

		he = he->next->next;
	}
	while(he != v->he  && count < 100);

	/*do{
		if(he->next->flip == NULL)
		{
			return false;
		}
		count++;
		he = he->next->flip;
	}
	while(he != v->he && count < 100);
*/
	if(count == 100)
		printf("fuck\n");
	if(count == 3)
		return true;

	return false;
}
bool HalfEdgeMesh::collapse ( HE_Edge *e, HE_Vertex *v /*vert to be collapsed*/ )
{
	//v->when_col = 0;
	//if (!isSafeCollapse(e))
	//	return false;

	//printf("qef1 = \n");
	//e->he->v->qef.print_mathmatica();

	//printf("qef2 = \n");
	//e->he->next->next->v->qef.print_mathmatica();

	QEF_TYPE<NUM_TYPE, 3> qef;
	qef = e->he->v->qef;
	qef.combineSelf( e->he->next->next->v->qef );
	//printf("combineqef = \n");
	//qef.print_mathmatica();

	double guess [ 3 ], point [ 3 ];
	guess [ 0 ] = e->pos[0];
	guess [ 1 ] = e->pos[1];
	guess [ 2 ] = e->pos[2];

	qef.calcPoint ( guess, point );
	e->pos[0] = point [ 0 ];
	e->pos[1] = point [ 1 ];
	e->pos[2] = point [ 2 ];

	double error = qef.minimizerError ( point );
	//printf("tp = {%f, %f, %f};\n", point[0], point[1], point[2] );

	/*printf("tri1 = {");
	print_buildQEFs ( e->he->v );
	printf("};\n");

	printf("tri2 = {");
	print_buildQEFs ( e->he->next->next->v );
	printf("};\n");*/

	int i;
	HE_HalfEdge *he = e->he;

	//printf("p1 = {%f,%f,%f};\n", he->v->pos[0], he->v->pos[1], he->v->pos[2] );
	//printf("p2 = {%f,%f,%f};\n", he->next->next->v->pos[0], he->next->next->v->pos[1], he->next->next->v->pos[2] );

	he->v->pos = e->pos;
	he->next->next->v->pos = e->pos;

	//printf("np = {%f,%f,%f};\n", e->pos[0], e->pos[1], e->pos[2] );

	if (he->v == v) //he needs to poinbt in direction of collapse
		he = he->flip;

	//Decide case of edge collapse
	if (e->he->flip == NULL)
	{
		//printf("b->b\n");
		//Case is boundary -> boundary

		/*if (is_3_val(e->he))
			return false;*/

		if (e->he->v == v)
		{
			
			//v collapses reverse he
			HE_HalfEdge *start = e->he;
			//start->v->boundary = true;
			start->v = e->he->prev->v;
			//start->v->boundary = true;
			start = start->next->flip;
			while (start != NULL)
			{
				start->v = e->he->prev->v;
				start = start->next->flip;
			}

			if( e->he->prev->flip != NULL)
				e->he->prev->v->he = e->he->prev->flip->prev;
			else
				e->he->prev->v->he = e->he->next->flip;

			
		}
		else
		{
			printf("shitty\n");
		}

		// triangle, must delete
			// fix half-edges
		if(e->he->next->flip != NULL)
			e->he->next->flip->flip = e->he->prev->flip;
		
		if(e->he->prev->flip != NULL)
			e->he->prev->flip->flip = e->he->next->flip;

			// fix wing-vertex   // Do nt get this... no flip
		if(e->he->prev->flip != NULL)
			e->he->next->v->he = e->he->prev->flip;
		else
			e->he->next->v->he =e->he->next->flip->prev;

		if(e->he->next->flip!= NULL)
			e->he->prev->e->he = e->he->next->flip;
		else
			e->he->prev->e->he = e->he->prev->flip;


	//		 fix edge
		if(e->he->next->flip != NULL)
			e->he->prev->e->he = e->he->next->flip;
		else
			e->he->prev->e->he = e->he->prev->flip;

		//if(e->he->prev->flip != NULL)
		if ( e->he->next->flip != NULL )
		{
			errQueue.remove ( e->he->next->flip->e );
			e->he->next->flip->e = e->he->prev->e;
		}

		//	 mark face as bad
			e->he->f->hole = true;
	}
	else if (/*e->he->next->flip == NULL || e->he->prev->flip == NULL ||*/ isboundary(he) )
	{
	//	printf("i->b\n");
		//Case is interior -> boundary
		HE_HalfEdge *curr = he->flip;


		//printf("set b\n");
		int count = 0;
		do
		{
			curr->v = he->v;
			curr = curr->next->flip;
			count++;
		} while (curr != NULL && curr != he->flip && count < 100);

		if(curr == NULL)
		{
			curr = he->prev;
			curr->v = he->v;
			while(curr->flip  != NULL)
			{
				curr = curr->flip->prev;
				curr->v = he->v;
			}
		}

#ifdef FLIP_VERT_POINT
		//he->v->he = he->flip->prev;   //shouldnt this be he->flip->prev
		// Note: v->he->v == v
		if ( he->next->flip == NULL )
		{
			he->v->he = he->prev->flip->prev; // note this face must exist due to previous check
		}
		else
		{
			he->v->he = he->next->flip;
		}
#else
		// fix new vertex DONT GET THIS LINE
		he->v->he = he->flip->prev->flip;   //shouldnt this be he->flip->prev
#endif

		// triangle, must delete
		// fix half-edges
	//	if(he->prev->flip!= NULL && )
		if ( he->next->flip != NULL )
		{
			he->next->flip->flip = he->prev->flip;
		}
		if ( he->prev->flip != NULL )
		{
			he->prev->flip->flip = he->next->flip;
		}

		// fix wing-vertex  
#ifdef FLIP_VERT_POINT
		if ( he->prev->flip != NULL )
		{
			he->next->v->he = he->prev->flip;
		}
		else
		{
			he->next->v->he = he->next->flip->prev;
		}
#else
		he->next->v->he = he->next->flip;
#endif
		// fix edge
		if ( he->next->flip == NULL )
		{
			he->next->e->he = he->prev->flip;
		}
		else
		{
			he->next->e->he = he->next->flip;
		}
		if ( he->prev->flip != NULL )
		{
			errQueue.remove ( he->prev->flip->e );
			he->prev->flip->e = he->next->e;
		}
		//

		// mark face as bad
		he->f->hole = true;

		// triangle, must delete
		// fix half-edges

		if ( he->flip->next->flip != NULL )
		{
			he->flip->next->flip->flip = he->flip->prev->flip;
		}
		if ( he->flip->prev->flip != NULL )
		{
			he->flip->prev->flip->flip = he->flip->next->flip;
		}

#ifdef FLIP_VERT_POINT
		if ( he->flip->prev->flip != NULL )
		{
			he->flip->next->v->he = he->flip->prev->flip;
		}
		else
		{
			he->flip->next->v->he = he->flip->next->flip->prev;
		}
//		he->flip->next->v->he = he->flip->next->flip->prev;
#else
		he->flip->next->v->he = he->flip->next->flip;
#endif
		
		// fix edge
		if ( he->flip->next->flip == NULL )
		{
			he->flip->next->e->he = he->flip->prev->flip;
		}
		else
		{
			he->flip->next->e->he = he->flip->next->flip;
		}
		if ( he->flip->prev->flip != NULL )
		{
			errQueue.remove ( he->flip->prev->flip->e );
			he->flip->prev->flip->e = he->flip->next->e;
		}

		// fix edge
	//	he->flip->prev->e->he = he->flip->next->flip;
	//	he->flip->next->flip->e = he->flip->prev->e;

		// mark face as bad
		he->flip->f->hole = true;
		
	}
	else
	{
	//	printf("i->i\n");
		//Case is interior -> interior
		HE_HalfEdge *curr = he->flip;

		int count = 0;
		do
		{
			curr->v = he->v;
			curr = curr->next->flip;

			
			count++;
		} while (curr != he->flip && count < 1000);
		if(count >900)
			printf("count = %d\n",count);
		//i -> b   .ok.
		//b -> b
#ifdef FLIP_VERT_POINT
		//he->v->he = he->flip->prev;   //shouldnt this be he->flip->prev
		//he->v->he = he->flip->next->flip;
		he->v->he = he->next->flip;
#else
		// fix new vertex DONT GET THIS LINE
		he->v->he = he->flip->prev->flip;   //shouldnt this be he->flip->prev
#endif
		//i -> b   .X. he->flip->prev->flip == null
		//b -> b

		// triangle, must delete
		// fix half-edges
		he->next->flip->flip = he->prev->flip;
		//i -> b   .X. he->next->flip == null
		//b -> b
		he->prev->flip->flip = he->next->flip;
		//i -> b   .ok. he->next->flip == null
		//b -> b

		// fix wing-vertex   // Do nt get this... no flip
#ifdef FLIP_VERT_POINT
		he->next->v->he = he->prev->flip;
		//he->next->v->he = he->next;
#else
		he->next->v->he = he->next->flip;
#endif
		//i -> b  .X.  he->next->flip == null  should be he->prev //Should be flip???
		//b -> b

		// fix edge
		he->next->e->he = he->next->flip;
		//i -> b .X. he->next->flip == null  should be he->prev->flip
		//b -> b
		errQueue.remove ( he->prev->flip->e );
		he->prev->flip->e = he->next->e;
		//i -> b  .ok.
		//b -> b

		// mark face as bad
		he->f->hole = true;
		//i -> b  
		//b -> b
		/////////////////////////////////////////////////////////////
		// triangle, must delete
		// fix half-edges
		he->flip->next->flip->flip = he->flip->prev->flip;
		//i -> b  .ok. he->flip->prev->flip ==null
		//b -> b
		he->flip->prev->flip->flip = he->flip->next->flip;
		//i -> b  .X. he->flip->prev->flip == null
		//b -> b
		// fix wing-vertex   //same flip problem?
#ifdef FLIP_VERT_POINT
	//	he->flip->next->v->he = he->flip->next;
		he->flip->next->v->he = he->flip->prev->flip;
#else
		he->flip->next->v->he = he->flip->next->flip;
#endif
		//i -> b  .ok.
		//b -> b
		// fix edge
		he->flip->next->e->he = he->flip->next->flip;
		//i -> b  .ok.
		//b -> b
		errQueue.remove ( he->flip->prev->flip->e );
		he->flip->prev->flip->e = he->flip->next->e;
		//i -> b  .X. he->flip->prev->flip == null
		//b -> b
		// mark face as bad
		he->flip->f->hole = true;
		//i -> b  .ok.
		//b -> b
		
	}
	currentVerts--;
	return true;
}
#endif

// pointer to e no longer valid after call to this function
ExpandOp *HalfEdgeMesh::collapse ( HE_Edge *e, int* &head, int &f1, int &f2, int &v )
{
	ExpandOp *rvalue;
//	int i;
//	//HE_HalfEdge *he = e->he;
////	HE_HalfEdge *debugHe = &heData[0];
//	HE_HalfEdge *he = min(e->he, e->he->flip);
//
//
//	//** count faces touching both vertices
//	HE_HalfEdge *curr = he->flip;
//	int numFaces = 0;
//	do
//	{
//		numFaces++;
//		curr = curr->next->flip;
//	} while ( curr != he->flip );
//	int oldFaces = numFaces - 2;
//	curr = he;
//	do
//	{
//		numFaces++;
//		curr = curr->next->flip;
//	} while ( curr != he );
//	numFaces -= 4; // 2 deleted faces counted twice
//
//	if (numFaces < 3)
//	{
//		return 0;
//	}
//
//	v = he->flip->v - &vertexData[0];
//	f1 = he->f - &faceData[0]; // face containing wing1
//	f2 = he->flip->f - &faceData[0]; // face contianing wing2
//	
//#ifdef VERT_TYPE_BLENDING
//	head -= 19 + numFaces*4;
//#else
//#ifdef AFFINE_WEIGHTS
//	head -= 13 + numFaces*4;
//#else
//	head -= 10 + numFaces;
//#endif
//#endif
//	rvalue = (ExpandOp *)head;
//	rvalue->renumSize = oldFaces;
//	rvalue->neighborSize = numFaces;
//	rvalue->vertInd = he->v - &vertexData[0]; // expanding vertex
//	rvalue->wing1 = he->next->v - &vertexData[0];
//	rvalue->wing2 = he->flip->next->v - &vertexData[0];
//	rvalue->p1 [ 0 ] = he->v->pos [ 0 ];
//	rvalue->p1 [ 1 ] = he->v->pos [ 1 ];
//	rvalue->p1 [ 2 ] = he->v->pos [ 2 ];
//	rvalue->p2 [ 0 ] = he->flip->v->pos [ 0 ];
//	rvalue->p2 [ 1 ] = he->flip->v->pos [ 1 ];
//	rvalue->p2 [ 2 ] = he->flip->v->pos [ 2 ];
//
//	/******** store expansion topology *********/ 
//	{
//		HE_ExpandOp *op = expandOps + currentVerts - 1;
//
//		// store old position in expand op
//		op->pos[0] = he->v->pos;
//		op->pos[1] = he->flip->v->pos;
//		op->he[0] = he->next->flip;
//		op->he[1] = he->flip->next->flip;
//
//#ifdef VERT_TYPE_BLENDING
//		rvalue->type1 = he->v->type;
//		rvalue->type2 = he->flip->v->type;
//#endif
//
//		// store removed objects in expand op
//		op->hec[0][0] = op->he[0]->flip;
//		op->hec[0][1] = op->he[0]->flip->next;
//		op->hec[0][2] = op->he[0]->flip->next->next;
//		op->hec[1][0] = op->he[1]->flip;
//		op->hec[1][1] = op->he[1]->flip->next;
//		op->hec[1][2] = op->he[1]->flip->next->next;
//		op->e[0] = he->prev->e;
//		op->e[1] = he->flip->prev->e;
//		op->ec = he->e;
//		op->f[0] = he->f;
//		op->f[1] = he->flip->f;
//		op->v = he->flip->v;
//	}
//
//	// update vertex position/qefs/weights
//#ifdef VERT_TYPE_BLENDING
//#if 1
//	////** simply average type values
//	he->v->type = (he->v->type + he->flip->v->type) / 2;
//#else
//	//** calculate blending value based on distance
//	{
//		vect3f p1(he->v->pos[0], he->v->pos[1], he->v->pos[2]);
//		vect3f p2(he->flip->v->pos[0], he->flip->v->pos[1], he->flip->v->pos[2]);
//		vect3f x(e->pos[0], e->pos[1], e->pos[2]);
//
//		vect3f d = p2 - p1;
//		float len = d.length();
//		d /= len;
//		float n = d * (x-p1);
//		float a = n/len;
//		if (a < 0)
//			a = 0;
//		else if (a > 1)
//			a = 1;
//
//		he->v->type = he->v->type * (1-a) + he->flip->v->type*a;
//	}
//#endif
//#endif
//	he->v->pos = e->pos;
//	//printf("pos collapse %f %f %f\n", e->pos[0], e->pos[1], e->pos[2]);
//
//	combineQEFs ( he->v, he->flip->v );
//	//printf("qef %f %f %f %f %f %f %f %f %f %f\n", he->v->qef.data[0], he->v->qef.data[1], he->v->qef.data[2], he->v->qef.data[3], he->v->qef.data[4], he->v->qef.data[5], he->v->qef.data[6], he->v->qef.data[7], he->v->qef.data[8], he->v->qef.data[9]);
//
//	curr = he->next->next->flip;
//	int ind = 0;
//	do
//	{
//		rvalue->getNeigh() [ ind++ ] = curr->f - &faceData[0];
//		curr = curr->next->next->flip;
//	} while ( curr != he->flip->next );
//
//	curr = he->flip->next->next->flip;
//	do
//	{
//		rvalue->getNeigh() [ ind++ ] = curr->f - &faceData[0];
//		curr = curr->next->next->flip;
//	} while ( curr != he->next );
//
//	// common to all cases, we must make edges point from old vertex to new one
//	curr = he->flip;
//	do
//	{
//		curr->v = he->v;
//		curr = curr->next->flip;
//	} while ( curr != he->flip );
//	// fix new vertex
//	he->v->he = he->flip->prev->flip;
//
//	// fix half-edges
//	he->next->flip->flip = he->prev->flip;
//	he->prev->flip->flip = he->next->flip;
//	// fix wing-vertex
//	he->next->v->he = he->next->flip;
//	// fix edge
//	he->next->e->he = he->next->flip;
//	he->prev->flip->e = he->next->e;
//	// mark face as bad
//	he->f->hole = true;
//
//	// delete extra edge
//	errQueue.remove ( he->prev->e );
//
//	// fix half-edges
//	he->flip->next->flip->flip = he->flip->prev->flip;
//	he->flip->prev->flip->flip = he->flip->next->flip;
//	// fix wing-vertex
//	he->flip->next->v->he = he->flip->next->flip;
//	// fix edge
//	he->flip->next->e->he = he->flip->next->flip;
//	he->flip->prev->flip->e = he->flip->next->e;
//	// mark face as bad
//	he->flip->f->hole = true;
//
//	// delete extra edge
//	errQueue.remove ( he->flip->prev->e );
//
//	// set remaining edges to dirty
//	HE_HalfEdge *start, *current;
//	start = he->v->he;
//	current = start;
//	do
//	{
//		current->e->dirty = true;
//		if ( current->e->index == -3 )
//		{
//			// insert back into priority queue
//			errQueue.insert ( current->e, minimizeEdge ( current->e ) );
//		}
//		if ( current->next->e->index == -3 )
//		{
//			// insert back into priority queue
//			errQueue.insert ( current->next->e, minimizeEdge ( current->next->e ) );
//		}
//		current = current->next->next->flip;
//	} while ( current != start );
//
//	currentVerts--;
//
//#ifdef AFFINE_WEIGHTS // position in terms of weights
//	ExpandOp *op = rvalue;
//
//	current = start;
//	for (int i = 0; i < op->neighborSize; i++)
//	{
//		current->f->he = current; // to ensure correct ordering of the faces found from the face later in the function
//		current = current->next->next->flip;
//	}
//
//	// find vertices that are in the 1-ring
//	start = faceData[op->getNeigh()[0]].he;
//	current = start;
//	vect3f p_cen(0,0,0);
//	vect3f pnorm(0,0,0);
//	float p_area = 0;
//	float pa[100];
//	op->vertCenter = start->flip->v - &vertexData[0];
//	for (int i = 0; i < op->neighborSize; i++)
//	{
//		int idx = current->v - &vertexData[0];
//		op->getVert()[i] = idx;
//		op->getWeight1()[i] = 0;
//		op->getWeight2()[i] = 0;
//
//		assert(current->flip->v - &vertexData[0] == op->vertCenter);
//		current = current->next->next->flip;
//	}
//
//	// find centroid
//	for (int i = 0; i < op->neighborSize; i++)
//	{
//		int neigh = op->getNeigh()[i];
//
//		vect3f p0, p1, p2;
//		vec3Copy(p0, *(vect3d*)&(faceData[neigh].he->v->pos));
//		vec3Copy(p1, *(vect3d*)&(faceData[neigh].he->next->v->pos));
//		vec3Copy(p2, *(vect3d*)&(faceData[neigh].he->next->next->v->pos));
//
//		vect3f norm = (p1-p0)%(p2-p0);
//		pnorm += norm;
//		pa[i] = norm.length()*.5;
//		p_area += pa[i];
//
//		p_cen += (p0 + p1 + p2)*pa[i];
//	}
//	pnorm /= sqrt ( pnorm.length ( ) );
//	p_cen /= 3*p_area;
//
//	// find weights of split vertices
//	float ATA[9];
//	for (int i = 0; i < 9; i++)
//	{
//		ATA[i] = 0;
//	}
//
//	for (int i = 0; i < op->neighborSize; i++)
//	{
//		int neigh = op->getNeigh()[i];
//
//		vect3f p0, p1, p2;
//		vec3Copy(p0, *(vect3d*)&(faceData[neigh].he->v->pos));
//		vec3Copy(p1, *(vect3d*)&(faceData[neigh].he->next->v->pos));
//		vec3Copy(p2, *(vect3d*)&(faceData[neigh].he->next->next->v->pos));
//
//		float p00 = p0[0] - p_cen[0];
//		float p01 = p0[1] - p_cen[1];
//		float p02 = p0[2] - p_cen[2];
//		float p10 = p1[0] - p_cen[0];
//		float p11 = p1[1] - p_cen[1];
//		float p12 = p1[2] - p_cen[2];
//		float p20 = p2[0] - p_cen[0];
//		float p21 = p2[1] - p_cen[1];
//		float p22 = p2[2] - p_cen[2];
//
//		ATA[0] += ((1.0)/(24.0)*(((p00)*(((2.0)*(p00))+(p10)+(p20)))+((p10)*((p00)+((2.0)*(p10))+(p20)))+((p20)*((p00)+(p10)+((2.0)*(p20))))))*pa[i];
//		ATA[1] += ((1.0)/(24.0)*(((p00)*(((2.0)*(p01))+(p11)+(p21)))+((p10)*((p01)+((2.0)*(p11))+(p21)))+((p20)*((p01)+(p11)+((2.0)*(p21))))))*pa[i];
//		ATA[2] += ((1.0)/(24.0)*(((p00)*(((2.0)*(p02))+(p12)+(p22)))+((p10)*((p02)+((2.0)*(p12))+(p22)))+((p20)*((p02)+(p12)+((2.0)*(p22))))))*pa[i];
//		ATA[3] += ((1.0)/(24.0)*(((p01)*(((2.0)*(p00))+(p10)+(p20)))+((p11)*((p00)+((2.0)*(p10))+(p20)))+((p21)*((p00)+(p10)+((2.0)*(p20))))))*pa[i];
//		ATA[4] += ((1.0)/(24.0)*(((p01)*(((2.0)*(p01))+(p11)+(p21)))+((p11)*((p01)+((2.0)*(p11))+(p21)))+((p21)*((p01)+(p11)+((2.0)*(p21))))))*pa[i];
//		ATA[5] += ((1.0)/(24.0)*(((p01)*(((2.0)*(p02))+(p12)+(p22)))+((p11)*((p02)+((2.0)*(p12))+(p22)))+((p21)*((p02)+(p12)+((2.0)*(p22))))))*pa[i];
//		ATA[6] += ((1.0)/(24.0)*(((p02)*(((2.0)*(p00))+(p10)+(p20)))+((p12)*((p00)+((2.0)*(p10))+(p20)))+((p22)*((p00)+(p10)+((2.0)*(p20))))))*pa[i];
//		ATA[7] += ((1.0)/(24.0)*(((p02)*(((2.0)*(p01))+(p11)+(p21)))+((p12)*((p01)+((2.0)*(p11))+(p21)))+((p22)*((p01)+(p11)+((2.0)*(p21))))))*pa[i];
//		ATA[8] += ((1.0)/(24.0)*(((p02)*(((2.0)*(p02))+(p12)+(p22)))+((p12)*((p02)+((2.0)*(p12))+(p22)))+((p22)*((p02)+(p12)+((2.0)*(p22))))))*pa[i];
//	}
//
//	// find inverse of lagrange multiplied ATA matrix
//	float C[9];
//	float D[3];
//	{
//		const float &a = ATA[0];
//		const float &b = ATA[1];
//		const float &c = ATA[2];
//		const float &d = pnorm[0];
//		//const float &e = ATA[3];
//		const float &f = ATA[4];
//		const float &g = ATA[5];
//		const float &h = pnorm[1];
//		//const float &i = ATA[6];
//		//const float &j = ATA[7];
//		const float &k = ATA[8];
//		const float &l = pnorm[2];
//
//		float det=1.0/(d*d*g*g-2*c*d*g*h+c*c*h*h-d*d*f*k+2*b*d*h*k-a*h*h*k+2*c*d*f*l-2*b*d*g*l-2*b*c*h*l+2*a*g*h*l+b*b*l*l-a*f*l*l);
//
//		C[0]=(-h*h*k+2*g*h*l-f*l*l)*det;
//		C[1]=(d*h*k-d*g*l-c*h*l+b*l*l)*det;
//		C[2]=(-d*g*h+c*h*h+d*f*l-b*h*l)*det;
//		C[3]=(d*h*k-d*g*l-c*h*l+b*l*l)*det;
//		C[4]=(-d*d*k+2*c*d*l-a*l*l)*det;
//		C[5]=(d*d*g-c*d*h-b*d*l+a*h*l)*det;
//		C[6]=(-d*g*h+c*h*h+d*f*l-b*h*l)*det;
//		C[7]=(d*d*g-c*d*h-b*d*l+a*h*l)*det;
//		C[8]=(-d*d*f+2*b*d*h-a*h*h)*det;
//
//		D[0]=(d*g*g-c*g*h-d*f*k+b*h*k+c*f*l-b*g*l)*det;
//		D[1]=(-c*d*g+c*c*h+b*d*k-a*h*k-b*c*l+a*g*l)*det;
//		D[2]=(c*d*f-b*d*g-b*c*h+a*g*h+b*b*l-a*f*l)*det;
//	}
//
//	// first point
//	for (int point = 0; point < 2; point++)
//	{
//		vect3f x;
//		float *weights;
//		if (point == 0)
//		{
//			x = *(vect3f*)op->p1 - p_cen;
//			weights = op->getWeight1();
//			op->normWeight1 = D[0]*x[0] + D[1]*x[1] + D[2]*x[2];
//		}
//		else
//		{
//			x = *(vect3f*)op->p2 - p_cen;
//			weights = op->getWeight2();
//			op->normWeight2 = D[0]*x[0] + D[1]*x[1] + D[2]*x[2];
//		}
//
//		float xC0 = x[0]*C[0] + x[1]*C[3] + x[2]*C[6];
//		float xC1 = x[0]*C[1] + x[1]*C[4] + x[2]*C[7];
//		float xC2 = x[0]*C[2] + x[1]*C[5] + x[2]*C[8];
//		float e = 1;
//
//		for (int i = 0; i < op->neighborSize; i++)
//		{
//			int neigh = op->getNeigh()[i];
//
//			int idx1 = i;
//			int idx2 = (i+1)%op->neighborSize;
//
//			vect3f p0, p1, p2;
//			vec3Copy(p0, *(vect3d*)&(faceData[neigh].he->v->pos));
//			vec3Copy(p1, *(vect3d*)&(faceData[neigh].he->next->v->pos));
//			vec3Copy(p2, *(vect3d*)&(faceData[neigh].he->next->next->v->pos));
//
//			float p00 = p0[0] - p_cen[0];
//			float p01 = p0[1] - p_cen[1];
//			float p02 = p0[2] - p_cen[2];
//			float p10 = p1[0] - p_cen[0];
//			float p11 = p1[1] - p_cen[1];
//			float p12 = p1[2] - p_cen[2];
//			float p20 = p2[0] - p_cen[0];
//			float p21 = p2[1] - p_cen[1];
//			float p22 = p2[2] - p_cen[2];
//		
//			float a0 = (p10*xC0 + p20*xC0 + p11*xC1 + p21*xC1 + p12*xC2 + p22*xC2 + 2*(p00*xC0 + p01*xC1 + p02*xC2)) * pa[i]/24;
//			float a1 = (p00*xC0 + p20*xC0 + p01*xC1 + p21*xC1 + p02*xC2 + p22*xC2 + 2*(p10*xC0 + p11*xC1 + p12*xC2)) * pa[i]/24;
//			float a2 = (p00*xC0 + p10*xC0 + p01*xC1 + p11*xC1 + p02*xC2 + p12*xC2 + 2*(p20*xC0 + p21*xC1 + p22*xC2)) * pa[i]/24;
//
//			weights[idx1] += a0;
//			weights[idx2] += a1;
//
//			e -= a0 + a1 + a2;
//		}
//
//		for (int i = 0; i < op->neighborSize; i++)
//		{
//			int idx1 = i;
//			int idx2 = (i+1)%op->neighborSize;
//
//			weights[idx1] += e*pa[i]/(3*p_area);
//			weights[idx2] += e*pa[i]/(3*p_area);
//		}
//	}
//#endif
//
//	assert ( rvalue != NULL );
	return rvalue;
}


void HalfEdgeMesh::expand()
{
	HE_HalfEdge *current, *start;
	if (currentVerts >= totalVerts)
		return;

	HE_ExpandOp *op = expandOps + currentVerts;

	//** reconnect all the pieces back together
	assert(op->he[0]->v->he->flip->v == op->he[0]->v);
	assert(op->he[1]->v == op->he[0]->v);

	// expand first wing
	op->f[0]->hole = false;

	if(op->hec[0][0]->flip != NULL)
		op->hec[0][0]->flip = op->he[0];
	
	if(op->hec[0][1]->flip != NULL)
		op->hec[0][1]->flip = op->he[0]->flip;
	
	if(op->hec[0][0]->flip != NULL)
		op->hec[0][0]->flip->flip = op->hec[0][0];
	
	if(op->hec[0][1]->flip != NULL)
		op->hec[0][1]->flip->flip = op->hec[0][1];

	op->hec[0][1]->e = op->e[0]; // this should be redundant
	
	if(op->hec[0][1]->flip != NULL)
		op->hec[0][1]->flip->e = op->e[0];

	// expand second wing
	op->f[1]->hole = false;

	if(op->hec[1][0]->flip != NULL)
		op->hec[1][0]->flip = op->he[1];

	if(op->hec[1][1]->flip != NULL)
	op->hec[1][1]->flip = op->he[1]->flip;
	op->hec[1][0]->flip->flip = op->hec[1][0];
	op->hec[1][1]->flip->flip = op->hec[1][1];
	op->hec[1][1]->e = op->e[1]; // this should be redundant
	op->hec[1][1]->flip->e = op->e[1];

	// connect half edges back to the new vertex
	current = start = op->hec[0][0]->next;
	do {
		assert(current->flip->flip == current);
		current->v = op->v;
		current = current->next->flip;
	} while (current != start);

	// connect half edges back to the old vertex (this step should be redundant)
	/*current = start = op->hec[1][1];
	do {
		assert(current->flip->flip == current);
		current->v = op->he[0]->v;
		current = current->next->flip;
	} while (current != start);*/

	op->he[0]->v->he = op->hec[0][0];

	assert(op->v->he->flip->v == op->v);
	assert(op->he[0]->v->he->flip->v == op->he[0]->v);

	// set positions
	op->he[0]->v->pos = op->pos[0];
	op->v->pos = op->pos[1];

	currentVerts++;
}

vect2d rot(vect2d me, double amt)
{
	vect2d t;
	t[0] = me[0] * cos(amt) - me[1] * sin(amt);
	t[1] = me[1] * cos(amt) + me[0] * sin(amt);
	return t;
}

#define REG_BOUND_SPACING
void HalfEdgeMesh::perform_floaters()
{
	//Setup collapsable
	for (int i = 0; i < totalVerts; i++)
	{
		if (collapsed[i] == -1)
			vertexData[i].collapsable = true;
	}
	int bsize = col_lev[col_lev.size() - 1].boundaries[chart_boundary_ind].size();
	//Get interior and boundary size
	int interiorSize = currentVerts - bsize; // current_boundary_size[chart_boundary_ind];
	printf("perform floaters HM interior %d, boundary %d\n", interiorSize, bsize);// current_boundary_size[chart_boundary_ind]);

	//Init data structures
	int row = interiorSize;
	int col = interiorSize;
	MathMatrix::SparseMatrix<double> AtA(row, col);
	AtA.zero();

	double *Xu = new double[row];
	double *Xv = new double[row];
	double *Bu = new double[row];
	double *Bv = new double[row];

	for (int i = 0; i < interiorSize; i++)
	{
		Xu[i] = 0;
		Xv[i] = 0;
		Bu[i] = 0;
		Bv[i] = 0;
	}

	//Build Iso triangles, and get area
	bool do_general = false;
	//do_general = true;
	double tarea = 0;
	vector<vect<3, vect2d>> isos;
	for (int it = 0; it < totalFaces; it++)
	{
		if (!faceData[it].hole)
		{
			//QUESTION.... SHOULD I ADD HOLE POLYS
			vect<3, vect2d> flat;

			// first vertex is straight in x, with length of the edge
			vect3d p[3];
			p[0] = faceData[it].he->v->pos;
			p[1] = faceData[it].he->next->v->pos;
			p[2] = faceData[it].he->next->next->v->pos;

			if(faceData[it].he->v->index == faceData[it].he->next->v->index || faceData[it].he->v->index == faceData[it].he->next->next->v->index || faceData[it].he->next->next->v->index == faceData[it].he->next->v->index)
				printf("wtf???\n");
			vect3d X, Y, Z;
			X = p[1] - p[0];
			X.normalize();
			Z = X % (p[2] - p[0]);
			Z.normalize();

			Y = Z % X;

			// store
			flat[0].set(0, 0);
			flat[1].set((p[1] - p[0]) * X, 0);
			flat[2].set((p[2] - p[0]) * X, (p[2] - p[0]) * Y);

			if (flat[1][0] * flat[2][1] < .00000001)
			{
			//	printf("iso shit\n");
			}
			if(!_finite(flat[1][0] * flat[2][1]))
			{
				printf("%f, %f, %f\n", faceData[it].he->v->pos[0] , faceData[it].he->v->pos[1] , faceData[it].he->v->pos[2] );
				printf("%f, %f, %f\n", faceData[it].he->next->v->pos[0] , faceData[it].he->next->v->pos[1] , faceData[it].he->next->v->pos[2] );
				printf("%f, %f, %f\n", faceData[it].he->next->next->v->pos[0] , faceData[it].he->next->next->v->pos[1] , faceData[it].he->next->next->v->pos[2] );
				printf("%d, %d, %d\n", faceData[it].he->v->index, faceData[it].he->next->v->index, faceData[it].he->next->next->v->index);
				printf("shit\n");
				do_general = true;
			}
			else
				tarea += flat[1][0] * flat[2][1];
			
			
			isos.push_back(flat);
		}
	}
	//Get boundary length
	double boundaryLength = 0;
	vector<int> curbound;
	for (int i = 0; i < boundaries[chart_boundary_ind].size(); i++)
	{
		while (! (collapsed[boundaries[chart_boundary_ind][i]->index] == -1) )
		{
			i++;
			if (i == boundaries[chart_boundary_ind].size())
			{
				break;
			}
		}
		if (i < boundaries[chart_boundary_ind].size())
		{
			curbound.push_back(boundaries[chart_boundary_ind][i]->index);
		}

	}

	//SANITY
	if (curbound.size() != current_boundary_size[chart_boundary_ind])
	{
		printf("not finding full outside boundary\n");
	}

	for (int i = 0; i < curbound.size(); i++)
		boundaryLength += (vertexData[curbound[i]].pos - vertexData[curbound[(i + 1) % curbound.size()]].pos).length();

	//Sertup Circle Boundary
	double radius = sqrt(tarea / 3.14159265359);;
	printf("area = %f\n", tarea);
	printf("bound length = %.10g\n", boundaryLength);
	vector<vect2d> newbound;
	vect2d t;
	t[1] = radius;
	t[0] = 0;
	vertexData[curbound[0]].uvpos = t;
	newbound.push_back(t);

	/*for(int i = 0; i < totalVerts; i++)
	{
		vertexData[i].uvpos[0] = 0;
		vertexData[i].uvpos[1] = 0;
	}*/
	for (int i = 1; i < curbound.size(); i++)
	{
		double l = (vertexData[curbound[i - 1]].pos - vertexData[curbound[i]].pos).length();
#ifdef REG_BOUND_SPACING
		double rotme = (boundaryLength/curbound.size() * 2 * 3.14159265359) / boundaryLength;
#else
		double rotme = (l * 2 * 3.14159265359) / boundaryLength;
#endif
		t = rot(newbound[i - 1], -1*rotme);
		vertexData[curbound[i]].uvpos = t;
		newbound.push_back(t);
	}


	vector<int> indchange;
	int vInd = 0;
	for (int i = 0; i < totalVerts; i++)
	{
//		if (!vertexData[i].collapsable ||  vertexData[i].which_boundary == chart_boundary_ind)
			if (!(collapsed[i] == -1) ||  vertexData[i].which_boundary == chart_boundary_ind)
		{
			indchange.push_back(-1);
		}
		else
		{
			indchange.push_back(vInd);
			vInd++;
		}
	}

	
	do_general = true;
	if(do_general)
	{
		/*vector<vector<int>> bverts;
		for (map<vect2i, int>::iterator it = g.charts.edge_count.begin(); it != g.charts.edge_count.end(); ++it)
		{
			vect2i t;
			t[0] = it->first.v[0];
			t[1] = it->first.v[1];

			bverts[t[0]].push_back( t[1] );
			bverts[t[1]].push_back( t[0] );
		}*/
		int dumbcount = 0;
		for (int i = 0; i < totalVerts; i++)
		{
			if(collapsed[i] != -1)
				dumbcount++;
		}
		printf("dumbcount = %d\n",dumbcount);
		vInd = 0;
		for (int i = 0; i < totalVerts; i++)
		{
			if(vertexData[i].boundary || collapsed[i] != -1)
			{
				//Does not effect A
			}
			else
			{
				//Loop through boundary
				int cbound = 0;
				HE_HalfEdge *start = vertexData[i].he;
				HE_HalfEdge *now = start;
				do
				{
					if(now->prev->v->boundary)
					{
						Bu[vInd] += vertexData[now->prev->v->index].uvpos[0];
						Bv[vInd] += vertexData[now->prev->v->index].uvpos[1];
					}
					else
					{
						AtA[vInd][indchange[now->prev->v->index]] = -1;
					}
					now = now->next->flip;
					cbound++; 
				}while(now != start);


				AtA[vInd][vInd] = cbound;
				//A(vInd, vInd) = bverts[i].size();
				//float d = 1.0f /*/  (float)bverts[i].size()*/;
				//for(int j = 0; j < bverts[i].size(); j++)
				//{
				//	if(indchange[bverts[i][j]] != -1)//Make A
				//		A(vInd, indchange[bverts[i][j]]) = -1*d;
				//	else
				//	{
				//		//Make B
				//		/*Bu(vInd, 0) += d * mesh.positions[bverts[i][j]][0];
				//		Bv(vInd, 0) += d * mesh.positions[bverts[i][j]][1];*/
				//		Bu(vInd, 0) += d * mesh.tex_coords[bverts[i][j]][0];
				//		Bv(vInd, 0) += d * mesh.tex_coords[bverts[i][j]][1];
				//	}
				//}

				vInd++;
			}
	}

	}
	else
	{
		int dumbcount = 0;
		for (int i = 0; i < totalVerts; i++)
		{
			if(collapsed[i] != -1)
				dumbcount++;
		}
		printf("dumbcountF= %d\n",dumbcount);
		for (int i = 0; i < totalFaces; i++)
		{
			if (!faceData[i].hole)
			{
				HE_HalfEdge *he = faceData[i].he;
				for (int j = 0; j < 3; j++)
				{
					HE_HalfEdge *cur = he->prev;

					if (cur->v->which_boundary == chart_boundary_ind ||
						(cur->v->which_boundary == chart_boundary_ind && he->v->which_boundary == chart_boundary_ind) ||
						(cur->v->which_boundary == chart_boundary_ind && he->next->v->which_boundary == chart_boundary_ind))
					{
						//Edge is boundary edge and do not contribute to the system
					}
					else
					{
						vect3d a = he->v->pos - cur->v->pos;
						vect3d b = he->next->v->pos - cur->v->pos;

						//sum contribution to diagonal
						double ang = (acos(b.dot(a) / (a.length()*b.length()))) / 2.0f;

						//create non diagonal entry... Need to normalize
						if (he->v->which_boundary == chart_boundary_ind)
						{
							AtA[indchange[cur->v->index]][indchange[cur->v->index]] += tan(ang) / a.length();// + .000001;
							//boundary verts contribute to B
							Bu[indchange[cur->v->index]] += ((tan(ang) / a.length()) * he->v->uvpos[0]);
							Bv[indchange[cur->v->index]] += ((tan(ang) / a.length()) * he->v->uvpos[1]);
							/*if (!isfinite(((tan(ang) / a.length()) * he->v->uvpos[1])))
								printf("a not finite\n");*/
						}
						else
						{
							AtA[indchange[cur->v->index]][indchange[he->v->index]] += (-1 * tan(ang) / a.length());
							AtA[indchange[cur->v->index]][indchange[cur->v->index]] += (tan(ang) / a.length());
						}

						if (he->next->v->which_boundary == chart_boundary_ind)
						{
							//boundary verts contribute to B
							AtA[indchange[cur->v->index]][indchange[cur->v->index]] += tan(ang) / b.length();// + .000001;
							Bu[indchange[cur->v->index]] += ((tan(ang) / b.length()) * he->next->v->uvpos[0]);
							Bv[indchange[cur->v->index]] += ((tan(ang) / b.length()) * he->next->v->uvpos[1]);
							/*if (!isfinite(((tan(ang) / b.length()) * he->next->v->uvpos[1])))
								printf("b not finite\n");*/
						}
						else
						{
							//boundary verts contribute to B
							AtA[indchange[cur->v->index]][indchange[he->next->v->index]] += (-1 * tan(ang) / b.length());
							AtA[indchange[cur->v->index]][indchange[cur->v->index]] += (tan(ang) / b.length());
						}
					}

					he = he->next;
					cur = he->next;
				}
			}
		}
	}

	//Now do for imaginary faces
	///////////////////////////////////////////
	//for (int i = 0; i < boundaries.size(); i++)
	//{
	//	if (i != chart_boundary_ind)
	//	{
	//		//May wanna make a list the simply loop through
	//		vector<int> innerbound;
	//		for (int j = 0; j < boundaries[i].size();j++)
	//		{
	//			if (boundaries[i][j]->collapsable)
	//				innerbound.push_back( boundaries[i][j]->index );// was j
	//		}

	//		//Determine ordering of fake tris
	//		bool flip = false;

	//		vect2d a = faceData[0].he->next->v->uvpos - faceData[0].he->v->uvpos;
	//		vect2d b = faceData[0].he->prev->v->uvpos - faceData[0].he->v->uvpos;

	//		float area = a[0] * b[1] - a[1] * b[0];

	//		vect3i fake;
	//		fake[0] = innerbound[0];
	//		fake[1] = innerbound[1];
	//		fake[2] = innerbound[2];

	//		a = vertexData[fake[1]].uvpos - vertexData[fake[0]].uvpos;
	//		b = vertexData[fake[2]].uvpos - vertexData[fake[0]].uvpos;

	//		float area2 = a[0] * b[1] - a[1] * b[0];
	//		if ((area < 0 && area2 < 0) || (area >= 0 && area2 >= 0))
	//			flip = false;
	//		else
	//			flip = true;

	//		//Loop through... finding groups of 3
	//		for (int j = 0; j < innerbound.size(); j++)
	//		{
	//			if (flip)
	//			{
	//				printf("\!n");
	//				fake[1] = innerbound[(j + 2) % innerbound.size()];
	//				fake[2] = innerbound[(j + 1) % innerbound.size()];
	//				fake[0] = innerbound[j];

	//				vect3d aa = vertexData[fake[2]].pos - vertexData[fake[0]].pos;
	//				vect3d bb = vertexData[fake[1]].pos - vertexData[fake[0]].pos;

	//				//sum contribution to diagonal
	//				double ang = (acos(bb.dot(aa) / (aa.length()*bb.length()))) / 2.0f;

	//				{
	//					AtA[indchange[fake[0]]][indchange[fake[2]]] += (-1 * tan(ang) / aa.length());
	//					AtA[indchange[fake[0]]][indchange[fake[0]]] += (tan(ang) / aa.length());
	//				}

	//				{
	//					AtA[indchange[fake[0]]][indchange[fake[1]]] += (-1 * tan(ang) / bb.length());
	//					AtA[indchange[fake[0]]][indchange[fake[0]]] += (tan(ang) / bb.length());
	//				}
	//			}
	//			else
	//			{
	//				fake[2] = innerbound[(j + 2) % innerbound.size()];
	//				fake[1] = innerbound[(j + 1) % innerbound.size()];
	//				fake[0] = innerbound[j];

	//				vect3d aa = vertexData[fake[1]].pos - vertexData[fake[0]].pos;
	//				vect3d bb = vertexData[fake[2]].pos - vertexData[fake[0]].pos;

	//				//sum contribution to diagonal
	//				double ang = (acos(bb.dot(aa) / (aa.length()*bb.length()))) / 2.0f;

	//				{
	//					AtA[indchange[fake[0]]][indchange[fake[1]]] += (-1 * tan(ang) / aa.length());
	//					AtA[indchange[fake[0]]][indchange[fake[0]]] += (tan(ang) / aa.length());
	//				}

	//				{
	//					AtA[indchange[fake[0]]][indchange[fake[2]]] += (-1 * tan(ang) / bb.length());
	//					AtA[indchange[fake[0]]][indchange[fake[0]]] += (tan(ang) / bb.length());
	//				}
	//			}

	//			
	//			
	//		}

	//		//int j = 0; 
	//		//while (!boundaries[i][j]->collapsable)
	//		//	j++;

	//		////find starting edge
	//		//HE_HalfEdge *he = boundaries[i][j]->he;
	//	}
	//}

	//Solve conjugant gradients
	MathMatrix::SparseMatrix<double> AtATranspose = AtA.transpose();
	LList<MathMatrix::SparseData<double> > **rowData = AtATranspose.getSparseData();

	for (int j = 0; j < AtATranspose.getNumRows(); j++)
	{
		LList<MathMatrix::SparseData<double> > *myPtr = rowData[j];

		while (myPtr != NULL)
		{
			Xu[j] += myPtr->getData().data * Bu[myPtr->getData().c];
			Xv[j] += myPtr->getData().data * Bv[myPtr->getData().c];
			myPtr = myPtr->getNext();
		}
	}

	AtA = AtATranspose * AtA;

	for (int j = 0; j < row; j++)
	{
		Bu[j] = Xu[j];
		Bv[j] = Xv[j];
		/*if (!isfinite(Bv[j]))
			printf("Xv not finite\n");*/
		Xu[j] = Xv[j] = 0;

	//	AtA[j][j] += .000001;
	}

	

	/*AtA.conjugateGradients(Bu, Xu, 200000, 0.00000001);
	AtA.conjugateGradients(Bv, Xv, 200000, 0.00000001);
*/
	AtA.conjugateGradients(Bu, Xu, 1000000, 0.00000001);
	AtA.conjugateGradients(Bv, Xv, 1000000, 0.00000001);


	//Final solution
	for (int i = 0; i < indchange.size(); i++)
	{
		if (indchange[i] != -1)
		{
		//	if (!isfinite(Xv[indchange[i]]))
		//		printf("y! not finite (%f, %f)\n", vertexData[i].uvpos[1], Xv[indchange[i]]);
			
			if (Xv[indchange[i]] > 10000)
				printf("y! not finite (%f, %f)\n", vertexData[i].uvpos[1], Xv[indchange[i]]);
			//Not boundary
			vertexData[i].uvpos[0] = Xu[indchange[i]];
			vertexData[i].uvpos[1] = Xv[indchange[i]];
			
		}
	}

	//Clean up
	delete[] Xu;
	delete[] Bu;

	delete[] Xv;
	delete[] Bv;


}
