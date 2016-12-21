#include "global.h"
#include <math.h>
#include "def_common.h"
#include "SparseMatrix.h"

//#define GENERAL
Global g;

#define MATHPI 3.14159265359
#define Radius 10
vect2d rotate(vect2d me, double amt)
{
	vect2d t;
	t[0] = me[0]*cos(amt) - me[1]*sin(amt);
	t[1] = me[1]*cos(amt) + me[0]*sin(amt);
	return t;
}

double get_angle(vect3d v1, vect3d v2)
{
	return acos( (v1.dot(v2) )/( v1.length() * v2.length() ) );
}
vect2d Global::calc_new_vert( vect3d p, vector<int> f, bool interior )
{
	vector<double> weights;
	weights.clear();

	for (int i = 0; i < f.size(); i++ )
	{
		weights.push_back(0);
	}

	vect2d newp;
	newp[0] = 0; newp[1] = 0;

	if(interior)
	{
		for(int i = 0; i < f.size(); i++)
		{
			int p1 = (i+1)%f.size();
			int m1 = (i-1 + f.size())%f.size();
			//printf("(%d, %d, %d), %d\n", m1,i,p1,f.size());

			double ap1 = get_angle( halfmesh.vertexData[f[i]].pos - halfmesh.vertexData[f[p1]].pos, p - halfmesh.vertexData[f[p1]].pos );
			double am1 = get_angle( halfmesh.vertexData[f[i]].pos - halfmesh.vertexData[f[m1]].pos, p - halfmesh.vertexData[f[m1]].pos );

			weights[i] = (atan( ap1 ) + atan( am1 )) / 2.0f;
		}

		double sumw = 0;
		for(int i = 0; i < f.size(); i++)
		{
			sumw += weights[i];
			newp += weights[i] * halfmesh.vertexData[f[i]].uvpos;
		}
		newp /= sumw;
	}
	else
	{
		double w[2][4];
		for(int i = 0; i < 4; i++)
		{
			w[0][i] = 0;
			w[1][i] = 0;
		}

		double sumw = 0;
		for(int i = 0; i < f.size()-1 /*stop at boundary*/; i++)
		{
			newp[0] += atan( get_angle( p - halfmesh.vertexData[f[i+1]].pos, halfmesh.vertexData[f[i]].pos - p - halfmesh.vertexData[f[i+1]].pos  ) )* halfmesh.vertexData[f[i+1]].uvpos[0] +
				1.0*halfmesh.vertexData[i+1].uvpos[1] +
				atan( get_angle( p - halfmesh.vertexData[f[i]].pos, halfmesh.vertexData[f[i+1]].pos - p - halfmesh.vertexData[f[i]].pos  ) ) * halfmesh.vertexData[i].uvpos[0] - 
				1.0* halfmesh.vertexData[i].uvpos[1];

			newp[1] += -1*halfmesh.vertexData[i+1].uvpos[0] +
				atan( get_angle( p - halfmesh.vertexData[f[i+1]].pos, halfmesh.vertexData[f[i]].pos - p - halfmesh.vertexData[f[i+1]].pos  ) )*halfmesh.vertexData[i+1].uvpos[1] + 
				1*halfmesh.vertexData[i].uvpos[0] + 
				atan( get_angle( p - halfmesh.vertexData[f[i]].pos, halfmesh.vertexData[f[i+1]].pos - p - halfmesh.vertexData[f[i]].pos  ) )*halfmesh.vertexData[i].uvpos[1];

			/*w[0][0] += atan( get_angle( p - halfmesh.vertexData[f[i+1]].pos, halfmesh.vertexData[f[i]].pos - p - halfmesh.vertexData[f[i+1]].pos  ) );
			w[0][1] += 1;
			w[0][2] += atan( get_angle( p - halfmesh.vertexData[f[i]].pos, halfmesh.vertexData[f[i+1]].pos - p - halfmesh.vertexData[f[i]].pos  ) );
			w[0][3] += -1;

			w[1][0] += -1;
			w[1][1] += atan( get_angle( p - halfmesh.vertexData[f[i+1]].pos, halfmesh.vertexData[f[i]].pos - p - halfmesh.vertexData[f[i+1]].pos  ) );
			w[1][2] += 1;
			w[1][3] += atan( get_angle( p - halfmesh.vertexData[f[i]].pos, halfmesh.vertexData[f[i+1]].pos - p - halfmesh.vertexData[f[i]].pos  ) );*/

			sumw += atan( get_angle( p - halfmesh.vertexData[f[i+1]].pos, halfmesh.vertexData[f[i]].pos - p - halfmesh.vertexData[f[i+1]].pos  ) ) + 
				    atan( get_angle( p - halfmesh.vertexData[f[i]].pos, halfmesh.vertexData[f[i+1]].pos - p - halfmesh.vertexData[f[i]].pos  ) );
		}

		newp /= sumw;
	}

	return newp;
}

void Global::convert_halfmesh()
{
	int ind = 0;
	for (int i = 0; i < halfmesh.totalVerts; i++)
	{
		if (halfmesh.collapsed[i] == -1)
		{
			halfmesh.vertexData[i].recalc = ind;
			ind++;
			smesh.positions.push_back( halfmesh.vertexData[i].pos );
			smesh.tex_coords.push_back( halfmesh.vertexData[i].uvpos );
		}

	}

	int i;
	int faceind = 0;
	for(i = 0; i < halfmesh.totalFaces; i++)
	{
		/*if(i == 297)
		{
			printf("\n");
		}*/
		if(!halfmesh.faceData[i].hole)
		{
			vect3i t;
			t[0] = halfmesh.faceData[i].he->v->recalc;
			t[1] = halfmesh.faceData[i].he->next->v->recalc;
			t[2] = halfmesh.faceData[i].he->next->next->v->recalc;

			smesh.indices_pos.push_back(t);
			smesh.indices_tex.push_back(t);
		}
	}

	printf("charts\n");
	g.scharts.merge(smesh);
	printf("boundary\n");
	g.get_boundary_order(&smesh, &scharts,&sbo);
	printf("iso\n");
	g.perform_isometric_flattening(&smesh, &siso_tris);
	printf("floaters\n");
	g.perform_floater_param(&smesh, &scharts, &sbo);
	
}

void Global::perform_floater_param(MeshData *m, MeshChart *c, vector<int> *bo)
{
	//Setup opt
	int interiorSize = m->tex_coords.s - bo->size();
	printf("interior %d, boundary %d\n", interiorSize, bo->size());

	int row= interiorSize;
	int col = interiorSize;
	MathMatrix::SparseMatrix<double> AtA ( row, col );
	AtA.zero();

	/*MathMatrix::Matrix<double> Xu(row, 1);
	MathMatrix::Matrix<double> Xv(row, 1);
	MathMatrix::Matrix<double> Bu(row, 1);
	MathMatrix::Matrix<double> Bv(row, 1);
*/

	double *Xu = new double[row];
	double *Xv = new double[row];
	double *Bu = new double[row];
	double *Bv = new double[row];

	for(int i =0; i < interiorSize; i++)
	{
		Xu[i] = 0;
		Xv[i] = 0;
		Bu[i] = 0;
		Bv[i] = 0;
	}

	int *boundary = new int[m->tex_coords.s];
	for(int i = 0; i < m->tex_coords.s; i++)
		boundary[i] = 0;

	vector<vector<int>> bverts;
	for(int i = 0; i < m->positions.s; i++)
	{
		vector<int> t;
		bverts.push_back(t);
	}

	double boundaryLength = 0;
	for (map<vect2i, int>::iterator it = c->edge_count.begin(); it != c->edge_count.end(); ++it)
	{
		vect2i t;
		t[0] = it->first.v[0];
		t[1] = it->first.v[1];

		if (it->second == 1)
		{
			boundary[t[0]]=1;
			boundary[t[1]]=1;

			boundaryLength += (m->positions[t[0]]-m->positions[t[1]]).length() ;
		}

		bverts[t[0]].push_back( t[1] );
		bverts[t[1]].push_back( t[0] );
	}

	vector<int> indchange;
	int vInd = 0;
	for (int i = 0; i < m->positions.s; i++)
	{
		if(boundary[i]==1)
		{
			indchange.push_back(-1);
		}
		else
		{
			indchange.push_back(vInd);
			vInd++;
		}
	}

	double tarea = 0;
	for (int it = 0; it < m->indices_pos.s; it++)
	{
		vect<3, vect2d> flat;

		// first vertex is straight in x, with length of the edge
		vect3d p[3];
		p[0] = m->positions[m->indices_pos[it][0]];
		p[1] = m->positions[m->indices_pos[it][1]];
		p[2] = m->positions[m->indices_pos[it][2]];

		vect3d X, Y, Z;
		X = p[1] - p[0];
		X.normalize();
		Z = X % (p[2] - p[0]);
		Z.normalize();

		Y = Z % X;

		// store
		flat[0].set(0, 0);
		flat[1].set((p[1]-p[0]) * X, 0);
		flat[2].set((p[2]-p[0]) * X, (p[2]-p[0]) * Y);

		if(flat[1][0] * flat[2][1] < .00000001)
		{
			printf("iso shit\n");
		}
		tarea += flat[1][0] * flat[2][1];
	}

	printf("tarea floaters = %f\n", tarea);
	double radius = sqrt(tarea / 3.14159265359);;
	vector<vect2d> newbound;
	vect2d t;
	t[0] = radius;
	t[1] = 0;

	newbound.push_back(t);
	for(int i = 1; i < bo->size(); i++)
	{
		double l = (m->positions[bo->operator[](i-1)] - m->positions[bo->operator[](i)]).length();
		double rotme = (l*2*MATHPI) / boundaryLength;

		newbound.push_back(rotate(newbound[i-1], rotme) );
	}
	for(int i = 0; i < bo->size(); i++)
	{
		m->tex_coords[bo->operator[](i)] = newbound[i];
	}

#ifdef GENERAL
	vInd = 0;
	for (int i = 0; i < m->positions.s; i++)
	{
		if(boundary[i]==1)
		{
			//Does not effect A
		}
		else
		{
			AtA[vInd][vInd] = bverts[i].size();// 1.0000f;
			//A(vInd, vInd) = bverts[i].size();
			float d = 1.0f /*/  (float)bverts[i].size()*/;
			for(int j = 0; j < bverts[i].size(); j++)
			{
				if(indchange[bverts[i][j]] != -1)//Make A
					AtA[vInd][indchange[bverts[i][j]]] = -1*d;
				else
				{
					//Make B
					/*Bu(vInd, 0) += d * m->positions[bverts[i][j]][0];
					Bv(vInd, 0) += d * m->positions[bverts[i][j]][1];*/
					Bu[vInd] += d * m->tex_coords[bverts[i][j]][0];
					Bv[vInd] += d * m->tex_coords[bverts[i][j]][1];
				}
			}

			vInd++;
		}
	}
#else
	for(int i = 0; i < m->indices_pos.s; i++)
	{
		for(int j =0; j < 3; j++)
		{
			int vnn = (j+1)%3;
			int vn = (j+2)%3;

			if  (boundary[m->indices_pos[i][j]] == 1 || 
				(boundary[m->indices_pos[i][j]]==1 && boundary[m->indices_pos[i][vn]]==1) || 
				(boundary[m->indices_pos[i][j]]==1 && boundary[m->indices_pos[i][vnn]]==1))
			{
				//Edges are boundary edges and do not contribute to the system
			}
			else
			{
				//Get angle between vn and vnn
				vect3d a =m->positions[m->indices_pos[i][vn]]- m->positions[m->indices_pos[i][j]];
				vect3d b =m->positions[m->indices_pos[i][vnn]]- m->positions[m->indices_pos[i][j]];

				//sum contribution to diagonal
				double ang = (acos(b.dot(a) / (a.length()*b.length())))/2.0f;

				//create non diagonal entry... Need to normalize
				if(boundary[m->indices_pos[i][vn]]==1)
				{
					//printf("%d, %d, %f\n",indchange[m->indices_pos[i][j]],indchange[m->indices_pos[i][j]], AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][j]]])
					AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][j]]] += tan(ang) / a.length();// + .000001;
					//boundary verts contribute to B
					Bu[indchange[m->indices_pos[i][j]]] += ((tan(ang) / a.length()) * m->tex_coords[m->indices_pos[i][vn]][0]);
					Bv[indchange[m->indices_pos[i][j]]] += ((tan(ang) / a.length()) * m->tex_coords[m->indices_pos[i][vn]][1]);
				}
				else
				{
					AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][vn]]] += (-1*tan(ang) / a.length() );
					AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][j]]] += (tan(ang) / a.length());
				}

				if(boundary[m->indices_pos[i][vnn]] == 1)
				{
					//boundary verts contribute to B
					AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][j]]] +=  tan(ang) / b.length();// + .000001;
					Bu[indchange[m->indices_pos[i][j]]] += ((tan(ang) / b.length()) * m->tex_coords[m->indices_pos[i][vnn]][0]);
					Bv[indchange[m->indices_pos[i][j]]] += ((tan(ang) / b.length()) * m->tex_coords[m->indices_pos[i][vnn]][1]);
				}
				else
				{
					//boundary verts contribute to B
					AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][vnn]]] += (-1*tan(ang) / b.length());
					AtA[indchange[m->indices_pos[i][j]]][indchange[m->indices_pos[i][j]]] += (tan(ang) / b.length());
				}
			}
		}
	}
#endif
//	AtA.outputMathematica ( "ata.txt" );
//	MathMatrix::Matrix<double> But(row, 1);
//	MathMatrix::Matrix<double> Bvt(row, 1);

	//for(int i = 0; i < row; i++)
	//{
	//	But[i][0] = Bu[i];
	//	Bvt[i][0] = Bv[i];
	//}
//	But.outputMathematica("bu.txt" );
//	Bvt.outputMathematica("bv.txt" );
	
//	MathMatrix::Matrix<double> But(row, 1);
//    MathMatrix::Matrix<double> Bvt(row, 1);
// 
//    for(int i = 0; i < row; i++)
//    {
//            But[i][0] = Bu[i];
//            Bvt[i][0] = Bv[i];
//    }
//    MathMatrix::SparseMatrix<double> AtATranspose = AtA.transpose ( );
//      
//    But =AtATranspose * But;
//    Bvt = AtATranspose * Bvt;
//    AtA = AtATranspose * AtA;
//    for(int i = 0; i < row; i++)
//    {
//            Bu[i] = But[i][0];
//            Bv[i] = Bvt[i][0];
//    }
//      
////     But.outputMathematica ( "but.txt" );
////     Bvt.outputMathematica ( "bvt.txt" );
// 
//    AtA.conjugateGradients ( Bu, Xu, 100000, 0.0000001 );
//    AtA.conjugateGradients ( Bv, Xv, 100000, 0.0000001 );

	

	MathMatrix::SparseMatrix<double> AtATranspose = AtA.transpose ( );
       LList<MathMatrix::SparseData<double> > **rowData = AtATranspose.getSparseData ( );
 
       for ( int j = 0; j < AtATranspose.getNumRows ( ); j++ )
       {
              LList<MathMatrix::SparseData<double> > *myPtr = rowData [ j ];
 
              while ( myPtr != NULL )
              {
                     Xu [ j ] += myPtr->getData ( ).data * Bu [ myPtr->getData ( ).c ];
                     Xv [ j ] += myPtr->getData ( ).data * Bv [ myPtr->getData ( ).c ];
                     myPtr = myPtr->getNext ( );
              }
       }
       for ( int j = 0; j < row; j++ )
       {
              Bu [ j ] = Xu [ j ];
              Bv [ j ] = Xv [ j ];
              Xu [ j ] = Xv [ j ] = 0;
       }
 
       AtA = AtATranspose * AtA;
 
       AtA.conjugateGradients ( Bu, Xu, 1000000, 0.0000000001 );
	   AtA.conjugateGradients ( Bv, Xv, 1000000, 0.0000000001 );


	for(int i = 0; i < indchange.size(); i++)
	{
		if(indchange[i] != -1)
		{
			//Not boundary
			m->tex_coords[i][0] = Xu[indchange[i]];
			m->tex_coords[i][1] = Xv[indchange[i]];
		}
	}

	delete[] Xu;
	delete[] Bu;

	delete[] Xv;
	delete[] Bv;
	delete[] boundary;


// 
//AtA [ i ] [ j ] += 20;
// 
//AtA.conjugateGradients ( AtB, x, 100000, 0.000001 );
// 
//AtB is a MathMatrix::Matrix<double> (rows, 1);
//x is a MathMatrix::Matrix<double> (rows, 1);
//x holds an initial guess as to the solution.  You could set the boundary vertices to their constrained positions and everything else to zero for example.
//The next number is the maximum number of iterations in the conjugate gradient algorithm (set to rows or rows*2).
//The last number is the error tolerance.  0.000001 is good.
// 



//
//	double radius = sqrt(tarea / 3.14159265359);
//	//Set new boundary tex coordinates
//	vector<vect2d> newbound;
//	vect2d t;
//#ifdef USE_RAD_MEANVALUE
//	t[0] = radius;
//#else
//	t[0] = 1;
//#endif
//	t[1] = 0;
//	newbound.push_back(t);
//	for(int i = 1; i < boundaryorder.size(); i++)
//	{
//#ifdef LOAD_OBJ_TYPE2
//		double l = (mesh.positions[boundaryorder[i-1]] - mesh.positions[boundaryorder[i]]).length();
//#else
//		double l = (mesh.tex_coords[boundaryorder[i-1]] - mesh.tex_coords[boundaryorder[i]]).length();
//#endif
//		double rotme = (l*2*MATHPI) / boundaryLength;
//
//		newbound.push_back(rotate(newbound[i-1], rotme) );
//	}
//
//	for(int i = 0; i < boundaryorder.size(); i++)
//	{
//		mesh.tex_coords[boundaryorder[i]] = newbound[i]*SCALE_OVERALL;
//	}
//
//	//Setup indices for interior vertices
//	vector<int> indchange;
//	int vInd = 0;
//#ifdef LOAD_OBJ_TYPE2
//	for (int i = 0; i < mesh.positions.s; i++)
//#else
//	for (int i = 0; i < mesh.tex_coords.s; i++)
//#endif
//	{
//		if(boundary(i,0)==1)
//		{
//			indchange.push_back(-1);
//		}
//		else
//		{
//			indchange.push_back(vInd);
//			vInd++;
//		}
//	}
//
//	//Setup matrices
//	//A, Bu, Bv, U, V
//	vector<Triplet<double>> trips;
//	trips.resize(mesh.indices_pos.s * 3);
//
//	SparseMatrix<double> A;
//	A.resize(interiorSize, interiorSize);
//	//DenseTaucsMatrix Bu, Bv, U, V;
//
//	VectorXd Bu, Bv, U, V;
//	Bu.resize(interiorSize);
//	Bv.resize(interiorSize);
//
//	for(int i =0; i < interiorSize; i++)
//	{
//		Bu(i) = 0;
//		Bv(i) = 0;
//	}
//
//	double *sum;
//	sum  = new double[mesh.indices_pos.s];
//
//	for(int i = 0; i < mesh.indices_pos.s; i++ )
//	{
//		sum[i] = 1;
//	}
//#ifdef LOAD_OBJ_TYPE2
//	printf("Build Trips\n");
//	for(int i = 0; i < mesh.indices_pos.s; i++)
//	{
//		for(int j =0; j < 3; j++)
//		{
//			int vnn = (j+1)%3;
//			int vn = (j+2)%3;
///*
//			int vn = (j+1)%3;
//			int vnn = (j+2)%3;
//*/
//			if(boundary(mesh.indices_pos[i][j],0) == 1 || 
//				(boundary(mesh.indices_pos[i][j],0)==1 && boundary(mesh.indices_pos[i][vn],0)==1) || 
//				(boundary(mesh.indices_pos[i][j],0)==1 && boundary(mesh.indices_pos[i][vnn],0)==1))
//			{
//				//Edges are boundary edges and do not contribute to the system
//			}
//			else
//			{
//				//Get angle between vn and vnn
//				vect3d a =mesh.positions[mesh.indices_pos[i][vn]]- mesh.positions[mesh.indices_pos[i][j]];
//				vect3d b =mesh.positions[mesh.indices_pos[i][vnn]]- mesh.positions[mesh.indices_pos[i][j]];
//
//			/*	vect3d a =mesh.positions[mesh.indices_pos[i][vnn]]- mesh.positions[mesh.indices_pos[i][j]];
//				vect3d b =mesh.positions[mesh.indices_pos[i][vn]]- mesh.positions[mesh.indices_pos[i][j]];*/
//
//				//sum contribution to diagonal
//				double ang = (acos(b.dot(a) / (a.length()*b.length())))/2.0f;
//
//				//create non diagonal entry... Need to normalize
//				if(boundary(mesh.indices_pos[i][vn],0)==1)
//				{
//					//boundary verts contribute to B
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]],indchange[mesh.indices_pos[i][j]],
//						tan(ang) / a.length() + .000001));
//
//					Bu(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / a.length()) * mesh.tex_coords[mesh.indices_pos[i][vn]][0])/sum[mesh.indices_pos[i][j]];
//					Bv(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / a.length()) * mesh.tex_coords[mesh.indices_pos[i][vn]][1])/sum[mesh.indices_pos[i][j]];
//				}
//				else
//				{
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][vn]],
//						(-1*tan(ang) / a.length())/sum[mesh.indices_pos[i][j]]  ) );
//
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]],indchange[mesh.indices_pos[i][j]],
//						(tan(ang) / a.length())/sum[mesh.indices_pos[i][j] ]+ .000001));
//				}
//
//				if(boundary(mesh.indices_pos[i][vnn],0)==1)
//				{
//					//boundary verts contribute to B
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][j]],
//						tan(ang) / b.length()+ .000001));
//					Bu(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / b.length()) * mesh.tex_coords[mesh.indices_pos[i][vnn]][0])/sum[mesh.indices_pos[i][j]];
//					Bv(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / b.length()) * mesh.tex_coords[mesh.indices_pos[i][vnn]][1])/sum[mesh.indices_pos[i][j]];
//				}
//				else
//				{
//					//boundary verts contribute to B
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][vnn]],
//						(-1*tan(ang) / b.length())/sum[mesh.indices_pos[i][j]]));
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][j]],
//						(tan(ang) / b.length())/sum[mesh.indices_pos[i][j]]+ .000001));
//				}
//			}
//		}
//	}
//#else
//	printf("Build Trips\n");
//	for(int i = 0; i < mesh.indices_pos.s; i++)
//	{
//		for(int j =0; j < 3; j++)
//		{
//			int vnn = (j+1)%3;
//			int vn = (j+2)%3;
///*
//			int vn = (j+1)%3;
//			int vnn = (j+2)%3;
//*/
//			if(boundary(mesh.indices_tex[i][j],0) == 1 || 
//				(boundary(mesh.indices_tex[i][j],0)==1 && boundary(mesh.indices_tex[i][vn],0)==1) || 
//				(boundary(mesh.indices_tex[i][j],0)==1 && boundary(mesh.indices_tex[i][vnn],0)==1))
//			{
//				//Edges are boundary edges and do not contribute to the system
//			}
//			else
//			{
//				//Get angle between vn and vnn
//				vect3d a =mesh.positions[mesh.indices_pos[i][vn]]- mesh.positions[mesh.indices_pos[i][j]];
//				vect3d b =mesh.positions[mesh.indices_pos[i][vnn]]- mesh.positions[mesh.indices_pos[i][j]];
//
//			/*	vect3d a =mesh.positions[mesh.indices_pos[i][vnn]]- mesh.positions[mesh.indices_pos[i][j]];
//				vect3d b =mesh.positions[mesh.indices_pos[i][vn]]- mesh.positions[mesh.indices_pos[i][j]];*/
//
//				//sum contribution to diagonal
//				double ang = (acos(b.dot(a) / (a.length()*b.length())))/2.0f;
//
//				//create non diagonal entry... Need to normalize
//				if(boundary(mesh.indices_pos[i][vn],0)==1)
//				{
//					//boundary verts contribute to B
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]],indchange[mesh.indices_pos[i][j]],
//						tan(ang) / a.length() + .000001));
//
//					Bu(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / a.length()) * mesh.tex_coords[mesh.indices_pos[i][vn]][0])/sum[mesh.indices_pos[i][j]];
//					Bv(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / a.length()) * mesh.tex_coords[mesh.indices_pos[i][vn]][1])/sum[mesh.indices_pos[i][j]];
//				}
//				else
//				{
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][vn]],
//						(-1*tan(ang) / a.length())/sum[mesh.indices_pos[i][j]]  ) );
//
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]],indchange[mesh.indices_pos[i][j]],
//						(tan(ang) / a.length())/sum[mesh.indices_pos[i][j] ]+ .000001));
//				}
//
//				if(boundary(mesh.indices_pos[i][vnn],0)==1)
//				{
//					//boundary verts contribute to B
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][j]],
//						tan(ang) / b.length()+ .000001));
//					Bu(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / b.length()) * mesh.tex_coords[mesh.indices_pos[i][vnn]][0])/sum[mesh.indices_pos[i][j]];
//					Bv(indchange[mesh.indices_pos[i][j]]) += ((tan(ang) / b.length()) * mesh.tex_coords[mesh.indices_pos[i][vnn]][1])/sum[mesh.indices_pos[i][j]];
//				}
//				else
//				{
//					//boundary verts contribute to B
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][vnn]],
//						(-1*tan(ang) / b.length())/sum[mesh.indices_pos[i][j]]));
//					trips.push_back(Triplet<double>(indchange[mesh.indices_pos[i][j]], indchange[mesh.indices_pos[i][j]],
//						(tan(ang) / b.length())/sum[mesh.indices_pos[i][j]]+ .000001));
//				}
//			}
//		}
//	}
//#endif
//
//	/*for(int i = 0; i < trips.size(); i++)
//	{
//		trips[i]. += .000001;
//	}*/
//	A.setFromTriplets(trips.begin(), trips.end());
//
//	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solverA;
//
//	printf("SOLVE\n");
//	solverA.compute(A);
//
//	printf("SOLVE U\n");
//	U = solverA.solve(Bu);
//
//	printf("SOLVE V\n");
//	V = solverA.solve(Bv);
//
//	for(int i = 0; i < indchange.size(); i++)
//	{
//		if(indchange[i] != -1)
//		{
//			//Not boundary
//			mesh.tex_coords[i][0] = U(indchange[i],0);
//			mesh.tex_coords[i][1] = V(indchange[i],0);
//		}
//	}

}

Global::Global()
{
	debug_seamless = false;;
	which_seamless = 0;
	draw_full_mesh = false;
	remove_tri_index = 0;
	ringJ = ringI = 0;
	iteration_count = 0;
	total_timing = timing = 0;

	init();

	init_draw();

	init_type = 0;
	max_param = 1;

	fully_expanded = false;

	seamless = 0;
}

void Global::init_draw()
{
	int i,j;

	// set rotation to identity
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
			rotation[i*4 + j] = (i == j ? 1 : 0);
	}
	//View setup
	focus(.5,.5,0);
	zoom =300;

	//Drawing Setup
	bg_color[0] = 1;
	bg_color[1] = 1;
	bg_color[2] = 1;
}


void j_function_grad(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2, int ind)
{
	//printf("grad");
	if(ind == 0)
	{
		printf("{%.10g, %.10g}\n",
			(2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*(((2.0)*(u0))+(-(u1))+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u1)+(-(u2))))+(((v0)+(-(v2)))*((v1)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+(((L)*(L)*(L))*((-((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*((v0)+(-(v2)))*((v1)+(-(v2))))))+(((u1)*(u1))*((v0)+(-(v2))))+(-(u0)*((u1)+(-(u2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u1)*(u2)*((-(v1))+(v2))))*(xx)*((y)*(y)))+(-((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((-((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((xx)*(xx))*((y)*(y)))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((y)*(y)*(y)*(y)))))))
			,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((2.0)*((u1)+(-(u2)))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((u1)+(-(u2)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((2.0)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((v0)+(-(v2))))+((L)*(((-2.0)*(v0))+(v1)+(v2))*(xx))+(((v0)+(-(v1)))*(((xx)*(xx))+((y)*(y)))))))));
	}
	else if(ind == 1)
	{
		printf("{%.10g, %.10g}\n",
		(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+((-2.0)*((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+((2.0)*((L)*(L)*(L))*((((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*(((v0)+(-(v2)))*((v0)+(-(v2)))))))+((u1)*(u2)*((v0)+(-(v2))))+(((u0)*(u0))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u0)*(((u1)*((-(v0))+(v2)))+((u2)*(((-3.0)*(v0))+((2.0)*(v1))+(v2))))))*(xx)*((y)*(y)))+((2.0)*((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+((-2.0)*((L)*(L))*((v0)+(-(v1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))))
		,(2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*((v0)+(-(v2)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((u0)+(-(u2)))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L)*(L))*(((-2.0)*((u0)*(u0)*(u0)))+((2.0)*((u0)*(u0))*((u1)+((2.0)*(u2))))+((u1)*(((2.0)*((u2)*(u2)))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+(-(u0)*(((2.0)*(u2)*(((2.0)*(u1))+(u2)))+(((v0)+(-(v2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+((u2)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+(((v0)+(-(v1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((u0)+(-(u1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y)))))));
	}
	else
	{
		printf("{%.10g, %.10g}\n",(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((L)*((u0)+(-(u2))))+(((-(u0))+(u1))*(xx)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y)))))+((2.0)*((-(v0))+(v1))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((-(v0))+(v1))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))))
	, (2.0/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*((((v0)+(-(v1)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(-((L)*(L)*(L))*((u0)+(-(u2)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L))*(((2.0)*((u0)*(u0)*(u0)))+((-2.0)*((u1)*(u1))*(u2))+((-2.0)*((u0)*(u0))*(((2.0)*(u1))+(u2)))+(-(u2)*(((v0)+(-(v1)))*((v0)+(-(v1)))))+((u0)*(((2.0)*(u1)*((u1)+((2.0)*(u2))))+(((v0)+(-(v1)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+(-(u1)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+((L)*((((v0)+(-(v2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((xx)*(xx))*((y)*(y)))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((y)*(y)*(y)*(y))))))));
	}
	/*printf("grad = {%.10g, %.10g, %.10g, %.10g,%.10g, %.10g}\n", (2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*(((2.0)*(u0))+(-(u1))+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u1)+(-(u2))))+(((v0)+(-(v2)))*((v1)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+(((L)*(L)*(L))*((-((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*((v0)+(-(v2)))*((v1)+(-(v2))))))+(((u1)*(u1))*((v0)+(-(v2))))+(-(u0)*((u1)+(-(u2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u1)*(u2)*((-(v1))+(v2))))*(xx)*((y)*(y)))+(-((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((-((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((xx)*(xx))*((y)*(y)))+(((-(v0))+(v1))*((-((u0)+(-(u1)))*((u1)+(-(u2))))+(-((v0)+(-(v1)))*((v1)+(-(v2)))))*((y)*(y)*(y)*(y)))))))
	,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((2.0)*((u1)+(-(u2)))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((u1)+(-(u2)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((2.0)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((v0)+(-(v2))))+((L)*(((-2.0)*(v0))+(v1)+(v2))*(xx))+(((v0)+(-(v1)))*(((xx)*(xx))+((y)*(y))))))))
	,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*((u0)+(-(u2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(xx))+((-2.0)*((L)*(L)*(L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((v0)+(-(v2)))*((y)*(y)))+((2.0)*((L)*(L)*(L))*((((v0)+(-(v1)))*(((u2)*(u2))+((2.0)*(((v0)+(-(v2)))*((v0)+(-(v2)))))))+((u1)*(u2)*((v0)+(-(v2))))+(((u0)*(u0))*(((2.0)*(v0))+(-(v1))+(-(v2))))+((u0)*(((u1)*((-(v0))+(v2)))+((u2)*(((-3.0)*(v0))+((2.0)*(v1))+(v2))))))*(xx)*((y)*(y)))+((2.0)*((u0)+(-(u1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+((-2.0)*((L)*(L))*((v0)+(-(v1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))))
	,(2.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((L)*((v0)+(-(v2)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(((L)*(L)*(L)*(L))*((u0)+(-(u2)))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L)*(L))*(((-2.0)*((u0)*(u0)*(u0)))+((2.0)*((u0)*(u0))*((u1)+((2.0)*(u2))))+((u1)*(((2.0)*((u2)*(u2)))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+(-(u0)*(((2.0)*(u2)*(((2.0)*(u1))+(u2)))+(((v0)+(-(v2)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+((u2)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+(((v0)+(-(v1)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2)))))*(((xx)*(xx))+((y)*(y))))+(((L)*(L))*((u0)+(-(u1)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y))*(((xx)*(xx))+((y)*(y))))))
	,(1.0/(((L)))/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*(((-2.0)*(L)*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((L)*((u0)+(-(u2))))+(((-(u0))+(u1))*(xx)))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y)))))+((2.0)*((-(v0))+(v1))*((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))+((-2.0)*((-(v0))+(v1))*(((((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2)))*(((u1)*(v0))+(-(u2)*(v0))+(-(u0)*(v1))+((u2)*(v1))+((u0)*(v2))+(-(u1)*(v2))))+(((L)*(L))*((y)*(y))))*((((L)*(L))*((((u0)+(-(u2)))*((u0)+(-(u2))))+(((v0)+(-(v2)))*((v0)+(-(v2))))))+((-2.0)*(L)*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx))+(((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*(((xx)*(xx))+((y)*(y))))))))
	, (2.0/(((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))))/(((y)))*((((v0)+(-(v1)))*((((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2))))*(((u2)*((-(v0))+(v1)))+((u1)*((v0)+(-(v2))))+((u0)*((-(v1))+(v2)))))*(xx))+(-((L)*(L)*(L))*((u0)+(-(u2)))*(((u0)*(u0))+((u1)*(u2))+(-(u0)*((u1)+(u2)))+(((v0)+(-(v1)))*((v0)+(-(v2)))))*((y)*(y)))+(((L)*(L))*(((2.0)*((u0)*(u0)*(u0)))+((-2.0)*((u1)*(u1))*(u2))+((-2.0)*((u0)*(u0))*(((2.0)*(u1))+(u2)))+(-(u2)*(((v0)+(-(v1)))*((v0)+(-(v1)))))+((u0)*(((2.0)*(u1)*((u1)+((2.0)*(u2))))+(((v0)+(-(v1)))*(((2.0)*(v0))+(-(v1))+(-(v2))))))+(-(u1)*((v0)+(-(v1)))*((v0)+(-(v2)))))*(xx)*((y)*(y)))+((L)*((((v0)+(-(v2)))*((((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))*(((u2)*((v0)+(-(v1))))+((u0)*((v1)+(-(v2))))+((u1)*((-(v0))+(v2))))))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((xx)*(xx))*((y)*(y)))+(((-(u0))+(u1))*((((u0)+(-(u1)))*((u0)+(-(u1))))+(((v0)+(-(v1)))*((v0)+(-(v1)))))*((y)*(y)*(y)*(y)))))))
	);*/
}
void Global::print_debug_info(int ind)
{

	HE_HalfEdge *start =g.halfmesh.vertexData [ ind ].he, *curr;
		bool isbound = g.halfmesh.vertexData[ ind ].boundary;

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
				printf("GradJ[{{%.10g, %.10g}\n, {%.10g, %.10g}\n, {%.10g, %.10g}},\n %.10g, %.10g, %.10g]\n",
					g.halfmesh.vertexData[t[0]].uvpos[0],g.halfmesh.vertexData[t[0]].uvpos[1],
					g.halfmesh.vertexData[t[1]].uvpos[0],g.halfmesh.vertexData[t[1]].uvpos[1],
					g.halfmesh.vertexData[t[2]].uvpos[0],g.halfmesh.vertexData[t[2]].uvpos[1],
					g.iso_tris[faceind][1][0],g.iso_tris[faceind][2][0],g.iso_tris[faceind][2][1]
					
					);

				j_function_grad(g.iso_tris[faceind][1][0],g.iso_tris[faceind][2][0],g.iso_tris[faceind][2][1],
					g.halfmesh.vertexData[t[0]].uvpos[0],g.halfmesh.vertexData[t[0]].uvpos[1],
					g.halfmesh.vertexData[t[1]].uvpos[0],g.halfmesh.vertexData[t[1]].uvpos[1],
					g.halfmesh.vertexData[t[2]].uvpos[0],g.halfmesh.vertexData[t[2]].uvpos[1],
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
				int faceind = curr->f->index;


				g.iso_tris[faceind][1][0];

				printf("GradJ[{{%.10g, %.10g}\n, {%.10g, %.10g}\n, {%.10g, %.10g}},\n %.10g, %.10g, %.10g]\n",
					g.halfmesh.vertexData[t[0]].uvpos[0],g.halfmesh.vertexData[t[0]].uvpos[1],
					g.halfmesh.vertexData[t[1]].uvpos[0],g.halfmesh.vertexData[t[1]].uvpos[1],
					g.halfmesh.vertexData[t[2]].uvpos[0],g.halfmesh.vertexData[t[2]].uvpos[1],
					g.iso_tris[faceind][1][0],g.iso_tris[faceind][2][0],g.iso_tris[faceind][2][1]
					
					);

					j_function_grad(g.iso_tris[faceind][1][0],g.iso_tris[faceind][2][0],g.iso_tris[faceind][2][1],
					g.halfmesh.vertexData[t[0]].uvpos[0],g.halfmesh.vertexData[t[0]].uvpos[1],
					g.halfmesh.vertexData[t[1]].uvpos[0],g.halfmesh.vertexData[t[1]].uvpos[1],
					g.halfmesh.vertexData[t[2]].uvpos[0],g.halfmesh.vertexData[t[2]].uvpos[1],
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

}
void Global::init()
{
	draw_simplified = false;
	mleft = false;
	mright = false;
	mmiddle = false;
	printInfo = false;
	draw_halfmesh = false;

	gradDraw = new double[mesh.tex_coords.s*2];
	globalboundarygrad = new double[mesh.tex_coords.s*2];
	interiorgrad = new double[mesh.tex_coords.s*2];

	for(int i = 0; i < mesh.tex_coords.s*2; i++)
	{
		globalboundarygrad[i] = 0;
		interiorgrad[i] = 0;
		gradDraw[i] = 0;
	}

	srand(101);
	for (int i = 0; i < 10; i++)
		rand();

	draw_uv01 = true;
	draw_grad = false;
	draw_both_grads = false;
	display_error = false;
}

void Global::perform_isometric_flattening(MeshData *m,  vector<vector<vect2d>> *isot)
{
	isot[0].clear();
	printf("Isometric flattening\n");
	vect3d X, Y, Z;
	double tarea = 0;// flat[1][0] * flat[2][1];
	for (int it = 0; it < m->indices_pos.s; it++)
	{
		vector<vect2d> flat;
		vect2d t1,t2,t3;

		flat.push_back(t1);
		flat.push_back(t2);
		flat.push_back(t3);
		// first vertex is straight in x, with length of the edge
		vect3d p[3];
		p[0] = m->positions[m->indices_pos[it][0]];
		p[1] = m->positions[m->indices_pos[it][1]];
		p[2] = m->positions[m->indices_pos[it][2]];

		X = p[1] - p[0];
		X.normalize();
		Z = X % (p[2] - p[0]);
		Z.normalize();

		Y = Z % X;

		//// scale triangle by its initial scale (to get approximately the same size mesh as the input)
		double amb = 1;
		X *= amb;
		Y *= amb;
		
		// store
		flat[0].set(0, 0);
		flat[1].set((p[1]-p[0]) * X, 0);
		flat[2].set((p[2]-p[0]) * X, (p[2]-p[0]) * Y);

		 isot->push_back(flat);
		//area of tris
		tarea += flat[1][0] * flat[2][1];
	}
	printf("Iso area %.10g\n", tarea);
}

void Global::get_boundary_order(MeshData *m, MeshChart *mc,  vector<int> *bo)
{
	vector<vector<int>> bverts;
	vector<vect2i> boundaryedges;

	for(int i = 0; i < m->positions.s; i++)
	{
		vector<int> t;
		bverts.push_back(t);
	}

	for (map<vect2i, int>::iterator it =mc->edge_count.begin(); it != mc->edge_count.end(); ++it)
	{
		vect2i t;
		t[0] = it->first.v[0];
		t[1] = it->first.v[1];
		if (it->second == 1)//Is boundary
		{
			boundaryedges.push_back(t);
		}

		bverts[t[0]].push_back( t[1] );
		bverts[t[1]].push_back( t[0] );
	}

	bo->clear();

	int e = boundaryedges[0][1];
	int s = boundaryedges[0][0];
	int st = boundaryedges[0][0];

	bo->push_back(s);
	while (e != s)
	{
		bo->push_back(e);
		vector<vect2i> t;
		for(int i = 0; i < boundaryedges.size(); i++)
		{
			if(boundaryedges[i][0] == e || boundaryedges[i][1] == e)
			{
				t.push_back(boundaryedges[i]);
			}
		}

		for(int i = 0; i < t.size(); i++)
		{
			if(t[i][0] == e && t[i][1] != st)
			{
				e = t[i][1];
				st = t[i][0];
				i = t.size();
			}
			else if(t[i][1] == e && t[i][0] != st)
			{
				e = t[i][0];
				st = t[i][1];
				i = t.size();
			}
			else
			{
				//Prev edge do nothing
			}
		}
	}
}

void Global::center_mesh()
{
	vect3d cen;
	cen[0] = 0; cen[1] = 0; cen[2] = 0;

	vect2d cen2;
	cen2[0] = 0; cen2[1] = 0;

	for(int i = 0; i < mesh.positions.s; i++)
		cen += mesh.positions[i];

	for(int i = 0; i < mesh.tex_coords.s; i++)
		cen2 += mesh.tex_coords[i];

	vect3d cent;
	cent[0] = .5; cent[1] = .5; cent[2] = .5;
	vect3d diff = cent - cen / mesh.positions.s;

	vect2d cen2t;
	cen2t[0] = .5; cen2t[1] = .5;
	vect2d diff2 = cen2t - cen2 / mesh.tex_coords.s;

	for(int i = 0; i < mesh.positions.s; i++)
		mesh.positions[i] = mesh.positions[i] + diff;

	for(int i = 0; i < mesh.tex_coords.s; i++)
		mesh.tex_coords[i] = mesh.tex_coords[i] + diff2;
}

void Global::normalize_mesh()
{
	vect3d minS, maxS;
	minS = mesh.positions[0];
	maxS = mesh.positions[0];

	vect2d min2S, max2S;
	min2S = mesh.tex_coords[0];
	max2S = mesh.tex_coords[0];

	for(int i = 1; i < mesh.positions.s; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(mesh.positions[i][j] < minS[j])
				minS[j] = mesh.positions[i][j];

			if(mesh.positions[i][j] > maxS[j])
				maxS[j] = mesh.positions[i][j];
		}
	}

	for(int i = 1; i < mesh.tex_coords.s; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			if(mesh.tex_coords[i][j] < min2S[j])
				min2S[j] = mesh.tex_coords[i][j];
			
			if(mesh.tex_coords[i][j] > max2S[j])
				max2S[j] = mesh.tex_coords[i][j];
		}
	}

	double scale = -999999999;
	for(int i = 0; i < 3; i++)
	{
		if(maxS[i] - minS[i] > scale)
			scale = maxS[i] - minS[i];
	}

	vect3d cent;
	cent[0] = .5; cent[1] = .5; cent[2] = .5;
	for(int i = 0; i < mesh.positions.s; i++)
	{
		mesh.positions[i] = mesh.positions[i]/scale + cent;
	}

	for(int i = 0; i < 2; i++)
	{
		if(max2S[i] - min2S[i] > scale)
			scale = max2S[i] - min2S[i];
	}
	vect2d cent2;
	cent2[0] = .5; cent2[1] = .5;
	for(int i = 0; i < mesh.tex_coords.s; i++)
	{
		mesh.tex_coords[i] = mesh.tex_coords[i]/scale + cent2;
	}
}

void Global::flip_uv( MeshData *m )
{
/*
	for(int i = 0; i < mesh.indices_pos.size(); i++)
	{
		int t =	mesh.indices_pos[i][1];
		mesh.indices_pos[i][1] = mesh.indices_pos[i][2];
		 mesh.indices_pos[i][2] = t;

		t = mesh.indices_tex[i][1];
		mesh.indices_tex[i][1] = mesh.indices_tex[i][2];
		mesh.indices_tex[i][2] = t;
	}
	*/
	for(int i = 0; i < m->tex_coords.s; i++)
	{
		m->tex_coords[i][0] *= -1;
	}

}

double Global::find_ang(vect2d v, vect2d w)
{
	vect3d vv;
	vv[0] = v[0]; vv[1] = v[1]; vv[2] = 0;

	vect3d ww;
	ww[0] = w[0]; ww[1] = w[1]; ww[2] = 0;

	return vv.cross(ww)[2] >= 0 ?  acos( v.dot(w) / (v.length()*w.length() ) ) : 2*MATH_PI - acos( v.dot(w) / (v.length()*w.length() ) );
}
double Global::tri_error(double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2)
{
	double area = u1*v0 - u2*v0 - u0*v1 + u2*v1 + u0*v2 - u1*v2;
	double u01 = u0 - u1;
	double v01 = v0 - v1;
	double vecU = L*(u0 - u2) - u01*xx;
	double vecV = L*(v0 - v2) - v01*xx;
	return ( area * area + L * L * y * y ) * ( vecU * vecU + vecV * vecV + (u01 * u01 + v01 * v01) * y * y ) / ( L * y * area * area );
}
void Global::build_half_mesh()
{
	halfmesh.init( &mesh, &charts );
}

void Global::save_mesh(string name)
{
	printf("saving %s\n", name.c_str());

	FILE *f = fopen(name.c_str(), "wb");

	for (int i = 0; i < mesh.positions.s; i++)
	{
		fprintf(f, "v %f %f %f\n", mesh.positions[i][0], mesh.positions[i][1], mesh.positions[i][2]);
	}
	
	for (int i = 0; i < g.mesh.tex_coords.s; i++)
	{
		fprintf(f, "vt %f %f\n", mesh.tex_coords[i][0], mesh.tex_coords[i][1]);
	}
	
	for (int i = 0; i < mesh.indices_pos.s; i++)
	{
		fprintf(f, "f %d/%d %d/%d %d/%d\n", mesh.indices_pos[i][0]+1, mesh.indices_tex[i][0]+1,
											mesh.indices_pos[i][1]+1, mesh.indices_tex[i][1]+1,
											mesh.indices_pos[i][2]+1, mesh.indices_tex[i][2]+1);
	}

	fclose(f);
}

void Global::save_mesh2(string name)
{
	printf("saving %s\n", name.c_str());

	FILE *f = fopen(name.c_str(), "wb");

	for (int i = 0; i < mesh.positions.s; i++)
	{
		fprintf(f, "v %f %f %f %f %f\n", mesh.positions[i][0], mesh.positions[i][1], mesh.positions[i][2], mesh.tex_coords[i][0], mesh.tex_coords[i][1]);
	}
	
	for (int i = 0; i < g.mesh.indices_pos.s; i++)
	{
		fprintf(f, "f %d %d %d\n", mesh.indices_pos[i][0]+1,mesh.indices_pos[i][1]+1,mesh.indices_pos[i][2]+1);
	}

	fclose(f);
}