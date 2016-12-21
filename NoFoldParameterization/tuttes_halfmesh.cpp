#include "HalfEdge.h"
#include "math.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/Core>

#include <fstream>
#include <sstream>
using namespace Eigen;

//Get angle
vect2d rota(vect2d me, double amt)
{
	vect2d t;
	t[0] = me[0]*cos(amt) - me[1]*sin(amt);
	t[1] = me[1]*cos(amt) + me[0]*sin(amt);
	return t;
}

void HalfEdgeMesh::tuttes_embedding_chart(int c)
{
	int interiorSize;
	double boundaryLength;

	//Get boundary length and use to normalize for arc length param
	boundaryLength = 0;
	
	double tarea = 0;
	// set flattened positions
	printf("Isometric flattening\n");

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

			tarea += 0.5 * flat[1][0] * flat[2][1];
			
			isos.push_back(flat);
		}
	}

	printf("1\n");
	//for (unsigned int it = 0; it < isos.size(); it++)
	//{
	//	tarea += .5 * (isos[it][1][0] * isos[it][2][1]);
	//}

	for(unsigned int i = 0; i < this->boundaries[c].size(); i++)
	{
		boundaryLength += (this->boundaries[c][i]->pos - this->boundaries[c][(i+1)%boundaries[c].size()]->pos).length();
	}

	vector<vector<int>> bverts;
	for(unsigned int i = 0; i < vertexData.size(); i++)
	{
		vector<int> t;
		bverts.push_back(t);
	}


	printf("2\n");
	for(unsigned int i = 0; i < faceData.size(); i++)
	{
		bverts[faceData[i].vi[0]].push_back(faceData[i].vi[1]);
		bverts[faceData[i].vi[1]].push_back(faceData[i].vi[0]);

		bverts[faceData[i].vi[2]].push_back(faceData[i].vi[1]);
		bverts[faceData[i].vi[1]].push_back(faceData[i].vi[2]);

		bverts[faceData[i].vi[0]].push_back(faceData[i].vi[2]);
		bverts[faceData[i].vi[2]].push_back(faceData[i].vi[0]);

		if(vertexData[faceData[i].vi[0]].boundary && vertexData[faceData[i].vi[1]].boundary)
		{
			bverts[faceData[i].vi[0]].push_back(faceData[i].vi[1]);
			bverts[faceData[i].vi[1]].push_back(faceData[i].vi[0]);
		}

		if(vertexData[faceData[i].vi[1]].boundary && vertexData[faceData[i].vi[2]].boundary)
		{
			bverts[faceData[i].vi[1]].push_back(faceData[i].vi[2]);
			bverts[faceData[i].vi[2]].push_back(faceData[i].vi[1]);
		}

		if(vertexData[faceData[i].vi[0]].boundary && vertexData[faceData[i].vi[2]].boundary)
		{
			bverts[faceData[i].vi[0]].push_back(faceData[i].vi[2]);
			bverts[faceData[i].vi[2]].push_back(faceData[i].vi[0]);
		}
	}

	double radius = sqrt(tarea / 3.14159265359);

	vector<vect2d> newbound;
	vect2d t;
	zero_param[0] = t[0] = radius;
	zero_param[1] = t[1] = 0;
	newbound.push_back(t);

	printf("3\n");
	pie2_param = rota(zero_param, 3.14159265359 / 2.0f);

	boundaryLength = boundaries[c].size(); // sahil modified boundary
	for(unsigned int i = 1; i < boundaries[c].size(); i++)
	{
		// more sahil modifications to boundary
		double l = 1; //(boundaries[c][i-1]->pos - boundaries[c][i]->pos).length();
		double rotme = (l*2*3.14159265359) / boundaryLength;

		newbound.push_back(rota(newbound[i-1], rotme) );
	}


	for(unsigned int i = 0; i < boundaries[c].size(); i++)
	{
		boundaries[c][i]->uvpos = newbound[i];

		// WTF? - Sahil
//		boundaries[c][i]->uvpos[0] += 3*c*radius;
	}


	printf("4\n");
	//Setup indices for interior vertices
	interiorSize = 0;
	vector<int> indchange;
	int vInd = 0;
	int boundarySize = 0;
	for (unsigned int i = 0; i < this->vertexData.size(); i++)
	{
		if(vertexData[i].boundary || vertexData[i].chart_ind != c)
		{
			indchange.push_back(-1);
			boundarySize++;
		}
		else
		{
			indchange.push_back(vInd);
			vInd++;
			interiorSize++;
		}
	}


	//Setup matrices
	//A, Bu, Bv, U, V
	vector<Triplet<double>> trips;
	//trips.resize(vertexData.size() * 3);

	printf("5\n");
	SparseMatrix<double> A;
	A.resize(interiorSize, interiorSize);
	//DenseTaucsMatrix Bu, Bv, U, V;

	VectorXd Bu, Bv, U, V;
	Bu.resize(interiorSize);
	//Bu.resize(boundarySize + vertexData.size());
	Bu.setZero();
	Bv.resize(interiorSize);
	//Bv.resize(boundarySize + vertexData.size());
	Bu.setZero();
	

	U.resize(interiorSize);
	U.setZero();
	V.resize(interiorSize);
	V.setZero();

/*	double theta = 0;
	for(int i = vertexData.size(); i < boundarySize + vertexData.size(); i++)
	{
		Bu(i) = radius * cos(theta);
		Bv(i) = radius * sin(theta);

		theta += 2 * PI / boundarySize;
	}

	IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::ofstream file("test.txt");
	file << Bu.format(CleanFmt);
	file.close()*/;

	printf("6\n");
	for(int i = 0; i < interiorSize; i++)
	{
		Bu(i) = 0;
		Bv(i) = 0;
	}

	vInd = 0;
	for (unsigned int i = 0; i < vertexData.size(); i++)
	{
		if(vertexData[i].boundary || vertexData[i].chart_ind != c)
		{
			//Does not effect A
		}
		else
		{
			trips.push_back(Triplet<double>(vInd,vInd,
						.5*bverts[i].size()));
			if ( vInd >= interiorSize )
			{
				printf ( "error %d\n", vInd );
			}

			float d = .5f /*/  (float)bverts[i].size()*/;
			for(unsigned int j = 0; j < bverts[i].size(); j++)
			{
				if ( vInd >= interiorSize )
				{
					printf ( "error %d\n", indchange[bverts[i][j]] );
				}
					if(indchange[bverts[i][j]] != -1)//Make A
						trips.push_back(Triplet<double>(vInd,indchange[bverts[i][j]],
							-1*d));
				else
				{
					Bu(vInd, 0) += .5 * vertexData[bverts[i][j]].uvpos[0];
					Bv(vInd, 0) += .5 * vertexData[bverts[i][j]].uvpos[1];
				}
			}

				vInd++;
		}
	}

	printf("7\n");
	A.setFromTriplets(trips.begin(), trips.end());
	printf("8\n");
//	A.makeCompressed();
	SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solverA;
	printf("SOLVE\n");
//	solverA.compute(A);
	solverA.analyzePattern(A);
	solverA.factorize(A);

	printf("SOLVE U\n");
	U = solverA.solve(Bu);

	printf("SOLVE V\n");
	V = solverA.solve(Bv);

	IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::ofstream file("test.txt");
	file << Bu.format(CleanFmt);
	file.close();
//	printf("%f, %d %d", Bu(1), Bu.size(), Bv.size());

	for(unsigned int i = 0; i < indchange.size(); i++)
	{
		if(indchange[i] != -1)
		{
			//Not boundary
			vertexData[i].uvpos[0] = U(indchange[i],0);
			vertexData[i].uvpos[1] = V(indchange[i],0);
		}
	}

	/*std::ifstream uvpos("cubeuv.txt");
	for(int i = 0; i < vertexData.size(); i++)
	{
		string line;
		getline(uvpos, line);
		line = line.substr(1, line.size() - 2);
		std::istringstream sstream(line);

		string num;
		getline(sstream, num, ',');
		float u = stof(num);
		getline(sstream, num, ',');
		float v = stof(num);

		vertexData[i].uvpos.v[0] = u;
		vertexData[i].uvpos.v[1] = v;
	}*/
}