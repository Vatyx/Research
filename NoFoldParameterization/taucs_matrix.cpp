#include "taucs_matrix.h"
#include <fstream>
#include <stdio.h>

SparseTaucsMatrix::~SparseTaucsMatrix()
{
	clear();
}

void SparseTaucsMatrix::clear()
{
	if (L) taucs_ccs_free(L);
	if (Aod) taucs_ccs_free(Aod);

	L = 0;
	Aod = 0;
	
	delete [] A.colptr;
	delete [] A.rowind;
	delete [] A.values.d;

	A.colptr = 0;
	A.rowind = 0;
	A.values.d = 0;

	free(perm);
	free(invperm);

	perm = 0;
	invperm = 0;
}

void SparseTaucsMatrix::save_ccs(string fn)
{
	// build arrays
	Array<__int32> colptr;
	Array<__int32> rowind;
	Array<double> vals;

	colptr.push_back(0);

	int elems = 0;
	for (int c = 0; c < data.size(); c++)
	{
		for (map<int, double>::iterator it  = data[c].begin(); it != data[c].end(); ++it)
		{
			rowind.push_back(it->first);
			vals.push_back(it->second);
			elems++;
		}
		colptr.push_back(elems);
	}

	// write arrays
	ofstream out(fn.c_str(), ios::binary);

	__int32 n = data.size();
	out.write((char*)&n, 4);

	out.write((char*)&colptr[0], 4*colptr.size());
	out.write((char*)&rowind[0], 4*rowind.size());
	out.write((char*)&vals[0], sizeof(double)*vals.size());

	out.close();
}

void SparseTaucsMatrix::make_ccs()
{
	A.m = A.n = data.size();
	A.flags = TAUCS_TRIANGULAR | TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;

	// calc size of arrays and allocate
	int num_vals = 0;
	for (int c = 0; c < data.size(); c++)
		num_vals += data[c].size();

	delete [] A.colptr;
	delete [] A.rowind;
	delete [] A.values.d;

	A.colptr = new int [data.size()+1];
	A.rowind = new int [num_vals];
	A.values.d = new double [num_vals];

	// build arrays
	A.colptr[0] = 0;

	int elem = 0;
	for (int c = 0; c < data.size(); c++)
	{
		for (map<int, double>::iterator it  = data[c].begin(); it != data[c].end(); ++it)
		{
			A.rowind[elem] = it->first;
			A.values.d[elem] = it->second;

			elem++;
		}
		A.colptr[c+1] = elem;
	}
}

void SparseTaucsMatrix::factorize()
{
	taucs_logfile("stdout");

	taucs_ccs_order(&A, &perm, &invperm, "metis");
	Aod = taucs_ccs_permute_symmetrically(&A, perm, invperm);

	void *F = taucs_ccs_factor_llt_mf(Aod);
	L = taucs_supernodal_factor_to_ccs(F);
	taucs_supernodal_factor_free(F);
}

void SparseTaucsMatrix::solve(DenseTaucsMatrix &X, DenseTaucsMatrix &B)
{
	if (!L)
		return;

	X.resize(B.rows, B.cols);

	double *Xod = new double [L->n];
	double *Bod = new double [L->n];

	for (int i = 0; i < B.cols; i++)
	{
		double *x = &X.data[L->n*i];
		double *b = &B.data[L->n*i];

		taucs_vec_permute(L->n, TAUCS_DOUBLE, b, Bod, perm);

		taucs_ccs_solve_llt(L, Xod, Bod); 
		
		taucs_vec_ipermute(L->n, TAUCS_DOUBLE, Xod, x, perm);
	}

	delete [] Xod;
	delete [] Bod;
}

void SparseTaucsMatrix::solve_cg(DenseTaucsMatrix &X, DenseTaucsMatrix &B)
{
	//taucs_logfile("stdout");

	if (X.cols != B.cols || X.rows != B.rows)
	{
		X.resize(B.rows, B.cols);
		for (int i = 0; i < B.rows*B.cols; i++)
			X.data[i] = 0;
	}

	for (int i = 0; i < B.cols; i++)
	{
		double *x = &X.data[B.rows*i];
		double *b = &B.data[B.rows*i];
		
		taucs_conjugate_gradients(&A, 0, 0, x, b, 1000, 1e-4);
	}
}

void SparseTaucsMatrix::save_matlab(string fn, string var_name)
{
	FILE *f = fopen(fn.c_str(), "wb");
	
	Array<int> rr, cc;
	Array<double> vv;

	for (int r = 0; r < data.size(); r++) {
		for (map<int, double>::iterator it = data[r].begin(); it != data[r].end(); ++it) {
			int c = it->first;
			double v = it->second;

			if (abs(v) < 1e-9)
				continue;
			
			rr.push_back(r+1);
			cc.push_back(c+1);
			vv.push_back(v);
			if (c != r) {
				rr.push_back(c+1);
				cc.push_back(r+1);
				vv.push_back(v);
			}
		}
	}
	
	fprintf(f, "mat_rows = [");
	for (int i = 0; i < rr.s; i++) {
		fprintf(f, " %d", rr[i]);
	}
	fprintf(f, "];\n");
	
	fprintf(f, "mat_cols = [");
	for (int i = 0; i < rr.s; i++) {
		fprintf(f, " %d", cc[i]);
	}
	fprintf(f, "];\n");
	
	fprintf(f, "mat_vals = [");
	for (int i = 0; i < rr.s; i++) {
		fprintf(f, " %18.16g", vv[i]);
	}
	fprintf(f, "];\n");
	fprintf(f, "%s = sparse(mat_rows, mat_cols, mat_vals, %d, %d);\n", var_name.c_str(), data.size(), data.size());
	fclose(f);
}

void SparseTaucsMatrix::save_mathematica(string fn)
{
	FILE *f = fopen(fn.c_str(), "wb");

	fprintf(f, "SparseArray[{");
	bool first = true;
	for (int r = 0; r < data.size(); r++)
	{
		for (map<int, double>::iterator it = data[r].begin(); it != data[r].end(); ++it)
		{
			int c = it->first;
			double v = it->second;

			if (!first)
				fprintf(f, ",");
			else
				first = false;

			fprintf(f, "{%d,%d}->%f", r+1, c+1, v);
			if (c != r)
				fprintf(f, ",{%d,%d}->%f", c+1, r+1, v);
		}
	}
	fprintf(f, "}]");

	fclose(f);
}

void DenseTaucsMatrix::save_mathematica(string fn)
{
	FILE *f = fopen(fn.c_str(), "wb");

	fprintf(f, "{");
	for (int r = 0; r < rows; r++)
	{
		fprintf(f, "{");
		for (int c = 0; c < cols; c++)
		{
			fprintf(f, "%f", (*this)(r,c));

			if (c < cols-1)
				fprintf(f, ",");
		}
		fprintf(f, "}");

		if (r < rows-1)
			fprintf(f, ",");
	}
	fprintf(f, "}");

	fclose(f);
}

void save_dense_matrix(DenseTaucsMatrix &x, string fn)
{
	ofstream out(fn.c_str(), ios::binary);

	__int32 n = x.rows;
	out.write((char*)&n, 4);
	__int32 m = x.cols;
	out.write((char*)&m, 4);

	out.write((char*)x.data, 8*n*m);

	out.close();
}

void load_dense_matrix(DenseTaucsMatrix &x, string fn)
{
	ifstream in(fn.c_str(), ios::binary);

	in.read((char*)&x.rows, 4);
	in.read((char*)&x.cols, 4);

	x.data = new double [x.rows*x.cols];
	in.read((char*)x.data, 8*x.rows*x.cols);

	in.close();
}
