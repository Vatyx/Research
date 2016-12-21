#ifndef MATRIX_H
#define MATRIX_H

//#include "../stdafx.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>

namespace MathMatrix
{

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

template <class Type>
static void jacobi ( Type **u, Type *d, Type **v, int n )
{
	int j, iq, ip, i, k;
	Type tresh, theta, tau, t, sm, s, h, g, c;
	Type *b = new Type [ n ];
	Type *z = new Type [ n ];
	Type **a;

	a = new Type *[ n ];
	for ( i = 0; i < n; i++ )
	{
		a [ i ] = new Type [ n ];
		for ( j = 0; j < n; j++ )
		{
			a [ i ] [ j ] = u [ i ] [ j ];
		}
	}

	for ( ip = 0; ip < n; ip++ ) 
	{
		for ( iq = 0; iq < n; iq++ )
		{
			v [ ip ] [ iq ] = 0.0f;
		}
		v [ ip ] [ ip ] = 1.0f;
	}

	for ( ip = 0; ip < n; ip++ )
	{
		b [ ip ] = a [ ip ] [ ip ];
		d [ ip ] = b [ ip ];
		z [ ip ] = 0.0f;
	}

	for ( i = 1; i <= 100; i++ )
	{
		sm = 0.0f;
		for ( ip = 0; ip < n - 1; ip++ )
		{
			for ( iq = ip + 1; iq < n; iq++ )
			{
				sm += (Type)fabs ( a [ ip ] [ iq ] );
			}
		}

		if ( sm == 0.0f )
		{
			// sort the stupid things and transpose

			// transpose
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					a [ i ] [ j ] = v [ j ] [ i ];
				}
			}

			// stupid sort n^2... however n should always be small
			// bubble sort
			// largest first
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n - i - 1; j++ )
				{
					if ( d [ j ] < d [ j + 1 ] )
//					if ( fabs ( d [ j ] ) < fabs ( d [ j + 1 ] ) )
					{
						sm = d [ j ];
						d [ j ] = d [ j + 1 ];
						d [ j + 1 ] = sm;

						for ( k = 0; k < n; k++ )
						{
							sm = a [ j ] [ k ];
							a [ j ] [ k ] = a [ j + 1 ] [ k ];
							a [ j + 1 ] [ k ] = sm;
						}
					}
				}
			}

			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					v [ i ] [ j ] = a [ i ] [ j ];
				}
				delete[] a [ i ];
			}
			delete[] a;
			delete[] b;
			delete[] z;

			return;
		}

		if ( i < 4 )
		{
			tresh = 0.2f * sm / ( n * n );
		}
		else
		{
			tresh = 0.0f;
		}

		for ( ip = 0; ip < n - 1; ip++ )
		{
			for ( iq = ip + 1; iq < n; iq++ ) 
			{
				g = 100.0f * (Type)fabs ( a [ ip ] [ iq ] );
				if ( i > 4 && ( fabs ( d [ ip ] ) + g ) == fabs ( d [ ip ] )
					&& ( fabs ( d [ iq ] ) + g ) == fabs ( d [ iq ] ) )
				{
					a [ ip ] [ iq ] = 0.0f;
				}
				else
				{
					if ( fabs ( a [ ip ] [ iq ] ) > tresh )
					{
						h = d [ iq ] - d [ ip ];
						if ( ( fabs ( h ) + g ) == fabs ( h ) )
						{
							t = ( a [ ip ] [ iq ] ) / h;
						}
						else
						{
							theta = 0.5f * h / ( a [ ip ] [ iq ] );
							t = 1.0f / (Type)( fabs ( theta ) + sqrt ( 1.0f + theta * theta ) );
							if ( theta < 0.0f ) 
							{
								t = -1.0f * t;
							}
						}

						c = 1.0f / (Type)sqrt ( 1 + t * t );
						s = t * c;
						tau = s / ( 1.0f + c );
						h = t * a [ ip ] [ iq ];
						z [ ip ] -= h;
						z [ iq ] += h;
						d [ ip ] -= h;
						d [ iq ] += h;
						a [ ip ] [ iq ] = 0.0f;
						for ( j = 0; j <= ip - 1; j++ )
						{
							ROTATE ( a, j, ip, j, iq )
						}
						for ( j = ip + 1; j <= iq - 1; j++ )
						{
							ROTATE ( a, ip, j, j, iq )
						}
						for ( j = iq + 1; j < n; j++ )
						{
							ROTATE ( a, ip, j, iq, j )
						}
						for ( j = 0; j < n; j++ )
						{
							ROTATE ( v, j, ip, j, iq )
						}
					}
				}
			}
		}

		for ( ip = 0; ip < n; ip++ )
		{
			b [ ip ] += z [ ip ];
			d [ ip ] = b [ ip ];
			z [ ip ] = 0.0f;
		}
	}

	assert ( false );
	// transpose
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			a [ i ] [ j ] = v [ j ] [ i ];
		}
	}

	// stupid sort n^2... however n should always be small
	// bubble sort
	// largest first
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n - i - 1; j++ )
		{
			if ( d [ j ] < d [ j + 1 ] )
//					if ( fabs ( d [ j ] ) < fabs ( d [ j + 1 ] ) )
			{
				sm = d [ j ];
				d [ j ] = d [ j + 1 ];
				d [ j + 1 ] = sm;

				for ( k = 0; k < n; k++ )
				{
					sm = a [ j ] [ k ];
					a [ j ] [ k ] = a [ j + 1 ] [ k ];
					a [ j + 1 ] [ k ] = sm;
				}
			}
		}
	}

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			v [ i ] [ j ] = a [ i ] [ j ];
		}
		delete[] a [ i ];
	}
	delete[] b;
	delete[] z;

	printf ( "too many iterations in jacobi\n" );
}

template <class Type, int n>
static void jacobi ( Type **u, Type *d, Type **v )
{
	int j, iq, ip, i, k;
	Type tresh, theta, tau, t, sm, s, h, g, c;
	Type b [ n ];
	Type z [ n ];
	Type a [ n ] [ n ];

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			a [ i ] [ j ] = u [ i ] [ j ];
		}
	}

	for ( ip = 0; ip < n; ip++ ) 
	{
		for ( iq = 0; iq < n; iq++ )
		{
			v [ ip ] [ iq ] = 0.0f;
		}
		v [ ip ] [ ip ] = 1.0f;
	}

	for ( ip = 0; ip < n; ip++ )
	{
		b [ ip ] = a [ ip ] [ ip ];
		d [ ip ] = b [ ip ];
		z [ ip ] = 0.0f;
	}

	for ( i = 1; i <= 100; i++ )
	{
		sm = 0.0f;
		for ( ip = 0; ip < n - 1; ip++ )
		{
			for ( iq = ip + 1; iq < n; iq++ )
			{
				sm += (Type)fabs ( a [ ip ] [ iq ] );
			}
		}

		if ( sm == 0.0f )
		{
			// sort the stupid things and transpose

			// transpose
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					a [ i ] [ j ] = v [ j ] [ i ];
				}
			}

			// stupid sort n^2... however n should always be small
			// bubble sort
			// largest first
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n - i - 1; j++ )
				{
					if ( d [ j ] < d [ j + 1 ] )
//					if ( fabs ( d [ j ] ) < fabs ( d [ j + 1 ] ) )
					{
						sm = d [ j ];
						d [ j ] = d [ j + 1 ];
						d [ j + 1 ] = sm;

						for ( k = 0; k < n; k++ )
						{
							sm = a [ j ] [ k ];
							a [ j ] [ k ] = a [ j + 1 ] [ k ];
							a [ j + 1 ] [ k ] = sm;
						}
					}
				}
			}

			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					v [ i ] [ j ] = a [ i ] [ j ];
				}
			}

			return;
		}

		if ( i < 4 )
		{
			tresh = 0.2f * sm / ( n * n );
		}
		else
		{
			tresh = 0.0f;
		}

		for ( ip = 0; ip < n - 1; ip++ )
		{
			for ( iq = ip + 1; iq < n; iq++ ) 
			{
				g = 100.0f * (Type)fabs ( a [ ip ] [ iq ] );
				if ( i > 4 && ( fabs ( d [ ip ] ) + g ) == fabs ( d [ ip ] )
					&& ( fabs ( d [ iq ] ) + g ) == fabs ( d [ iq ] ) )
				{
					a [ ip ] [ iq ] = 0.0f;
				}
				else
				{
					if ( fabs ( a [ ip ] [ iq ] ) > tresh )
					{
						h = d [ iq ] - d [ ip ];
						if ( ( fabs ( h ) + g ) == fabs ( h ) )
						{
							t = ( a [ ip ] [ iq ] ) / h;
						}
						else
						{
							theta = 0.5f * h / ( a [ ip ] [ iq ] );
							t = 1.0f / (Type)( fabs ( theta ) + sqrt ( 1.0f + theta * theta ) );
							if ( theta < 0.0f ) 
							{
								t = -1.0f * t;
							}
						}

						c = 1.0f / (Type)sqrt ( 1 + t * t );
						s = t * c;
						tau = s / ( 1.0f + c );
						h = t * a [ ip ] [ iq ];
						z [ ip ] -= h;
						z [ iq ] += h;
						d [ ip ] -= h;
						d [ iq ] += h;
						a [ ip ] [ iq ] = 0.0f;
						for ( j = 0; j <= ip - 1; j++ )
						{
							ROTATE ( a, j, ip, j, iq )
						}
						for ( j = ip + 1; j <= iq - 1; j++ )
						{
							ROTATE ( a, ip, j, j, iq )
						}
						for ( j = iq + 1; j < n; j++ )
						{
							ROTATE ( a, ip, j, iq, j )
						}
						for ( j = 0; j < n; j++ )
						{
							ROTATE ( v, j, ip, j, iq )
						}
					}
				}
			}
		}

		for ( ip = 0; ip < n; ip++ )
		{
			b [ ip ] += z [ ip ];
			d [ ip ] = b [ ip ];
			z [ ip ] = 0.0f;
		}
	}

	assert ( false );
	// transpose
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			a [ i ] [ j ] = v [ j ] [ i ];
		}
	}

	// stupid sort n^2... however n should always be small
	// bubble sort
	// largest first
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n - i - 1; j++ )
		{
			if ( d [ j ] < d [ j + 1 ] )
//					if ( fabs ( d [ j ] ) < fabs ( d [ j + 1 ] ) )
			{
				sm = d [ j ];
				d [ j ] = d [ j + 1 ];
				d [ j + 1 ] = sm;

				for ( k = 0; k < n; k++ )
				{
					sm = a [ j ] [ k ];
					a [ j ] [ k ] = a [ j + 1 ] [ k ];
					a [ j + 1 ] [ k ] = sm;
				}
			}
		}
	}

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			v [ i ] [ j ] = a [ i ] [ j ];
		}
	}

	printf ( "too many iterations in jacobi\n" );
}

template <class Type>
static void matInverse ( Type **mat, Type **rvalue, int n, Type tolerance = 0.00000001 )
{
	// there is an implicit assumption that mat is symmetric and real
	// U and V in the SVD will then be the same matrix whose rows are the eigenvectors of mat
	// W will just be the eigenvalues of mat
	Type *w = new Type [ n ];
	Type **u;
	int i, j, k;

	u = new Type *[ n ];
	for ( i = 0; i < n; i++ )
	{
		u [ i ] = new Type [ n ];
	}

	jacobi<Type> ( mat, w, u, n );

	if ( w [ 0 ] == 0.0f )
	{
//		printf ( "error: largest eigenvalue is 0!\n" );
	}
	else
	{
		for ( i = 1; i < n; i++ )
		{
			if ( fabs ( (double)w [ i ] / w [ 0 ] ) < tolerance ) // / w [ 0 ] < TOLERANCE )
			{
					w [ i ] = 0;
			}
			else
			{
				w [ i ] = 1.0f / w [ i ];
			}
		}
		w [ 0 ] = 1.0f / w [ 0 ];
	}

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			rvalue [ i ] [ j ] = 0;
			for ( k = 0; k < n; k++ )
			{
				rvalue [ i ] [ j ] += w [ k ] * u [ k ] [ i ] * u [ k ] [ j ];
			}
		}
	}

	for ( i = 0; i < n; i++ )
	{
		delete[] u [ i ];
	}
	delete[] u;
}

#undef ROTATE







template <class type>
class Matrix
{
//friend Matrix<float> constructMatrix ( MathVector3 *charMap, int *indices, int n );

protected:
	type **data;
	int r, c;

public:
	Matrix ( ) : r ( 0 ), c ( 0 ), data ( NULL ) { }

	Matrix ( int rows, int columns )
	{
		r = rows;
		c = columns;

		data = new type *[ rows ];
		for ( int i = 0; i < rows; i++ )
		{
			data [ i ] = new type [ columns ];
		}
	}

	Matrix ( const Matrix<type> &mat )
	{
		int i;

		r = mat.r;
		c = mat.c;

		data = new type *[ r ];
		for ( i = 0; i < r; i++ )
		{
			data [ i ] = new type [ c ];
			memcpy ( data [ i ], mat.data [ i ], sizeof ( type ) * c );
		}
	}

	Matrix ( type **fill, int rows, int columns )
	{
		r = rows;
		c = columns;

		data = new type *[ rows ];
		for ( int i = 0; i < rows; i++ )
		{
			data [ i ] = new type [ columns ];
			memcpy ( data [ i ], fill [ i ], sizeof ( type ) * columns );
		}
	}

	void fillWith ( type **fill )
	{
		for ( int i = 0; i < rows; i++ )
		{
			memcpy ( data [ i ], fill [ i ], sizeof ( type ) * columns );
		}
	}

	virtual ~Matrix ( )
	{
		if ( data )
		{
			for ( int i = 0; i < r; i++ )
			{
				delete[] data [ i ];
			}
			delete[] data;
		}
	}


	virtual Matrix<type>& zero ( void )
	{
		int i, j;

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				data [ i ] [ j ] = 0;
			}
		}

		return *this;
	}


	virtual Matrix<type>& identity ( void )
	{
		int i;

		if ( r == c )
		{
			zero ( );

			for ( i = 0; i < r; i++ )
			{
				data [ i ] [ i ] = 1;
			}
		}

		return *this;
	}

	int getNumRows ( void ) const { return r; }
  
	int getNumCols ( void ) const { return c; }

	type **getData ( void ) { return data; }

	virtual int operator = ( const Matrix<type> &mat )
	{
		int i;

		if ( data )
		{
			for ( i = 0; i < r; i++ )
			{
				delete[] data [ i ];
			}
			delete[] data;
		}

		r = mat.r;
		c = mat.c;

		data = new type *[ r ];
		for ( i = 0; i < r; i++ )
		{
			data [ i ] = new type [ c ];
			memcpy ( data [ i ], mat.data [ i ], sizeof ( type ) * c );
		}

		return 0;
	}

	virtual Matrix<type> operator * ( Matrix<type> &mat ) const
	{
		assert ( c == mat.r );
		int i, j, k;
		Matrix<type> rvalue ( r, mat.c );

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < mat.c; j++ )
			{
				rvalue [ i ] [ j ] = 0;

				for ( k = 0; k < c; k++ )
				{
					rvalue [ i ] [ j ] += data [ i ] [ k ] * mat [ k ] [ j ];
				}
			}
		}

		return rvalue;
	}

	virtual Matrix<type> operator * ( type val ) const
	{
		Matrix<type> rvalue ( *this );

		rvalue *= val;

		return rvalue;
	}

	virtual void operator *= ( Matrix<type> &mat )
	{
		assert ( c == mat.r );

		*this = *this * mat;
	}

	virtual void operator *= ( type val )
	{
		int i, j;

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				data [ i ] [ j ] *= val;
			}
		}
	}

	virtual Matrix<type> operator + ( Matrix<type> &mat ) const
	{
		assert ( r == mat.r && c == mat.c );
		int i, j;
		Matrix<type> rvalue ( r, c );

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				rvalue [ i ] [ j ] = data [ i ] [ j ] + mat [ i ] [ j ];
			}
		}

		return rvalue;
	}

	virtual Matrix<type> operator + ( type val ) const
	{
		Matrix<type> rvalue ( *this );

		rvalue += val;

		return rvalue;
	}

	virtual void operator += ( Matrix<type> &mat )
	{
		assert ( r == mat.r && c == mat.c );

		*this = *this + mat;
	}

	virtual void operator += ( type val )
	{
		int i, j;

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				data [ i ] [ j ] += val;
			}
		}
	}

	virtual Matrix<type> operator - ( Matrix<type> &mat ) const
	{
		assert ( r == mat.r && c == mat.c );
		int i, j;
		Matrix<type> rvalue ( r, c );

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				rvalue [ i ] [ j ] = data [ i ] [ j ] - mat [ i ] [ j ];
			}
		}

		return rvalue;
	}

	virtual Matrix<type> operator - ( type val ) const
	{
		Matrix<type> rvalue ( *this );

		rvalue -= val;

		return rvalue;
	}

	virtual void operator -= ( Matrix<type> &mat )
	{
		assert ( r == mat.r && c == mat.c );

		*this = *this - mat;
	}

	virtual void operator -= ( type val )
	{
		int i, j;

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				data [ i ] [ j ] -= val;
			}
		}
	}

	virtual type *operator [] ( int index )
	{
		assert ( index >= 0 && index < r );

		return data [ index ];
	}

	virtual Matrix<type> transpose ( void ) const
	{
		Matrix<type> rvalue ( c, r );
		int i, j;

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < c; j++ )
			{
				rvalue [ j ] [ i ] = data [ i ] [ j ];
			}
		}

		return rvalue;
	}

	virtual void solve ( Matrix<type> &mat )
	{
		assert ( r == c && mat.r == r );
		int *pvts = new int [ r ];
		int temp;
		int i, j, k;

		for ( i = 0; i < r; i++ )
		{
			pvts [ i ] = i;
		}

		for ( i = 0; i < r; i++ )
		{
			int max = i;

			for ( j = i + 1; j < r; j++ )
			{
				if ( data [ pvts [ max ] ] [ i ] < data [ pvts [ j ] ] [ i ] )
				{
					max = j;
				}
			}

			temp = pvts [ i ];
			pvts [ i ] = pvts [ max ];
			pvts [ max ] = temp;

//			assert ( fabs ( data [ pvts [ i ] ] [ i ] ) >= 0.0000001 );

			for ( j = i + 1; j < c; j++ )
			{
				data [ pvts [ i ] ] [ j ] /= data [ pvts [ i ] ] [ i ];
			}
			for ( j = 0; j < mat.c; j++ )
			{
				mat [ pvts [ i ] ] [ j ] /= data [ pvts [ i ] ] [ i ];
			}
			data [ pvts [ i ] ] [ i ] = 1;

			for ( j = i + 1; j < r; j++ )
			{
				for ( k = i + 1; k < c; k++ )
				{
					data [ pvts [ j ] ] [ k ] -= data [ pvts [ i ] ] [ k ] * data [ pvts [ j ] ] [ i ];
				}

				for ( k = 0; k < mat.c; k++ )
				{
					mat [ pvts [ j ] ] [ k ] -= mat [ pvts [ i ] ] [ k ] * data [ pvts [ j ] ] [ i ];
				}
				data [ pvts [ j ] ] [ i ] = 0;
			}
		}

		// we should be upper triangular by now with 1 along diagonal... back substitute
		for ( i = r - 1; i > 0; i-- )
		{
			for ( j = i - 1; j >= 0; j-- )
			{
				for ( k = 0; k < mat.c; k++ )
				{
					mat [ pvts [ j ] ] [ k ] -= mat [ pvts [ i ] ] [ k ] * data [ pvts [ j ] ] [ i ];
				}
			}
		}

		Matrix<type> tempMat ( mat.r, mat.c );

		for ( i = 0; i < mat.r; i++ )
		{
			memcpy ( tempMat [ i ], mat [ pvts [ i ] ], sizeof ( type ) * mat.c );
		}
		for ( i = 0; i < mat.r; i++ )
		{
			memcpy ( mat [ i ], tempMat [ i ], sizeof ( type ) * mat.c );
		}
		delete[] pvts;
	}

	virtual Matrix<type> pseudoInverse ( type tolerance = 0.00000001 )
	{
		assert ( r == c );
		Matrix<type> rvalue ( r, r );

		matInverse<type> ( data, rvalue.data, r, tolerance );

		return rvalue;
	}

	virtual void outputMathematica ( char *filename )
	{
		FILE *fptr = fopen ( filename, "wt" );
		int i, j;

		if ( fptr == NULL )
		{
			fprintf ( stderr, "Error: could not open %s for writing\n", filename );
			return;
		}

		fprintf ( fptr, "{" );
		for ( i = 0; i < r - 1; i++ )
		{
			fprintf ( fptr, "{" );
			for ( j = 0; j < c - 1; j++ )
			{
				fprintf ( fptr, "%.3f,", data [ i ] [ j ] );
			}
			fprintf ( fptr, "%.3f},", data [ i ] [ j ] );
		}
		fprintf ( fptr, "{" );
		for ( j = 0; j < c - 1; j++ )
		{
			fprintf ( fptr, "%.3f,", data [ i ] [ j ] );
		}
		fprintf ( fptr, "%.3f}}", data [ i ] [ j ] );

		fclose ( fptr );
	}
};

template <class type>
class sVector : public Matrix<type>
{
public:
	sVector ( ) : Matrix<type> ( )
	{
	}

	sVector ( int n ) : Matrix<type> ( 1, n )
	{
	}

	sVector ( type *fill, int n ) : Matrix<type> ( 1, n )
	{
		int i;

		for ( i = 0; i < n; i++ )
		{
			data [ 0 ] [ i ] = fill [ i ];
		}
	}

	sVector ( const Matrix<type> &mat ) : Matrix<type> ( mat )
	{
	}

	virtual ~sVector ( void )
	{
	}

	int getNumRows ( void ) { return r; }

	int getNumCols ( void ) { return c; }

	int operator = ( const sVector<type> &v )
	{
		return Matrix<type>::operator = ( (const Matrix<type> &)v );
	}

	sVector<type> operator * ( Matrix<type> &mat )
	{
		return sVector<type> ( ((Matrix<type> *)this)->operator * ( mat ) );
	}

	sVector<type> operator * ( type val )
	{
		return sVector<type> ( ((Matrix<type> *)this)->operator * ( val ) );
	}

	void operator *= ( Matrix<type> &mat )
	{
		((Matrix<type> *)this)->operator *= ( mat );
	}

	void operator *= ( type val )
	{
		((Matrix<type> *)this)->operator *= ( val );
	}

	sVector<type> operator + ( sVector<type> &mat )
	{
		return sVector<type> ( ((Matrix<type> *)this)->operator + ( (Matrix<type> &)mat ) );
	}

	sVector<type> operator + ( type val )
	{
		return sVector<type> ( ((Matrix<type> *)this)->operator + ( val ) );
	}

	void operator += ( sVector<type> &mat )
	{
		((Matrix<type> *)this)->operator += ( (Matrix<type> &)mat );
	}

	void operator += ( type val )
	{
		((Matrix<type> *)this)->operator += ( val );
	}

	sVector<type> operator - ( sVector<type> &mat )
	{
		return sVector<type> ( ((Matrix<type> *)this)->operator - ( (Matrix<type> &)mat ) );
	}

	sVector<type> operator - ( type val )
	{
		return sVector<type> ( ((Matrix<type> *)this)->operator - ( val ) );
	}

	void operator -= ( sVector<type> &mat )
	{
		((Matrix<type> *)this)->operator -= ( (Matrix<type> &)mat );
	}

	void operator -= ( type val )
	{
		((Matrix<type> *)this)->operator -= ( val );
	}

	type length ( void )
	{
		int i;
		type rvalue = 0;

		for ( i = 0; i < c; i++ )
		{
			rvalue += data [ 0 ] [ i ] * data [ 0 ] [ i ];
		}

		return sqrt ( rvalue );
	}

	type &operator [] ( int index )
	{
		assert ( index >= 0 && index < c );

		return data [ 0 ] [ index ];
	}
};
/*
template <class type>
ostream& operator<< ( ostream& os, Matrix<type>& mat )
{
	int i, j;

	for ( i = 0; i < mat.getNumRows ( ); i++ )
	{
		for ( j = 0; j < mat.getNumCols ( ); j++ )
		{
			os << mat [ i ] [ j ] << " ";
		}
		os << endl;
	}

	return os;
}
*/
static char getChar ( FILE *fptr )
{
	char token;
	
	if ( fread ( &token, 1, 1, fptr ) != 1 )
	{
		// shit... end of file
		return 0;
	}

	while ( token == ' ' || token == '\r' ||
		token == '\n' || token == '\t' ) {

		if ( fread ( &token, 1, 1, fptr ) != 1 )
		{
			// shit... end of file
			return 0;
		}
	}
	
	return token;
}

template <class Type>
Matrix<Type> loadMatrix ( char *filename )
{
	FILE *fptr = fopen ( filename, "rt" );
	char buffer [ 256 ];
	int pos;
	char token;
	vector< vector<Type> > data;
	Matrix<Type> rvalue;
	int i, j;

	token = getChar ( fptr );
	if ( token != '{' )
	{
		cout << "Error: expected { but read " << token << endl;
	}

	token = getChar ( fptr );
	while ( token == '{' )
	{
		vector<Type> row;

		while ( token != '}' )
		{
			pos = 0;
			memset ( buffer, 0, 256 );

			token = getChar ( fptr );
			while ( token != ',' && token != '}' ) {
				buffer [ pos++ ] = token;
				token = getChar ( fptr );
			}
			buffer [ pos ] = '\0';

			row.push_back ( atof ( buffer ) );
		}
		// consume next token (either } or ,)
		token = getChar ( fptr );

		token = getChar ( fptr );

		data.push_back ( row );
	}

	rvalue = Matrix<Type> ( data.size ( ), data [ 0 ].size ( ) );
	for ( i = 0; i < rvalue.getNumRows ( ); i++ )
	{
		vector<Type> &tempRow = data [ i ];

		assert ( tempRow.size ( ) == rvalue.getNumCols ( ) );

		for ( j = 0; j < rvalue.getNumCols ( ); j++ )
		{
			rvalue [ i ] [ j ] = tempRow [ j ];
		}
	}

	return rvalue;
}
}
#endif