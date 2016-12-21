#include "QEF.h"
#include <stdlib.h>
#include <math.h>
#include <memory.h>

extern FILE *trace_file;

double crapD [ 4 ];
QEFNormal<double, 3> temp1 ( crapD );
QEFQR<double, 3> temp2 ( crapD );
QEFCenter<double, 3> temp5 ( crapD );

float crapF [ 4 ];
QEFNormal<float, 3> temp3 ( crapF );
QEFQR<float, 3> temp4 ( crapF );

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

// this is the worst kind of hack, but the VC++ compiler can't do templates correctly
// In the case below, the compiler seems to be treating the variable "n" as a dynamic
// variable and not using the code transformation the template defines.
/*
template <class Type, int n>
void jacobi ( Type u[n][n], Type *d, Type v[n][n]);  <- Doesn't compile in VC++
*/

template <class Type, int n>
void jacobi ( ArrayWrapper<Type, n> &u, Type *d, ArrayWrapper<Type, n> &v )
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
			a [ i ] [ j ] = u.data [ i ] [ j ];
		}
	}

	for ( ip = 0; ip < n; ip++ ) 
	{
		for ( iq = 0; iq < n; iq++ )
		{
			v.data [ ip ] [ iq ] = 0.0f;
		}
		v.data [ ip ] [ ip ] = 1.0f;
	}

	for ( ip = 0; ip < n; ip++ )
	{
		b [ ip ] = a [ ip ] [ ip ];
		d [ ip ] = b [ ip ];
		z [ ip ] = 0.0f;
	}

	for ( i = 1; i <= 50; i++ )
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
					a [ i ] [ j ] = v.data [ j ] [ i ];
				}
			}

			// stupid sort n^2... however n should always be small
			// bubble sort
			// largest first
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n - i - 1; j++ )
				{
					if ( fabs ( d [ j ] ) < fabs ( d [ j + 1 ] ) )
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
					v.data [ i ] [ j ] = a [ i ] [ j ];
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
							ROTATE ( v.data, j, ip, j, iq )
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

	printf ( "too many iterations in jacobi\n" );
	exit ( 1 );
}

template <class Type, int n>
int estimateRank ( ArrayWrapper<Type, n> &mat )
{
	Type w [ n ];
	ArrayWrapper<Type, n> u;
	int i;

/*
	mat [ 0 ] [ 0 ] = a [ 0 ] * a [ 0 ];
	mat [ 0 ] [ 1 ] = a [ 1 ] * a [ 0 ];
	mat [ 0 ] [ 2 ] = a [ 2 ] * a [ 0 ];
	mat [ 1 ] [ 0 ] = a [ 1 ] * a [ 0 ];
	mat [ 1 ] [ 1 ] = a [ 1 ] * a [ 1 ] + a [ 4 ] * a [ 4 ];
	mat [ 1 ] [ 2 ] = a [ 1 ] * a [ 2 ] + a [ 4 ] * a [ 5 ];
	mat [ 2 ] [ 0 ] = a [ 2 ] * a [ 0 ];
	mat [ 2 ] [ 1 ] = a [ 1 ] * a [ 2 ] + a [ 5 ] * a [ 4 ];
	mat [ 2 ] [ 2 ] = a [ 2 ] * a [ 2 ] + a [ 5 ] * a [ 5 ] + a [ 7 ] * a [ 7 ];
*/
	jacobi<Type, n> ( mat, w, u );

	if ( w [ 0 ] == 0.0f )
	{
		return 0;
	}
	else
	{
		for ( i = 1; i < n; i++ )
		{
			if ( w [ i ] < 0.1f )
			{
				return i;
			}
		}

		return n;
	}
}

template <class Type, int n>
void nullspace ( ArrayWrapper<Type, n> &mat, Type *rvalue )
{
	Type w [ n ];
	ArrayWrapper<Type, n> u;
	int i;

	jacobi<Type, n> ( mat, w, u );

	// return last vector
	for ( i = 0; i < n; i++ )
	{
		rvalue [ i ] = u.data [ n - 1 ] [ i ];
	}
}

template <class Type, int n>
void matInverse ( ArrayWrapper<Type, n> &mat, ArrayWrapper<Type, n> &rvalue, Type tolerance )
{
	// there is an implicit assumption that mat is symmetric and real
	// U and V in the SVD will then be the same matrix whose rows are the eigenvectors of mat
	// W will just be the eigenvalues of mat
	Type w [ n ];
	ArrayWrapper<Type, n> u;
	int i, j, k;

	jacobi<Type, n> ( mat, w, u );

/*	if ( fabs ( w [ 0 ] ) < tolerance )
	{
		for ( i = 0; i < n; i++ )
		{
			for ( j = 0; j < n; j++ )
			{
				rvalue.data [ i ] [ j ] = 0;
			}
		}
//		printf ( "error: largest eigenvalue is 0!\n" );
		return;
	}
	else */
	{
		for ( i = 1; i < n; i++ )
		{
//			if ( w [ i ] / w [ 0 ] < 0.02f ) // / w [ 0 ] < TOLERANCE )
			if ( fabs ( w [ i ] / w [ 0 ] ) < tolerance )
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
			rvalue.data [ i ] [ j ] = 0;
			for ( k = 0; k < n; k++ )
			{
				rvalue.data [ i ] [ j ] += w [ k ] * u.data [ k ] [ i ] * u.data [ k ] [ j ];
			}
		}
	}
}

template <class Type, int n>
void calcPoint ( Type midpoint[], Type rvalue[], ArrayWrapper<Type, n> &mat, Type *b )
{
	Type newB [ n ];
	ArrayWrapper<Type, n> inv;
	int i, j;

	for ( i = 0; i < n; i++ )
	{
		newB [ i ] = b [ i ];
		for ( j = 0; j < n; j++ )
		{
			newB [ i ] -= mat.data [ i ] [ j ] * midpoint [ j ];
		}
	}

	::matInverse<Type, n> ( mat, inv, 0.0001 );

	for ( i = 0; i < n; i++ )
	{
		rvalue [ i ] = midpoint [ i ];
		for ( j = 0; j < n; j++ )
		{
			rvalue [ i ] += (Type)( inv.data [ j ] [ i ] * newB [ j ] );
		}
	}
}


template <class Type, int n>
QEFQR<Type, n>::QEFQR( Type *eqn )
{
	makeQEF ( eqn );
}

template <class Type, int n>
void QEFQR<Type, n>::combineSelf ( Type *eqn )
{
	int i, j, index;
	Type a, b, mag, temp;
	int rvalueN = 0;
	Type scale = 0;

	eqn [ n ] *= -1;

	for ( i = 0; i < n + 1; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;

		a = data [ index ];
		b = eqn [ i ];

		if ( fabs ( a ) > 0 || fabs ( b ) > 0 )
		{
			rvalueN = i + 1;

			mag = (Type)sqrt ( a * a + b * b );
			a /= mag;
			b /= mag;

			for ( j = 0; j < n + 1 - i; j++ )
			{
				temp = a * data [ index + j ] + b * eqn [ j + i ];
				eqn [ j + i ] = b * data [ index + j ] - a * eqn [ j + i ];
				data [ index + j ] = temp;
			}
		}
	}

	eqn [ n ] *= -1;

	// print the equation
	//fprintf(trace_file, "qef data = ");
	//for ( i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
	//{
	//	fprintf(trace_file, "%I64x ", *(__int64*)&data [ i ]);
	//	//fprintf(trace_file, "%g ", data [ i ]);
	//}
	//fprintf(trace_file, "\n");
}

template <class Type, int n>
void QEFQR<Type, n>::combineSelf ( QEF<Type, n> &qef )
{
	int i, j, k, index1, index2;
	Type a, b, mag, temp;
	QEFQR<Type, n> tempQEF = *((QEFQR<Type, n> *)&qef);
//	Type *tempQEF = qef.data; 

	for ( i = 0; i < n + 1; i++ ) // && i < qef.num; i++ )
	{
		int row = 0;
		index1 = ( ( 2 * n + 3 - i ) * i ) / 2;

		index2 = 0;

		for ( j = 0; j < i; j++ )
		{
			a = data [ index2 + j - row ];

			if ( fabs ( a ) > 0 )
			{
				row++;
				index2 = ( ( 2 * n + 3 - row ) * row ) / 2;
			}
		}

		for ( ; j < n + 1; j++ ) // && row < num + 1; j++ )
		{
			a = data [ index2 + j - row ];
			b = tempQEF.data [ index1 + j - i ];

			if ( fabs ( a ) > 0 || fabs ( b ) > 0 )
			{
				mag = (Type)sqrt ( a * a + b * b );
				a /= mag;
				b /= mag;

				for ( k = j; k < n + 1; k++ )
				{
					temp = a * data [ index2 + k - row ] + b * tempQEF.data [ index1 + k - i ];
					tempQEF.data [ index1 + k - i ] = b * data [ index2 + k - row ] - a * tempQEF.data [ index1 + k - i ];
					data [ index2 + k - row ] = temp;
				}

				row++;
/*
				if ( row > num )
				{
					num = row;
				}
*/
				index2 = ( ( 2 * n + 3 - row ) * row ) / 2;
			}
		}
	}
/*
	// destroys argument QEF
	for ( i = 0; i < n + 1; i++ )
	{
		index1 = ( ( 2 * n + 3 - i ) * i ) / 2;

		for ( j = 0; j <= i; j++ )
		{
			index2 = ( ( 2 * n + 3 - j ) * j ) / 2;

			a = data [ index1 ];
			b = tempQEF [ index2 + i - j ];

			if ( fabs ( a ) > 0.000001f || fabs ( b ) > 0.000001f )
			{
				mag = (Type)sqrt ( a * a + b * b );
				a /= mag;
				b /= mag;

				for ( k = i; k < n + 1; k++ )
				{
					temp = a * data [ index1 + k - i ] + b * tempQEF [ index2 + k - j ];
					tempQEF [ index2 + k - j ] = b * data [ index1 + k - i ] - a * tempQEF [ index2 + k - j ];
					data [ index1 + k - i ] = temp;
				}
			}
		}
	}
*/
//	num += qef.num;
}

template <class Type, int n>
void QEFQR<Type, n>::combineSelf ( Type **eqs, int num )
{
	int i;

	for ( i = 0; i < num; i++ )
	{
		combineSelf ( eqs [ i ] );
	}

//	this->num += num;
}

template <class Type, int n>
void QEFQR<Type, n>::makeQEF ( Type *eqn )
{
	int i;

	for ( i = 0; i < n + 1; i++ )
	{
		data [ i ] = eqn [ i ];
	}

	for ( ; i < ( ( n + 1 ) * ( n + 2 ) ) / 2; i++ )
	{
		data [ i ] = 0;
	}

//	num = 1;
}

template <class Type, int n>
Type QEFQR<Type, n>::calcError ( Type point[] ) const
{
	Type error;
	Type tempError;
	int i, j;
	int index;

	error = data [ ( n + 1 ) * ( n + 2 ) / 2 - 1 ] * data [ ( n + 1 ) * ( n + 2 ) / 2 - 1 ];

	for ( i = 0; i < n ; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;

		tempError = 0;
		for ( j = i; j < n; j++ )
		{
			tempError += data [ index + j - i ] * point [ j ];
		}
		tempError -= data [ index + n - i ];

		error += tempError * tempError;
	}

	return error;
}

template <class Type, int n>
Type QEFQR<Type, n>::minimizerError ( Type point[] ) const
{
//	return data [ ( n + 1 ) * ( n + 2 ) / 2 - 1 ] * data [ ( n + 1 ) * ( n + 2 ) / 2 - 1 ];
	return calcError ( point );
}

template <class Type, int n>
int QEFQR<Type, n>::estimateRank ( void ) const
{
	int i, j, k;
	int index;
	ArrayWrapper<Type, n> mat;

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			mat.data [ i ] [ j ] = 0;
		}
	}

	for ( k = 0; k < n; k++ )
	{
		index = ( ( 2 * n + 3 - k ) * k ) / 2;
		for ( i = k; i < n; i++ )
		{
			for ( j = k; j < n; j++ )
			{
				mat.data [ i ] [ j ] += data [ index + i - k ] * data [ index + j - k ];
			}
		}
	}

	return ::estimateRank<Type, n> ( mat );
}

template <class Type, int n>
Type *QEFQR<Type, n>::nullspace ( void ) const
{
	int i, j, k;
	int index;
	ArrayWrapper<Type, n> mat;
	Type *rvalue = new Type [ n ];

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			mat.data [ i ] [ j ] = 0;
		}
	}

	for ( k = 0; k < n; k++ )
	{
		index = ( ( 2 * n + 3 - k ) * k ) / 2;
		for ( i = k; i < n; i++ )
		{
			for ( j = k; j < n; j++ )
			{
				mat.data [ i ] [ j ] += data [ index + i - k ] * data [ index + j - k ];
			}
		}
	}

	::nullspace<Type, n> ( mat, rvalue );

	return rvalue;
}

template <class Type, int n>
void QEFQR<Type, n>::calcPoint ( Type midpoint[], Type rvalue[] ) const
{
	//// always collapse to the midpoint (debug)
	for ( int i = 0; i < n; i++ )
	{
		rvalue [ i ] = midpoint [ i ];
	}
	return;

	ArrayWrapper<Type, n> a;
	Type b [ n ];
	int i, j, k;
	int index;

	for ( i = 0; i < n; i++ )
	{
		b [ i ] = 0;
		for ( j = 0; j < n; j++ )
		{
			a.data [ i ] [ j ] = 0;
		}
	}

	for ( k = 0; k < n; k++ )
	{
		index = ( ( 2 * n + 3 - k ) * k ) / 2;
		for ( i = k; i < n; i++ )
		{
			for ( j = k; j < n; j++ )
			{
				a.data [ i ] [ j ] += data [ index + i - k ] * data [ index + j - k ];
			}
		}
	}

	for ( k = 0; k < n; k++ )
	{
		for ( i = 0; i <= k; i++ )
		{
			index = ( ( 2 * n + 3 - i ) * i ) / 2;
			b [ k ] += data [ index + k - i ] * data [ index + n - i ];
		}
	}

	
	::calcPoint<Type, n> ( midpoint, rvalue, a, b );
}

//template <class Type, int n>
//void QEFQR<Type, n>::calcPointOnBox ( Type midpoint[], Type rvalue[], const MathBoundingBox3 &box ) const
//{
//	calcPoint ( midpoint, rvalue );
//  for ( int i = 0; i < 3; i++ )
//  {
//    if ( rvalue [ i ] < box.minCorner [ i ] )
//    {
////      printf ( "clamping %.2f,%.2f,%.2f inside box %.1f,%.1f,%.1f ->%.1f,%.1f,%.1f\n", rvalue [ 0 ], rvalue [ 1 ], rvalue [ 2 ],
////        box.minCorner.x, box.minCorner.y, box.minCorner.z, box.maxCorner.x, box.maxCorner.y, box.maxCorner.z );
//      rvalue [ i ] = box.minCorner [ i ];
//    }
//    else if ( rvalue [ i ] > box.maxCorner [ i ] )
//    {
////      printf ( "clamping %.2f,%.2f,%.2f inside box %.1f,%.1f,%.1f ->%.1f,%.1f,%.1f\n", rvalue [ 0 ], rvalue [ 1 ], rvalue [ 2 ],
////        box.minCorner.x, box.minCorner.y, box.minCorner.z, box.maxCorner.x, box.maxCorner.y, box.maxCorner.z );
//      rvalue [ i ] = box.maxCorner [ i ];
//    }
//  }
//}

template <class Type, int n>
QEFNormal<Type, n>::QEFNormal( Type *eqn )
{
	makeQEF ( eqn );
}

template <class Type, int n>
void QEFNormal<Type, n>::combineSelf ( Type *eqn )
{
	int i, j;
	int index;

//	eqn [ n ] *= -1;
	index = 0;
	for ( i = 0; i < n + 1; i++ )
	{
		for ( j = i; j < n + 1; j++ )
		{
			data [ index ] += eqn [ i ] * eqn [ j ];
			index++;
		}
	}
//	eqn [ n ] *= -1;
}

template <class Type, int n>
void QEFNormal<Type, n>::combineSelf ( QEF<Type, n> &qef )
{
	int i;

	for ( i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
	{
		data [ i ] += qef.data [ i ];
	}
}

template <class Type, int n>
void QEFNormal<Type, n>::combineSelf ( Type **eqs, int num )
{
	int i;

	for ( i = 0; i < num; i++ )
	{
		combineSelf ( eqs [ i ] );
	}
}

template <class Type, int n>
void QEFNormal<Type, n>::makeQEF ( Type *eqn )
{
	int i, j;
	int index;

//	eqn [ n ] *= -1;
	index = 0;
	for ( i = 0; i < n + 1; i++ )
	{
		for ( j = i; j < n + 1; j++ )
		{
			data [ index ] = eqn [ i ] * eqn [ j ];
			index++;
		}
	}
//	eqn [ n ] *= -1;
}

template <class Type, int n>
int QEFNormal<Type, n>::estimateRank ( void ) const
{
	ArrayWrapper<Type, n> mat;
	int i, j;
	int index;

	for ( i = 0; i < n; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;
		for ( j = i; j < n; j++ )
		{
			mat.data [ i ] [ j ] = data [ index + j - i ];
			mat.data [ j ] [ i ] = mat.data [ i ] [ j ];
		}
	}

	return ::estimateRank<Type, n> ( mat );
}

template <class Type, int n>
Type *QEFNormal<Type, n>::nullspace ( void ) const
{
	ArrayWrapper<Type, n> mat;
	int i, j;
	int index;
	Type *rvalue = new Type [ n ];

	for ( i = 0; i < n; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;
		for ( j = i; j < n; j++ )
		{
			mat.data [ i ] [ j ] = data [ index + j - i ];
			mat.data [ j ] [ i ] = mat.data [ i ] [ j ];
		}
	}

	::nullspace<Type, n> ( mat, rvalue );

	return rvalue;
}

template <class Type, int n>
void QEFNormal<Type, n>::calcPoint ( Type midpoint[], Type rvalue[] ) const
{
	ArrayWrapper<Type, n> a;
	Type b [ n ];
	int i, j;
	int index;

	for ( i = 0; i < n; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;
		for ( j = i; j < n; j++ )
		{
			a.data [ i ] [ j ] = data [ index + j - i ];
			a.data [ j ] [ i ] = a.data [ i ] [ j ];
		}

		b [ i ] = -data [ index + n - i ];
	}

	::calcPoint<Type, n> ( midpoint, rvalue, a, b );
}

template <class Type, int n>
Type QEFNormal<Type, n>::calcError ( Type point[] ) const
{
	Type error = 0;
	int i, j;
	int index;

	// calc x^T A^T A x
	for ( i = 0; i < n; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;

		for ( j = i; j < n; j++ )
		{
			error += point [ i ] * data [ index + j - i ] * point [ j ];
		}

		for ( j = 0; j < i; j++ )
		{
			index = ( ( 2 * n + 3 - j ) * j ) / 2;

			error += point [ i ] * data [ index + i - j ] * point [ j ];
		}
	}

	// x^T A^T B
	for ( i = 0; i < n; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2 + n - i;
		error = error + 2 * ( data [ index ] * point [ i ] );
	}

	// B^T B
	error += data [ ( n + 1 ) * ( n + 2 ) / 2 - 1 ];

	return error;
}

template <class Type, int n>
Type QEFNormal<Type, n>::minimizerError ( Type point[] ) const
{
	Type error = 0;
	int i;
	int index;

	error = data [ ( n + 1 ) * ( n + 2 ) / 2 - 1 ];
	// x^T A^T B
	for ( i = 0; i < n; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2 + n - i;
		error += ( data [ index ] * point [ i ] );
	}

	return error;
}

//template <class Type, int n>
//void QEFNormal<Type, n>::calcPointOnBox ( Type midpoint[], Type rvalue[], const MathBoundingBox3 &box ) const
//{
//	calcPoint ( midpoint, rvalue );
//}

template <class Type, int n>
void QEFCenter<Type, n>::combineSelf ( Type *eqn ) 
{ 
	data [ 0 ] += 1;
	data [ 1 ] += eqn [ 0 ];
	data [ 2 ] += eqn [ 1 ];
	data [ 3 ] += eqn [ 2 ];
	data [ 4 ] += eqn [ 0 ] * eqn [ 0 ] + eqn [ 1 ] * eqn [ 1 ] + eqn [ 2 ] * eqn [ 2 ];
}

template <class Type, int n>
void QEFCenter<Type, n>::combineSelf ( QEF<Type, n> &qef )
{
	for ( int i = 0; i < qef.mydata.size ( ); i++ )
	{
		mydata.push_back ( qef.mydata [ i ] );
	}
	data [ 0 ] += qef.data [ 0 ];
	data [ 1 ] += qef.data [ 1 ];
	data [ 2 ] += qef.data [ 2 ];
	data [ 3 ] += qef.data [ 3 ];
	data [ 4 ] += qef.data [ 4 ];
}

template <class Type, int n>
void QEFCenter<Type, n>::combineSelf ( Type **eqs, int num )
{
	data [ 0 ] = num;
	data [ 1 ] = eqs [ 0 ] [ 0 ];
	data [ 2 ] = eqs [ 0 ] [ 1 ];
	data [ 3 ] = eqs [ 0 ] [ 2 ];
	data [ 4 ] = eqs [ 0 ] [ 0 ] * eqs [ 0 ] [ 0 ] + eqs [ 0 ] [ 1 ] * eqs [ 0 ] [ 1 ] + eqs [ 0 ] [ 2 ] * eqs [ 0 ] [ 2 ];
	for ( int i = 1; i < num; i++ )
	{
		data [ 1 ] += eqs [ i ] [ 0 ];
		data [ 2 ] += eqs [ i ] [ 1 ];
		data [ 3 ] += eqs [ i ] [ 2 ];
		data [ 4 ] += eqs [ i ] [ 0 ] * eqs [ i ] [ 0 ] + eqs [ i ] [ 1 ] * eqs [ i ] [ 1 ] + eqs [ i ] [ 2 ] * eqs [ i ] [ 2 ];
	}
}

template <class Type, int n>
void QEFCenter<Type, n>::makeQEF ( Type *eqn )
{
	vect3d t;
	t[0] = eqn[0];
	t[1] = eqn[1];
	t[2] = eqn[2];


	mydata.push_back ( t );

	data [ 0 ] = 1;
	data [ 1 ] = eqn [ 0 ];
	data [ 2 ] = eqn [ 1 ];
	data [ 3 ] = eqn [ 2 ];
	data [ 4 ] = eqn [ 0 ] * eqn [ 0 ] + eqn [ 1 ] * eqn [ 1 ] + eqn [ 2 ] * eqn [ 2 ];
}

template <class Type, int n>
void QEFCenter<Type, n>::calcPoint ( Type midpoint[], Type rvalue[] ) const
{
	rvalue[0] = midpoint[0];
	rvalue[1] = midpoint[1];
	rvalue[2] = midpoint[2];
	/*rvalue [ 0 ] = data [ 1 ] / data [ 0 ];
	rvalue [ 1 ] = data [ 2 ] / data [ 0 ];
	rvalue [ 2 ] = data [ 3 ] / data [ 0 ];*/
}

template <class Type, int n>
Type QEFCenter<Type, n>::calcError ( Type point[] ) const
{
	return data [ 0 ] * ( point [ 0 ] * point [ 0 ] + point [ 1 ] * point [ 1 ] + point [ 2 ] * point [ 2 ] ) - 2 * ( point [ 0 ] * data [ 1 ] + point [ 1 ] * data [ 2 ] + point [ 2 ] * data [ 3 ] ) + data [ 4 ];
}

template <class Type, int n>
Type QEFCenter<Type, n>::minimizerError ( Type point[] ) const
{
	return data [ 4 ] - ( data [ 1 ] * data [ 1 ] + data [ 2 ] * data [ 2 ] + data [ 3 ] * data [ 3 ] ) / data [ 0 ];
}

template <class Type, int n>
QEFCenter<Type, n>::QEFCenter( Type *eqn )
{
	makeQEF ( eqn );
}