#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "Matrix.h"
#include <iostream>
#include "LList.h"

using namespace std;

namespace MathMatrix
{

template <class type>
class SparseData
{
public:
	int c;
	type data;
};

template <class type>
class ColSelector 
{
protected:
	LList<SparseData<type> > **row;

public:
	ColSelector ( LList<SparseData<type> > **sel ) : row ( sel ) { }
	ColSelector ( ) : row ( NULL ) { }

	type &operator [] ( int index )
	{
		LList<SparseData<type> > *ptr = *row;

		if ( ptr == NULL || ptr->getData ( ).c > index )
		{
//			return 0;
			SparseData<type> newData;
			newData.c = index;
			newData.data = 0;

			LList<SparseData<type> > *newElement = new LList<SparseData<type> > ( newData, ptr );
			*row = newElement;
			return newElement->getData ( ).data;
		}

		while ( ptr->getNext ( ) != NULL && ptr->getNext ( )->getData ( ).c < index )
		{
			ptr = ptr->getNext ( );
		}

		if ( ptr->getData ( ).c == index )
		{
			return ptr->getData ( ).data;
		}

		if ( ptr->getNext ( ) == NULL || ptr->getNext ( )->getData ( ).c > index )
		{
//			return 0;
			SparseData<type> newData;
			newData.c = index;
			newData.data = 0;

			LList<SparseData<type> > *newElement = new LList<SparseData<type> > ( newData, ptr->getNext ( ) );
			ptr->setNext ( newElement );
			return newElement->getData ( ).data;
		}

		return ptr->getNext ( )->getData ( ).data;
	}

	type getValue ( int index ) const
	{
		LList<SparseData<type> > *ptr = *row;

		if ( ptr == NULL || ptr->getData ( ).c > index )
		{
			return 0;
		}

		while ( ptr->getNext ( ) != NULL && ptr->getNext ( )->getData ( ).c < index )
		{
			ptr = ptr->getNext ( );
		}

		if ( ptr->getData ( ).c == index )
		{
			return ptr->getData ( ).data;
		}

		if ( ptr->getNext ( ) == NULL || ptr->getNext ( )->getData ( ).c > index )
		{
			return 0;
		}

		return ptr->getNext ( )->getData ( ).data;
	}
};

template <class type>
class SparseMatrix
{
protected:
	LList<SparseData<type> > **rowData;
	int r, c;

public:
	SparseMatrix ( ) : r ( 0 ), c ( 0 ), rowData ( NULL ) { }

	SparseMatrix ( int rows, int columns )
	{
		int i;

		r = rows;
		c = columns;

//		data = NULL;
		rowData = new LList<SparseData<type> > *[ rows ];
		for ( i = 0; i < r; i++ )
		{
			rowData [ i ] = NULL;
		}
	}

	SparseMatrix ( const SparseMatrix<type> &mat )
	{
		int i;

		r = mat.getNumRows ( );
		c = mat.getNumCols ( );

//		data = NULL;
		rowData = new LList<SparseData<type> > *[ r ];

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = mat.rowData [ i ];
			rowData [ i ] = NULL;
			LList<SparseData<type> > *end = NULL;

			while ( ptr != NULL )
			{
				if ( end == NULL )
				{
					rowData [ i ] = new LList<SparseData<type> > ( ptr->getData ( ) );
					end = rowData [ i ];
				} 
				else
				{
					end->setNext ( new LList<SparseData<type> > ( ptr->getData ( ) ) );
					end = end->getNext ( );
				}

				ptr = ptr->getNext ( );
			}
		}
	}

	LList<SparseData<type> > **getSparseData ( void )
	{
		return rowData;
	}

	virtual ~SparseMatrix ( )
	{
		int i;

		if ( rowData )
		{
			for ( i = 0; i < r; i++ )
			{
				if ( rowData [ i ] )
				{
					rowData [ i ]->deleteList ( );
				}
			}
		}
		delete[] rowData;
	}

	int getNumRows ( void ) const
	{
		return r;
	}

	int getNumCols ( void ) const
	{
		return c;
	}

	virtual SparseMatrix<type>& zero ( void )
	{
		int i;

		for ( i = 0; i < r; i++ )
		{
			if ( rowData [ i ] )
			{
				rowData [ i ]->deleteList ( );
				rowData [ i ] = NULL;
			}
		}

		return *this;
	}


	virtual SparseMatrix<type>& identity ( void )
	{
		int i;

		if ( r == c )
		{
			zero ( );

			for ( i = 0; i < r; i++ )
			{
				SparseData<type> dat;

				dat.c = i;
				dat.data = 1;

				rowData [ i ] = new LList<SparseData<type> > ( dat );
			}
		}

		return *this;
	}

	virtual int operator = ( Matrix<type> &mat )
	{
		int i, j;

		if ( rowData )
		{
			zero ( );
		}

		r = mat.getNumRows ( );
		c = mat.getNumCols ( );

		rowData = new LList<SparseData<type> > *[ r ];
		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *end;
			rowData [ i ] = 0;

			for ( j = 0; j < c; j++ )
			{
				SparseData<type> dat;

				dat.c = j;
				dat.data = mat [ i ] [ j ];

				if ( j == 0)
				{
					rowData [ i ] = new LList<SparseData<type> > ( dat );
					end = rowData [ i ];
				}
				else
				{
					end->setNext ( new LList<SparseData<type> > ( dat ) );
					end = end->getNext ( );
				}
			}
		}

		return 0;
	}

	virtual int operator = ( const SparseMatrix<type> &mat )
	{
		int i;

		if ( rowData )
		{
			zero ( );
		}

		r = mat.r;
		c = mat.c;

		rowData = new LList<SparseData<type> > *[ r ];
		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = mat.rowData [ i ];
			LList<SparseData<type> > *end = NULL;
			rowData [ i ] = NULL;

			while ( ptr != NULL )
			{
				if ( end == NULL )
				{
					rowData [ i ] = new LList<SparseData<type> > ( ptr->getData ( ) );
					end = rowData [ i ];
				}
				else
				{
					end->setNext ( new LList<SparseData<type> > ( ptr->getData ( ) ) );
					end = end->getNext ( );
				}
				ptr = ptr->getNext ( );
			}
		}

		return 0;
	}

	virtual Matrix<type> operator * ( Matrix<type> &mat )
	{
		assert ( c == mat.getNumRows ( ) );
		int i, j;
		Matrix<type> rvalue ( r, mat.getNumCols ( ) );
		rvalue.zero ( );

		for ( i = 0; i < r; i++ )
		{
			for ( j = 0; j < mat.getNumCols ( ); j++ )
			{
				LList<SparseData<type> > *ptr = rowData [ i ];

				while ( ptr != NULL )
				{
					rvalue [ i ] [ j ] += ptr->getData ( ).data * mat [ ptr->getData ( ).c ] [ j ];
					ptr = ptr->getNext ( );
				}
			}
		}

		return rvalue;
	}

	virtual SparseMatrix<type> operator * ( SparseMatrix<type> &mat )
	{
		assert ( c == mat.getNumRows ( ) );
		int i;
		SparseMatrix<type> rvalue ( r, mat.getNumCols ( ) );

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *myPtr = rowData [ i ];

			while ( myPtr != NULL )
			{
				LList<SparseData<type> > *matPtr = mat.rowData [ myPtr->getData ( ).c ];
				while ( matPtr != NULL )
				{
					rvalue [ i ] [ matPtr->getData ( ).c ] += myPtr->getData ( ).data * matPtr->getData ( ).data;
					matPtr = matPtr->getNext ( );
				}
				myPtr = myPtr->getNext ( );
			}
		}

		return rvalue;
	}

	virtual SparseMatrix<type> operator * ( type val )
	{
		SparseMatrix<type> rvalue ( *this );

		rvalue *= val;

		return rvalue;
	}

	virtual void operator *= ( Matrix<type> &mat )
	{
		assert ( c == mat.getNumRows ( ) );

		*this = *this * mat;
	}

	virtual void operator *= ( type val )
	{
		int i;

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];
			while ( ptr != NULL )
			{
				ptr->getData ( ).data *= val;
				ptr = ptr->getNext ( );
			}
		}
	}

	virtual Matrix<type> operator + ( Matrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );
		int i;
		Matrix<type> rvalue ( mat );

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				rvalue [ i ] [ ptr->getData ( ).c ] += ptr->getData ( ).data;
				ptr = ptr->getNext ( );
			}
		}

		return rvalue;
	}

	virtual SparseMatrix<type> operator + ( SparseMatrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );
		int i;
		SparseMatrix<type> rvalue ( mat );

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				rvalue [ i ] [ ptr->getData ( ).c ] += ptr->getData ( ).data;
				ptr = ptr->getNext ( );
			}
		}

		return rvalue;
	}

	virtual SparseMatrix<type> operator + ( type val )
	{
		SparseMatrix<type> rvalue ( *this );

		rvalue += val;

		return rvalue;
	}

	virtual void operator += ( Matrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );

		*this = *this + mat;
	}

	virtual void operator += ( SparseMatrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );

		*this = *this + mat;
	}

	virtual void operator += ( type val )
	{
		int i;

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				ptr->getData ( ).data += val;
				ptr = ptr->getNext ( );
			}
		}
	}

	virtual Matrix<type> operator - ( Matrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );
		int i;
		Matrix<type> rvalue ( mat );

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				rvalue [ i ] [ ptr->getData ( ).c ] -= ptr->getData ( ).data;
				ptr = ptr->getNext ( );
			}
		}

		return rvalue;
	}

	virtual SparseMatrix<type> operator - ( SparseMatrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );
		int i;
		SparseMatrix<type> rvalue ( mat );

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				rvalue [ i ] [ ptr->getData ( ).c ] -= ptr->getData ( ).data;
				ptr = ptr->getNext ( );
			}
		}

		return rvalue;
	}

	virtual SparseMatrix<type> operator - ( type val )
	{
		SparseMatrix<type> rvalue ( *this );

		rvalue -= val;

		return rvalue;
	}

	virtual void operator -= ( Matrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );

		*this = *this - mat;
	}

	virtual void operator -= ( SparseMatrix<type> &mat )
	{
		assert ( r == mat.getNumRows ( ) && c == mat.getNumCols ( ) );

		*this = *this - mat;
	}

	virtual void operator -= ( type val )
	{
		int i;

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				ptr->getData ( ).data -= val;
				ptr = ptr->getNext ( );
			}
		}
	}

	virtual ColSelector<type> operator [] ( int index )
	{
		assert ( index >= 0 && index < r );

		return ColSelector<type> ( &rowData [ index ] );
	}

	virtual SparseMatrix<type> transpose ( void )
	{
		SparseMatrix<type> rvalue ( c, r );
		int i;

		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *myPtr = rowData [ i ];

			while ( myPtr != NULL )
			{
				SparseData<type> newData;
				LList<SparseData<type> > *rvaluePtr = rvalue.rowData [ myPtr->getData ( ).c ];
				LList<SparseData<type> > *newElement;

				newData.data = myPtr->getData ( ).data;
				newData.c = i;

				newElement = new LList<SparseData<type> > ( newData );
				if ( rvaluePtr == NULL || rvaluePtr->getData ( ).c > i )
				{
					newElement->setNext ( rvalue.rowData [ myPtr->getData ( ).c ] );

					rvalue.rowData [ myPtr->getData ( ).c ] = newElement;
				}
				else
				{
					while ( rvaluePtr->getNext ( ) != NULL && rvaluePtr->getNext ( )->getData ( ).c < i )
					{
						rvaluePtr = rvaluePtr->getNext ( );
					}
					newElement->setNext ( rvaluePtr->getNext ( ) );
					rvaluePtr->setNext ( newElement );
				}
				myPtr = myPtr->getNext ( );
			}
		}

		return rvalue;
	}

	Matrix<type> toFullMatrix ( void )
	{
		Matrix<type> rvalue ( r, c );
		int i;

		rvalue.zero ( );
		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				rvalue [ i ] [ ptr->getData ( ).c ] = ptr->getData ( ).data;
				ptr = ptr->getNext ( );
			}
		}

		return rvalue;
	}

	static SparseMatrix<type> toSparseMatrix ( Matrix<type> &mat )
	{
		SparseMatrix<type> rvalue ( mat.getNumRows ( ), mat.getNumCols ( ) );
		int i, j;

		for ( i = 0; i < mat.getNumRows ( ); i++ )
		{
			LList<SparseData<type> > *ptr = rvalue.rowData [ i ];

			for ( j = 0; j < mat.getNumCols ( ); j++ )
			{
				if ( mat [ i ] [ j ] != 0 )
				{
					SparseData<type> newData;
					newData.c = j;
					newData.data = mat [ i ] [ j ];

					if ( ptr == NULL )
					{
						rvalue.rowData [ i ] = new LList<SparseData<type> > ( newData );
						ptr = rvalue.rowData [ i ];
					}
					else
					{
						ptr->setNext ( new LList<SparseData<type> > ( newData ) );
						ptr = ptr->getNext ( );
					}
				}
			}
		}

		return rvalue;
	}

	void conjugateGradients ( Matrix<type> &b, Matrix<type> &x, int maxIter = 100 )
	{
		assert ( getNumRows ( ) == getNumCols ( ) );
		int i, j, k;

		Matrix<type> r = b - (*this) * x;
		Matrix<type> d = r;
		type *deltaNew, deltaOld, beta, alpha;
		Matrix<type> q;

		deltaNew = new type [ x.getNumCols ( ) ];
		type error = 0;
		for ( i = 0; i < x.getNumCols ( ); i++ )
		{
			deltaNew [ i ] = 0;
			for ( j = 0; j < getNumRows ( ); j++ )
			{
				deltaNew [ i ] += r [ j ] [ i ] * r [ j ] [ i ];
			}
			error += deltaNew [ i ];
		}
		if ( !_finite ( error ) )
		{
			printf ( "starting error\n" );
			return;
		}

		if ( error < 0.000001 )
		{
//			printf ( "CG return early because error too small\n" );
			return;
		}

		i = 0;
//		maxIter = maxIter > getNumRows ( ) ? maxIter : getNumRows ( );
		while ( i < maxIter )
		{
			q = (*this) * d;
			for ( j = 0; j < x.getNumCols ( ); j++ )
			{
				type temp = 0;
				for ( k = 0; k < getNumRows ( ); k++ )
				{
					temp += d [ k ] [ j ] * q [ k ] [ j ];
				}
				type alpha = deltaNew [ j ] / temp;
				if ( !_finite ( alpha ) )
				{
					printf ( "error\n" );
					return;
					continue;
				}
				for ( k = 0; k < getNumRows ( ); k++ )
				{
					x [ k ] [ j ] += d [ k ] [ j ] * alpha;
					assert ( _finite ( x [ k ] [ j ] ) );
				}
				if ( i % 50 != 0 )
				{
					for ( k = 0; k < getNumRows ( ); k++ )
					{
						r [ k ] [ j ] -= q [ k ] [ j ] * alpha;
					}
				}
			}
			if ( i % 50 == 0 )
			{
//				printf ( "%d of %d\n", i, maxIter );
				r = b - (*this) * x;
			}
			for ( j = 0; j < x.getNumCols ( ); j++ )
			{
				deltaOld = deltaNew [ j ];
				deltaNew [ j ] = 0;
				for ( k = 0; k < getNumRows ( ); k++ )
				{
					deltaNew [ j ] += r [ k ] [ j ] * r [ k ] [ j ];
				}
				beta = deltaNew [ j ] / deltaOld;
				for ( k = 0; k < getNumRows ( ); k++ )
				{
					d [ k ] [ j ] = r [ k ] [ j ] + d [ k ] [ j ] * beta;
				}
			}
			i++;
		}
	}

		void conjugateGradients ( type *b, type *x, int maxIter = 100, type epsilon = 0.000001 )
	{
		assert ( getNumRows ( ) == getNumCols ( ) );
		int i, j, k;

		type *r = new type [ getNumRows ( ) ];
		type *d = new type [ getNumRows ( ) ];
		type *q = new type [ getNumRows ( ) ];
//		Matrix<type> r = b - (*this) * x;
		for ( j = 0; j < getNumRows ( ); j++ )
		{
			r [ j ] = b [ j ];
			LList<SparseData<type> > *myPtr = rowData [ j ];

			while ( myPtr != NULL )
			{
				r [ j ] -= myPtr->getData ( ).data * x [ myPtr->getData ( ).c ];
				myPtr = myPtr->getNext ( );
			}
		}
		

		for ( i = 0; i < getNumRows ( ); i++ )
		{
			d [ i ] = r [ i ];
		}
		type deltaNew, deltaOld, beta, alpha;

		deltaNew = 0;
		for ( i = 0; i < getNumRows ( ); i++ )
		{
			deltaNew += r [ i ] * r [ i ];
		} 
		type delta0 = deltaNew;

		if ( deltaNew < epsilon )
		{
//			printf ( "CG return early because error too small\n" );
			return;
		}

		printf ( "iterating\n" );
		i = 0;
		while ( deltaNew > epsilon && i < maxIter )
		{
//			q = (*this) * d;
			for ( j = 0; j < getNumRows ( ); j++ )
			{
				q [ j ] = 0;
				LList<SparseData<type> > *myPtr = rowData [ j ];

				while ( myPtr != NULL )
				{
					q [ j ] += myPtr->getData ( ).data * d [ myPtr->getData ( ).c ];
					myPtr = myPtr->getNext ( );
				}
			}

			type temp = 0;
			for ( j = 0; j < getNumRows ( ); j++ )
			{
				temp += d [ j ] * q [ j ];
			}

			type alpha = deltaNew / temp;
			assert ( _finite ( alpha ) );
			for ( k = 0; k < getNumRows ( ); k++ )
			{
				x [ k ] += d [ k ] * alpha;
//				assert ( _finite ( x [ k ] ) );
			}
			if ( i % 50 != 0 )
			{
				for ( k = 0; k < getNumRows ( ); k++ )
				{
					r [ k ] -= q [ k ] * alpha;
				}
			} else
			{
//				r = b - (*this) * x;
				for ( j = 0; j < getNumRows ( ); j++ )
				{
					r [ j ] = b [ j ];
					LList<SparseData<type> > *myPtr = rowData [ j ];

					while ( myPtr != NULL )
					{
						r [ j ] -= myPtr->getData ( ).data * x [ myPtr->getData ( ).c ];
						myPtr = myPtr->getNext ( );
					}
				}
			}

			deltaOld = deltaNew;
			deltaNew = 0;
			for ( j = 0; j < getNumRows ( ); j++ )
			{
				deltaNew += r [ j ] * r [ j ];
			}
			beta = deltaNew / deltaOld;
			for ( k = 0; k < getNumRows ( ); k++ )
			{
				d [ k ] = r [ k ] + d [ k ] * beta;
			}
			i++;
		}
		printf ( "%d of %d\n", i, maxIter );
		printf("error %.10g \n", deltaNew);

		delete[] r;
		delete[] d;
		delete[] q;
	}


	// I stopped here
	virtual void solve ( Matrix<type> &mat )
	{
		assert ( false );
		Matrix<type> temp = toFullMatrix ( );

		temp.solve ( mat );
	}

	virtual Matrix<type> pseudoInverse ( void )
	{
//		assert ( false );
		Matrix<type> temp = toFullMatrix ( );

		return temp.pseudoInverse ( );
	}

	virtual void outputMathematica ( char *filename )
	{
		FILE *fptr = fopen ( filename, "wt" );
		bool first = true;
		int i;

		fprintf ( fptr, "SparseArray[{" );
		
		for ( i = 0; i < r; i++ )
		{
			LList<SparseData<type> > *ptr = rowData [ i ];

			while ( ptr != NULL )
			{
				if ( !first )
				{
					fprintf ( fptr, "," );
				}
				fprintf ( fptr, "{%d,%d}->%f", i + 1, ptr->getData ( ).c + 1, ptr->getData ( ).data );
				first = false;

				ptr = ptr->getNext ( );
			}
		}
		fprintf ( fptr, "},{%d,%d}]", getNumRows ( ), getNumCols ( ) );
		fclose ( fptr );
	}
};

template <class type>
ostream &operator<< ( ostream& os, SparseMatrix<type>& mat )
{
	int i, j;

	for ( i = 0; i < mat.getNumRows ( ); i++ )
	{
		for ( j = 0; j < mat.getNumCols ( ); j++ )
		{
			os << mat [ i ].getValue ( j ) << " ";
		}
		os << endl;
	}

	return os;
}

}
#endif