#ifndef QEF_H
#define QEF_H

#include <vector>
#include "vect.h"

#include <stdio.h>
//#include "BoundingBox3.h"

template <class Type, int n>
struct ArrayWrapper
{
	Type data [ n ] [ n ];
};

template <class Type, int n>
class QEF
{
public:
	Type data [ ( n + 1 ) * ( n + 2 ) / 2 ];
	std::vector<vect3d> mydata;
//	int num;

public:
	QEF ( void ) { }//num = 0; }

	void zero ( void )
	{
		int i;

		for ( i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
		{
			data [ i ] = 0;
		}
	}

	virtual void combineSelf ( Type *eqn ) = 0;

	virtual void combineSelf ( QEF<Type, n> &qef ) = 0;

	virtual void combineSelf ( Type **eqs, int num ) = 0;

	virtual void makeQEF ( Type *eqn ) = 0;

	virtual int estimateRank ( void ) const = 0;

	virtual void calcPoint ( Type midpoint[], Type rvalue[] ) const = 0;

	virtual Type calcError ( Type point[] ) const = 0;

	virtual Type minimizerError ( Type point[] ) const = 0;

	//virtual void calcPointOnBox ( Type midpoint[], Type rvalue[], const MathBoundingBox3 &box ) const = 0;

	virtual Type *nullspace ( void ) const = 0;

	virtual void print ( void )
	{
		int i;

//    printf ( "num is %d\n", num );
		for ( i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
		{
			printf ( "%f ", data [ i ] );
		}
		printf ( "\n" );
	}

	virtual void print_mathmatica ( void )
	{
		int i;

//    printf ( "num is %d\n", num );
		printf("{");
		for ( i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
		{
			printf ( "%.10g", data [ i ] );
			if( i < ( n + 1 ) * ( n + 2 ) / 2 - 1)
				printf(", ");
		}
		printf ( "}\n" );
	}
};



template <class Type, int n>
class QEFQR : public QEF<Type, n>
{
public:
	QEFQR ( void ) : QEF<Type, n> ( ) { }

	QEFQR ( Type *eqn );

	virtual void combineSelf ( Type *eqn );

	virtual void combineSelf ( QEF<Type, n> &qef );

	virtual void combineSelf ( Type **eqs, int num );

	virtual void makeQEF ( Type *eqn );

	virtual int estimateRank ( void ) const;

	virtual void calcPoint ( Type midpoint[], Type rvalue[] ) const;

	virtual Type calcError ( Type point[] ) const;

	virtual Type minimizerError ( Type point[] ) const;

	virtual Type *nullspace ( void ) const;

	//virtual void calcPointOnBox ( Type midpoint[], Type rvalue[], const MathBoundingBox3 &box ) const;
};

template <class Type, int n>
class QEFNormal : public QEF<Type, n>
{
public:
	QEFNormal ( void ) : QEF<Type, n> ( ) {}

	QEFNormal ( Type *eqn );

	virtual void combineSelf ( Type *eqn );

	virtual void combineSelf ( QEF<Type, n> &qef );

	virtual void combineSelf ( Type **eqs, int num );

	virtual void makeQEF ( Type *eqn );

	virtual int estimateRank ( void ) const;

	virtual void calcPoint ( Type midpoint[], Type rvalue[] ) const;

	virtual Type calcError ( Type point[] ) const;

	virtual Type minimizerError ( Type point[] ) const;

	virtual Type *nullspace ( void ) const;

	//virtual void calcPointOnBox ( Type midpoint[], Type rvalue[], const MathBoundingBox3 &box ) const;
};

template <class Type, int n>
class QEFCenter : public QEF<Type, n>
{

public:
	//
	QEFCenter ( void )  : QEF<Type, n>() { }
	QEFCenter ( Type *eqn );
	virtual void combineSelf ( Type *eqn );

	virtual void combineSelf ( QEF<Type, n> &qef );

	virtual void combineSelf ( Type **eqs, int num );

	virtual void makeQEF ( Type *eqn );

	virtual int estimateRank ( void ) const { return 0; }

	virtual void calcPoint ( Type midpoint[], Type rvalue[] ) const;

	virtual Type calcError ( Type point[] ) const;

	virtual Type minimizerError ( Type point[] ) const;

	//virtual void calcPointOnBox ( Type midpoint[], Type rvalue[], const MathBoundingBox3 &box ) const = 0;

	virtual Type *nullspace ( void ) const { return NULL; }
};

template <class Type, int n>
void matInverse ( ArrayWrapper<Type, n> &mat, ArrayWrapper<Type, n> &rvalue, Type tolerance = 0.000001 );

#endif