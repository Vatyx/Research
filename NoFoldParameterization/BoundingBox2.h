#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "vect.h"

class BoundingBox
{
public:
	vect2d minCorner, maxCorner;

	BoundingBox ( void ) { }
	BoundingBox ( vect2d minC, vect2d maxC ) : minCorner ( minC ), maxCorner ( maxC ) { }

	bool intersect (/* const*/ BoundingBox &box ) /*const*/
	{
		return !( minCorner[1] > box.maxCorner[1] || maxCorner[1] < box.minCorner[1] ||
			 minCorner[0] > box.maxCorner[0] || maxCorner[0] < box.minCorner[0] );
	}
	
	bool intersectLineSegment ( /*const*/ vect2d &p1, /*const*/ vect2d &p2 ) /*const*/
	{
		// if calling from IntervalTree 
		if ( minCorner[1] > max ( p1[1], p2[1] ) || maxCorner[1] < min ( p1[1], p2[1] ) ||
			 minCorner[0] > max ( p1[0], p2[0] ) || maxCorner[0] < min ( p1[0], p2[0] ) )
		{
			return false;
		}
		vect2d ortho;
		ortho[0]= p2[1] - p1[1]; ortho[1] = p1[0] - p2[0];
		vect2d box [ 3 ];
		box [ 0 ] = maxCorner;
		box [ 1 ][0] = minCorner[0];
		box [ 1 ][1] = maxCorner[1];
		box [ 2 ][0] = maxCorner[0];
		box [ 2 ] [1] = minCorner[1];

		double proj = ( minCorner - p1 ).dot ( ortho );
		if ( proj > 0 )
		{
			for ( int i = 0; i < 3; i++ )
			{
				proj = ( box [ i ] - p1 ).dot ( ortho );
				if ( proj <= 0 )
				{
					return true;
				}
			}
			return false;
		}
		else if ( proj < 0 )
		{
			for ( int i = 0; i < 3; i++ )
			{
				proj = ( box [ i ] - p1 ).dot ( ortho );
				if ( proj >= 0 )
				{
					return true;
				}
			}
			return false;
		}
		return true;
	}
};

#endif