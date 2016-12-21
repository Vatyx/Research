#ifndef POINT_GRID_H
#define POINT_GRID_H

#include "vect.h"
#include "ResizeArray.h"

//#define OLD_POINT_GRID

template <int sz>
class PointGrid
{
	vect2d minCorner, extents;
#ifdef OLD_POINT_GRID
	ResizeArray<int> grid [ sz ] [ sz ];
#else
	int head[sz][sz];
	int tail[sz][sz];
	int *data;
	int *ptr;
	int maxSize, nextFree;
#endif

public:
#ifdef OLD_POINT_GRID
	PointGrid ( ) {}
#else
	PointGrid()
	{
		maxSize = sz * sz * 2;
		data = new int[maxSize];
		ptr = new int[maxSize];
		memset(head, -1, sizeof (int)* sz*sz);
		memset(tail, -1, sizeof (int)* sz*sz);
		nextFree = 0;
	}

	~PointGrid()
	{
		delete[] data;
		delete[] ptr;
	}
#endif

	void setBoundingBox ( const vect2d &minP, const vect2d &maxP )
	{
		minCorner = minP;
		extents = ( maxP - minP ) * ( 1 + 0.0000001 );
	}

	void clear ( void )
	{
#ifdef OLD_POINT_GRID
		for ( int i = 0; i < sz; i++ )
		{
			for ( int j = 0; j < sz; j++ )
			{
				grid [ i ] [ j ].clear ( );
			}
		}
#else
		memset( head, -1, sizeof(int) * sz * sz );
		memset( tail, -1, sizeof(int) * sz * sz );
		nextFree = 0;
#endif
	}

	void addPoint(const vect2d &p, int index)
	{
		vect2d temp = ((p - minCorner) * sz) / extents;
		int x = (int)temp[0];
		int y = (int)temp[1];

#ifdef OLD_POINT_GRID
		grid [ x ] [ y ].insert ( index );
#else
		if ( nextFree > maxSize )
		{
			printf("major Error\n"); 
		}
		else
		{
			data[nextFree] = index;
			ptr[nextFree] = -1;
			if ( tail [ x ] [ y ] != -1 )
			{
				ptr[tail[x][y]] = nextFree;
			}
			else
			{
				head[x][y] = nextFree;
			}
			tail[x][y] = nextFree;
			nextFree++;
		}
#endif
	}

	void findPoints ( ResizeArray<int> &ans, vect2d v1, vect2d v2, double dist )
	{
		int minX, minY, maxX, maxY;

		minX = (int)(sz * ( min ( v1[0], v2[0] ) - dist - minCorner[0] ) / extents[0]);
		minY = (int)(sz * ( min ( v1[1], v2[1] ) - dist - minCorner[1] ) / extents[1]);
		maxX = (int)(sz * ( max ( v1[0], v2[0] ) + dist - minCorner[0] ) / extents[0]);
		maxY = (int)(sz * ( max ( v1[1], v2[1] ) + dist - minCorner[1] ) / extents[1]);
		if ( minX < 0 )
		{
			minX = 0;
		}
		if ( minY < 0 )
		{
			minY = 0;
		}
		if ( maxX >= sz )
		{
			maxX = sz - 1;
		}
		if ( maxY >= sz )
		{
			maxY = sz - 1;
		}

		// just uses bounding box... could make this tighter if necessary...
		for ( int i = minX; i <= maxX; i++ )
		{
			for ( int j = minY; j <= maxY; j++ )
			{
#ifdef OLD_POINT_GRID
				for ( int k = 0; k < grid [ i ] [ j ].length ( ); k++ )
//				for (int k = grid[i][j].length() - 1; k >= 0; k--)
				{
					ans.insert ( grid [ i ] [ j ].data [ k ] );
				}
#else
				int curr = head [ i ] [ j ];
				while (curr != -1)
				{
					ans.insert ( data [ curr ] );
					curr = ptr [ curr ];
				}
#endif
			}
		}
	}
};

#endif