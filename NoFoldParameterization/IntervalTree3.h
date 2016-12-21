#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include "ResizeArray.h"
#include "BoundingBox2.h"
#include <assert.h>

//#define OLD_GRID

template <int size>
class IntervalTree
{
	vect2d minCorner, extents;
#ifdef OLD_GRID
	ResizeArray<int> grid [ size ] [ size ];
	ResizeArray<bool> used;
#else
	int head[size][size];
	int *data;
	int *ptr;
	bool *used;
	int maxSize, nextFree, nextUsed; 
#endif

public:
#ifdef OLD_GRID
	IntervalTree ( ) {}
#else
	IntervalTree() 
	{
		maxSize = size * size * 40;
		data = new int[maxSize];
		ptr = new int[maxSize];
		used = new bool[size * size * 2]; // better not go over.. maybe need a parameter for the constructor
		memset(head, -1, sizeof (int)* size*size);
		memset(used, 0, sizeof(bool)* size * size * 2);
		nextFree = 0;
		nextUsed = 0;
	}

	~IntervalTree()
	{
		delete[] used;
	}
#endif

	void setBoundingBox ( const vect2d &minP, const vect2d &maxP )
	{
		minCorner = minP;
		extents = ( maxP - minP ) * ( 1 + 0.0000001 );
	}
	
	void clear ( void )
	{
#ifdef OLD_GRID
		for ( int i = 0; i < size; i++ )
		{
			for ( int j = 0; j < size; j++ )
			{
				grid [ i ] [ j ].clear ( );
			}
		}
		used.clear ( );
#else
		memset(head, -1, sizeof(int)* size * size);
		nextFree = 0;
		nextUsed = 0;
#endif
	}

	void insert ( BoundingBox &box )
	{
#ifdef OLD_GRID
		int index = used.length ( );
		used.insert ( false );
#else
		int index = nextUsed;
		used [ nextUsed++ ] = false;
#endif
		vect2d minBox = ( ( box.minCorner - minCorner ) * size ) / extents;
		vect2d maxBox = ( ( box.maxCorner - minCorner ) * size ) / extents;

		int minX = (int)minBox[0];
		int minY = (int)minBox[1];
		int maxX = (int)maxBox[0];
		int maxY = (int)maxBox[1];

		for ( int i = minX; i <= maxX; i++ )
		{
			for ( int j = minY; j <= maxY; j++ )
			{
#ifdef OLD_GRID
				grid [ i ] [ j ].insert ( index );
#else
				if (nextFree > maxSize)
				{
					printf("major Error\n");
				}
				else
				{
					data[nextFree] = index;
					ptr[nextFree] = head[i][j];
					head[i][j] = nextFree;
					nextFree++;
				}
#endif
			}
		}
	}

	void insertLineSegment ( vect2d &start, vect2d &end )
	{
#ifdef OLD_GRID
		int index = used.length ( );
		used.insert ( false );
#else
		int index = nextUsed;
		used[nextUsed++] = false;
#endif

		vect2d startP = ( ( start - minCorner ) * size ) / extents;
		vect2d endP = ( ( end - minCorner ) * size ) / extents;
		double nextXt, nextYt, deltaXt, deltaYt;

		int currX = (int)startP[0];
		int currY = (int)startP[1];
		int xIncr, yIncr;

		// potential optimization
		if ( (int)(endP [ 0 ]) == currX )
		{
			for ( int i = min ( currY, (int)(endP[1]) ); i <= max ( currY, (int)(endP[1]) ); i++ )
			{
#ifdef OLD_GRID
				grid [ currX ] [ i ].insert ( index );
#else
				if (nextFree > maxSize)
				{
					printf("major Error\n");
				}
				else
				{
					data[nextFree] = index;
					ptr[nextFree] = head[currX][i];
					head[currX][i] = nextFree;
					nextFree++;
				}
#endif
			}
			return;
		}
		else if ( (int)(endP [ 1 ] ) == currY )
		{
			for ( int i = min ( currX, (int)(endP[0]) ); i <= max ( currX, (int)(endP[0]) ); i++ )
			{
#ifdef OLD_GRID
				grid [ i ] [ currY ].insert ( index );
#else
				if (nextFree > maxSize)
				{
					printf("major Error\n");
				}
				else
				{
					data[nextFree] = index;
					ptr[nextFree] = head[i][currX];
					head[i][currX] = nextFree;
					nextFree++;
				}
#endif
			}
			return;
		}


		xIncr = ( endP [ 0 ] > startP [ 0 ] ) - ( startP [ 0 ] < endP [ 0 ] ); // sign bit (+1, -1, 0)
		if ( xIncr != 0 )
		{
			nextXt = ( currX + xIncr - startP[0] ) / ( endP [ 0 ] - startP [ 0 ] );
			deltaXt = ( xIncr ) / ( endP [ 0 ] - startP [ 0 ] );
		}
		else
		{
			nextXt = FLT_MAX;
		}
		yIncr = ( endP [ 1 ] > startP [ 1 ] ) - ( startP [ 1 ] < endP [ 1 ] ); // sign bit (+1, -1, 0)
		if ( yIncr != 0 )
		{
			nextYt = ( currY + yIncr - startP[1] ) / ( endP [ 1 ] - startP [ 1 ] );
			deltaYt = ( yIncr ) / ( endP [ 1 ] - startP [ 1 ] );
		}
		else
		{
			nextYt = FLT_MAX;
		}

		do
		{
#ifdef OLD_GRID
			grid [ currX ] [ currY ].insert ( index );
#else
			if (nextFree > maxSize)
			{
				printf("major Error\n");
			}
			else
			{
				data[nextFree] = index;
				ptr[nextFree] = head[currX][currY];
				head[currX][currY] = nextFree;
				nextFree++;
			}
#endif

			if ( nextXt < nextYt )
			{
				nextXt += deltaXt;
				currX += xIncr;
			}
			else
			{
				nextYt += deltaYt;
				currY += yIncr;
			}
		} while ( nextXt < 1 || nextYt < 1 );
	}

	void findIntersect ( BoundingBox &box, ResizeArray<int> &ans )
	{
		int minX, minY, maxX, maxY;

		minX = (int)(size * ( box.minCorner[0] - minCorner[0] ) / extents[0]);
		minY = (int)(size * ( box.minCorner[1] - minCorner[1] ) / extents[1]);
		maxX = (int)(size * ( box.maxCorner[0] - minCorner[0] ) / extents[0]);
		maxY = (int)(size * ( box.maxCorner[1] - minCorner[1] ) / extents[1]);
		// BEGIN: SHOULD NOT BE NECESSARY
/*		if ( minX < 0 )
		{
			minX = 0;
		}
		if ( minY < 0 )
		{
			minY = 0;
		}
		if ( maxX >= size )
		{
			maxX = size - 1;
		}
		if ( maxY >= size )
		{
			maxY = size - 1;
		}*/
		// END: SHOULD NOT BE NECESSARY

		for ( int i = minX; i <= maxX; i++ )
		{
			for (int j = minY; j <= maxY; j++)
			{
#ifdef OLD_GRID
				for (int k = 0; k < grid[i][j].length(); k++)
				{
					// eliminate duplicates
					if (!(used.data[grid[i][j].data[k]]))
					{
						used.data[grid[i][j].data[k]] = true;

						ans.insert(grid[i][j].data[k]);
					}
				}
#else
				int curr = head [ i ] [ j ];
				while (curr != -1)
				{
					if ( !used [ data [ curr ] ] )
					{
						used [ data[curr] ] = true;
						ans.insert ( data [ curr ] );
					}
					curr = ptr [ curr ];
				}
#endif
			}
		}

		for ( int i = 0; i < ans.length ( ); i++ )
		{
#ifdef OLD_GRID
			used.data [ ans.data [ i ] ] = false;
#else
			used [ ans.data [ i ] ] = false;
#endif
		}
	}

#ifdef FALSE // this code appears to be unused
	void findIntersectLineSegment ( vect2d &start, vect2d &end, ResizeArray<int> &ans )
	{
		vect2d startP = ( ( start - minCorner ) * size ) / extents;
		vect2d endP = ( ( end - minCorner ) * size ) / extents;
		double nextXt, nextYt, deltaXt, deltaYt;

		int currX = (int)startP[0];
		int currY = (int)startP[1];
		int xIncr, yIncr;

		if ( (int)(endP [ 0 ]) == currX )
		{
			for ( int i = min ( currY, (int)(endP[1]) ); i <= max ( currY, (int)(endP[1]) ); i++ )
			{
				for ( int k = 0; k < grid [ currX ] [ i ].length ( ); k++ )
				{
					// eliminate duplicates
					if ( !(used.data [ grid [ currX ] [ i ].data [ k ] ]) )
					{
						used.data [ grid [ currX ] [ i ].data [ k ] ] = true;

						ans.insert ( grid [ currX ] [ i ].data [ k ] );
					}
				}
			}
		}
		else if ( (int)(endP [ 1 ] ) == currY )
		{
			for ( int i = min ( currX, (int)(endP[0]) ); i <= max ( currX, (int)(endP[0]) ); i++ )
			{
				for ( int k = 0; k < grid [ i ] [ currY ].length ( ); k++ )
				{
					// eliminate duplicates
					if ( !(used.data [ grid [ i ] [ currY ].data [ k ] ]) )
					{
						used.data [ grid [ i ] [ currY ].data [ k ] ] = true;

						ans.insert ( grid [ i ] [ currY ].data [ k ] );
					}
				}
			}
		}
		else
		{
			xIncr = ( endP [ 0 ] > startP [ 0 ] ) - ( startP [ 0 ] < endP [ 0 ] ); // sign bit (+1, -1, 0)
			if ( xIncr != 0 )
			{
				nextXt = ( currX + xIncr - startP[0] ) / ( endP [ 0 ] - startP [ 0 ] );
				deltaXt = ( xIncr ) / ( endP [ 0 ] - startP [ 0 ] );
			}
			else
			{
				nextXt = FLT_MAX;
			}
			yIncr = ( endP [ 1 ] > startP [ 1 ] ) - ( startP [ 1 ] < endP [ 1 ] ); // sign bit (+1, -1, 0)
			if ( yIncr != 0 )
			{
				nextYt = ( currY + yIncr - startP[1] ) / ( endP [ 1 ] - startP [ 1 ] );
				deltaYt = ( yIncr ) / ( endP [ 1 ] - startP [ 1 ] );
			}
			else
			{
				nextYt = FLT_MAX;
			}

			do
			{
				for ( int k = 0; k < grid [ currX ] [ currY ].length ( ); k++ )
				{
					// eliminate duplicates
					if ( !(used.data [ grid [ currX ] [ currY ].data [ k ] ]) )
					{
						used.data [ grid [ currX ] [ currY ].data [ k ] ] = true;

						ans.insert ( grid [ currX ] [ currY ].data [ k ] );
					}
				}

				if ( nextXt < nextYt )
				{
					nextXt += deltaXt;
					currX += xIncr;
				}
				else
				{
					nextYt += deltaYt;
					currY += yIncr;
				}
			} while ( nextXt < 1 || nextYt < 1 );
		}

		for ( int i = 0; i < ans.length ( ); i++ )
		{
			used.data [ ans.data [ i ] ] = false;
		}
	}
#endif
};

#endif