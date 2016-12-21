#ifndef RESIZE_ARRAY_H
#define RESIZE_ARRAY_H

template <class Type>
class ResizeArray
{
public:
	Type *data;
	int size;
	int count;

	ResizeArray ( int defaultSize = 4 )
	{
		size = defaultSize;
		data = new Type [ size ];
		count = 0;
	}

	~ResizeArray ( void )
	{
		delete[] data;
	}

	void insert ( Type elem )
	{
		if ( count >= size )
		{
			Type *newData = new Type [ 2 * size ];
			for ( int i = 0; i < size; i++ )
			{
				newData [ i ] = data [ i ];
			}
			delete[] data;
			data = newData;
			size *= 2;
		}
		data [ count ] = elem;
		count++;
	}

	void clear ( void )
	{
		count = 0;
	}

	int length ( void ) const
	{
		return count;
	}
};

#endif