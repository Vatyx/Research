/********************************************************************
	created:	2001/05/11
	created:	11:5:2001   18:11
	filename: 	e:\my documents\cvs\final_aggression\util\adaptivearray.h
	file path:	e:\my documents\cvs\final_aggression\util
	file base:	adaptivearray
	file ext:	h
	author:		Frank Losasso

	Copyright(c) 2001  Frank Losasso, All Rights Reserved.
	
	purpose:	An array that grows when needed
*********************************************************************/

#ifndef ADAPTIVEARRAY_H
#define ADAPTIVEARRAY_H


//#include "../stdafx.h"


template <class Type>
class AdaptiveArray
{
private:
	Type* array;
	int size;
	int numElements;
	
	void doubleSize()
	{
		Type* temp = new Type[size*2];

		// copy the array content
		for (int i = 0 ; i < size; i++)
		{
			temp[i] = array[i];
		}

		delete[] array;
		array = temp;
		size *= 2;
	}

public:
	AdaptiveArray ( const AdaptiveArray<Type> &a )
	{
		int i;

		size = a.size;
		numElements = a.numElements;
		array = new Type [ size ];

		for ( i = 0; i < numElements; i++ )
		{
			array [ i ] = a.array [ i ];
		}
	}

	~AdaptiveArray()
	{
		delete[] array;
	}

	int operator = ( const AdaptiveArray<Type> &a )
	{
		copy ( &a );

		return 1;
	}

	void copy ( const AdaptiveArray<Type> *a )
	{
		int i;

		delete[] array;

		size = a->size;
		numElements = a->numElements;
		array = new Type [ size ];

		for ( i = 0; i < numElements; i++ )
		{
			array [ i ] = a->array [ i ];
		}
	}

	/**
	 * constructor which sets the initial size
	 *
	 * @param initialSize  the initial size of the array
	 */
	AdaptiveArray(int initialSize = 1)
	{
		array = new Type[initialSize];

		size = initialSize;

		numElements = 0;
	}


	/**
	 * function to insert element into the array
	 *
	 * @param element  the element to insert into the array
	 */
	int insertElement(Type element)
	{
		if (numElements >= size)
			doubleSize();

		array[numElements] = element;
		numElements++;

		return numElements - 1;
	}

	/**
	 * function that will return the index of the next empty element in the array
	 * NOTE: this function should only be used when Type is NOT a pointer, as the
	 *       pointer will be trash
	 *
	 * DO NOT USE THIS FUNCTION UNLESS YOU UNDERSTAND WHAT IT DOES... ONLY FOR USE WITH REAL TYPES!!!
	 *
	 * @returns the index of the new element that can be used
	 */
	int getNewElement()
	{
		if (numElements >= size)
			doubleSize();

		numElements++;

		return numElements-1;
	}

	/**
	 * removes an element from the array
	 *
	 * @param index  the index of the element to remove
	 */
	void removeElement(int index)
	{
		assert(index < numElements && index >=0);

		array[index] = array[numElements-1];
		numElements--;
	}

	void removeAndShift ( int index )
	{
		assert ( index < numElements && index >= 0 );
		int i;

		numElements--;
		for ( i = index; i < numElements; i++ )
		{
			array [ i ] = array [ i + 1 ];
		}
	}

	void removeData ( Type element )
	{
		int i;

		for ( i = 0; i < numElements; i++ )
		{
			if ( element == array [ i ] )
			{
				removeElement ( i );
				return;
			}
		}
	}

	void removeAndShiftData ( Type element )
	{
		int i;

		for ( i = 0; i < numElements; i++ )
		{
			if ( element == array [ i ] )
			{
				removeAndShift ( i );
				return;
			}
		}
	}

	/**
	 * Returns the tightly packet array which is keeped internally
	 *
	 * @param arraySize  the size of the array that will be returned
	 * @return the array of hashable objects
	 */
	Type* getArray(int& arraySize)
	{
		arraySize = numElements;
		return array;
	}

	/**
	 * Returns an element from the array
	 *
	 * @param index   the index into the array
	 * @return the element at index in the array
	 */
	Type& getElement(int index)
	{
		assert(index >= 0);
		assert(index < numElements);
		return array[index];
	}

	int memberOf ( Type &x ) const
	{
		int i;

		for ( i = 0; i < numElements; i++ )
		{
			if ( array [ i ] == x )
			{
				return i;
			}
		}

		return -1;
	}

	/**
	 * returns the number of elements in the array
	 *
	 * @return the number of elements
	 */
	int getSize() const
	{
		return numElements;
	}

	/**
	 * clears the contents of the array
	 */
	void clear()
	{
		numElements = 0;
	}

	/**
	 * Reverses the elements in this list
	 */
	void reverse ( void )
	{
		int i;
		Type swap;


		for ( i = 0; i < numElements / 2; i++ )
		{
			swap = array [ i ];
			array [ i ] = array [ numElements - 1 - i ];
			array [ numElements - 1 - i ] = swap;
		}
	}
};

#endif