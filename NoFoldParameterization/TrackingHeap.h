#ifndef TRACKING_HEAP_STRUCTURE_H
#define TRACKING_HEAP_STRUCTURE_H

/**
 * Templated heap data structure.  This should be used for a priority
 * queue.  Inserting and deleting the highest priority element are
 * always log(n).  The tree it keeps is always perfectly balanced.
 * The comparing function defines a minimizing heap
 *
 * @author Scott Schaefer
 */

#include "AdaptiveArray.h"
#include <stdio.h>
class TrackingHeap
{
public:
	// users must override this class
	class TrackingHeapData
	{
	public:
		int index;
	};

	// internal representation for an item
	struct HeapElement
	{
		TrackingHeapData *data;
		float sortKey;
	};

private:
	// stores the data for the heap
	AdaptiveArray<HeapElement> heapData;

	/**
	 * Walks up the tree from the given index 
	 * and enforces the heap property
	 */
	void up ( int index )
	{
		// test if we're the root (because we're done then)
		if ( index > 0 )
		{
			// not the root
			int p = parent ( index );

			if ( compare ( heapData.getElement ( p ).sortKey, heapData.getElement ( index ).sortKey ) )
			{
				swap ( p, index );
				up ( p );
			}
		}
	}

	/**
	 * Enforces the heap property walking down 
	 * the tree from the given index.
	 */
	void heapify ( int index )
	{
		int left = leftChild ( index );
		int right = rightChild ( index );
		int largest;

		if ( left < heapData.getSize ( ) && 
			compare ( heapData.getElement ( index ).sortKey, heapData.getElement ( left ).sortKey ) )
		{
			largest = left;
		} 
		else
		{
			largest = index;
		}

		if ( right < heapData.getSize ( ) &&
			compare ( heapData.getElement ( largest ).sortKey, heapData.getElement ( right ).sortKey ) )
		{
			largest = right;
		}

		if ( largest != index )
		{
			swap ( index, largest );
			heapify ( largest );
		}
	}

	/**
	 * Gives the index of the parent of the given node
	 */
	int parent ( int index )
	{
		return ( index - 1 ) / 2;
	}

	/**
	 * Gives the index of the left child of the given node
	 */
	int leftChild ( int index )
	{
		return 2 * index + 1;
	}

	/**
	 * Gives the index of the right child of the given node
	 */
	int rightChild ( int index )
	{
		return 2 * ( index + 1 );
	}

	/**
	 * Swaps two elements in the list
	 */
	void swap ( int i, int j )
	{
		HeapElement tmp;

		tmp = heapData.getElement ( i );
		heapData.getElement ( i ) = heapData.getElement ( j );
		heapData.getElement ( j ) = tmp;
		heapData.getElement ( i ).data->index = i;
		heapData.getElement ( j ).data->index = j;
	}

protected:
	/**
	 * Function to override to determine how the heap sorts items
	 */
	virtual bool compare ( float value1, float value2 )
	{
		return value1 > value2;
	}

public:
	/**
	 * Constructor
	 */
	TrackingHeap ( void ) {}

	/**
	 * Destructor
	 */
	virtual ~TrackingHeap ( void ) {}

	/**
	 * Inserts a piece of data into the heap associated 
	 * with the given sort key.
	 * 
	 * @param data the data to insert into the heap
	 * @param sortKey the priority of the data in the heap
	 */
	void insert ( TrackingHeapData *data, float sortKey )
	{
		HeapElement item;

		item.data = data;
		item.data->index = heapData.getSize ( );
		item.sortKey = sortKey;

		heapData.insertElement ( item );

		up ( heapData.getSize ( ) - 1 );
	}

	/**
	 * Serves the highest priority item out of the heap.
	 * 
	 * @param item the place to store the served item
	 *
	 * @return true if the operation succeeded or 
	 *         false if the heap was empty
	 */
	TrackingHeapData *serveNext ( void )
	{
		TrackingHeapData *item = NULL;

		if ( heapData.getSize ( ) == 0 )
		{
			return false;
		}

		item = heapData.getElement ( 0 ).data;

		heapData.removeElement ( 0 );
		if ( heapData.getSize ( ) > 0 )
		{
			heapData.getElement ( 0 ).data->index = 0;

			heapify ( 0 );
		}

		item->index = -2;
		return item;
	}

	/**
	 * @return whether or not the heap contains any items
	 */
	bool isEmpty ( void )
	{
		return heapData.getSize ( ) == 0;
	}

	/**
	 * Note that this removal is O(log(n))... typical heap is O(n)
	 *
	 * @param item the item to remove from the heap
	 *
	 * @return true if the item was found or false otherwise
	 */
	bool remove ( TrackingHeapData *item )
	{
		//assert ( item->index < heapData.getSize ( ) );

		if ( item->index >=0 )
		{
			heapData.removeElement ( item->index );
			if ( item->index < heapData.getSize ( ) ) // could have been the last element!!!
			{
				heapData.getElement ( item->index ).data->index = item->index;
				heapify ( item->index );
			}
			item->index = -2;
		}

		return true;
	}

	/**
	 * @return the priority associated with the next item to be
	 *         served from the heap
	 */
	float peekNextPriority ( void )
	{
		if ( isEmpty ( ) )
		{
			return 0;
		}

		return heapData.getElement ( 0 ).sortKey;
	}

	TrackingHeapData *getElement ( int ind )
	{
		return heapData.getElement ( ind ).data;
	}

	/**
	 * @return the number of items in the heap
	 */
	int getNumItems ( void )
	{
		return heapData.getSize ( );
	}

	/**
	 * This function is really for debug purposes.  It prints a
	 * list of all of the priorities in the heap in the order in
	 * which they appear in the array built internally by the heap
	 */
	void printPriorities ( void )
	{
		int i;

		for ( i = 0; i < heapData.getSize ( ); i++ )
		{
			printf ( "%f ", heapData.getElement ( i ).sortKey );
		}
		printf ( "\n" );
	}

	/**
	 * Clears the heap of all items stored in it.  The heap is
	 * empty after this operation.
	 */
	void clear ( void )
	{
		heapData.clear ( );
	}
};

#endif