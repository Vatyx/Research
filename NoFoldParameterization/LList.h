/********************************************************************
	created:	2001/10/31
	created:	31:10:2001   21:20
	filename: 	C:\Documents and Settings\Administrator\Desktop\openGL\util\LList.h
	file path:	C:\Documents and Settings\Administrator\Desktop\openGL\util
	file base:	LList
	file ext:	h
	author:		Frank Losasso
	
	purpose:	A linked list
*********************************************************************/

#ifndef LList_h
#define LList_h

//#include "../stdafx.h"
#include <stdio.h>


template <class Type>
class LList
{
private:
	LList<Type>*	next;
	LList<Type>*	previous;
	Type			data;

public:
  LList ( void )
  {
    next = NULL;
    previous = NULL;
  }

	LList(Type f)
	{
		next = NULL;
		previous = NULL;
		data = f;
	}

	LList(Type f, LList<Type>* nextElt)
	{
		next = nextElt;
		data = f;
	}

	~LList()
	{
	}

	// changed to iterative instead of recursive to prevent stack overflow -SS
	void deleteList()
	{
		LList<Type> *ptr;

		ptr = this;

		while ( ptr->next != NULL )
		{
			ptr->next->previous = ptr;
			ptr = ptr->next;
		}

		while ( ptr != this )
		{
			ptr = ptr->previous;
			delete ptr->next;
		}

		// must do this the backackwards way because this must be deleted last
		// if method was static, this wouldn't be a problem
		delete this;
	}

	void setData(Type& d)
	{
		data = d;
	}

	Type& getData()
	{
		return data;
	}

	void setNext(LList<Type>* a)
	{
		next = a;
	}

	void setPrevious(LList<Type>* prev)
	{
		previous = prev;
	}

	LList<Type>* getNext()
	{
		return next;
	}

	LList<Type>* getPrevious()
	{
		return previous;
	}
};







template <class Type>
class SLinkedList
{
private:
	LList<Type>*	current;

	LList<Type>*	first;

	LList<Type>*	last;

public:
	SLinkedList()
	{
		first = NULL;
		last = NULL;
		current = NULL;
	}

	void deleteAll()
	{
		first->deleteList();
	}

	void addElement(Type& data)
	{
		LList<Type>* lList = new LList<Type>(data);

		if (last == NULL)
		{
			last = lList;
			first = lList;
			current = lList;
		}
		last->setNext(lList);
		lList->setPrevious(last);
		last = lList;
		current = last;
	}

	Type getNext()
	{
		if (current == NULL)
			return NULL;
		if (current->getNext() == NULL)
		{
			return current->getData();
		}
		current = current->getNext();
		return current->getData();
	}

	Type getPrevious()
	{
		if (current == NULL)
			return NULL;
		if (current->getPrevious() == NULL)
		{
			return current->getData();
		}
		current = current->getPrevious();
		return current->getData();
	}
};






#endif
