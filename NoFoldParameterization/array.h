#pragma once
#include <assert.h>

template <class T>
struct Array
{
	Array(){data = 0; s = 0; c = 0;}
	Array(const Array<T> &a){data = 0; s = 0; c = 0; copy(a);}
	~Array(){clear();}

	T *data;
	int s, c; // size and capacity of the array

	void operator=(const Array<T> &a) {copy(a);}

	void copy(const Array<T> &a)
	{
		if (a.s)
		{
			resize(a.s);
			for (int i = 0; i < s; i++)
				data[i] = a.data[i];
		}
	}

	int size() const {return s;}
	bool empty(){return s == 0;}

	void resize(int s_new)
	{
		reserve(s_new);
		s = s_new;
	}

	void reserve(int c_new)
	{
		if (c_new < s)
			c_new = s;
		
		if (c_new > c)
		{
			T *d = new T[c_new];

			for (int i = 0; i < s; i++)
				d[i] = data[i];
			delete[] data;

			data = d;
			c = c_new;
		}
	}

	T &operator[](int i) const
	{
		assert(c >= s);
		assert(i >= 0 && i < s);
		return data[i];
	}
	
	T &get(int i) const
	{
		assert(c >= s);
		assert(i >= 0 && i < s);
		return data[i];
	}

	void clear()
	{
		delete[] data;
		data = 0;
		s = 0;
		c = 0;
	}

	void push_back(const T &v)
	{
		if (s >= c)
			reserve(c > 0 ? c * 2 : 1);

		data[s] = v;
		s++;
	}

	void push_back(const T *v, int n)
	{
		for (int i = 0; i < n; i++)
			push_back(v[i]);
	}

	void push_back(const Array<T> &v)
	{
		for (int i = 0; i < v.s; i++)
			push_back(v.data[i]);
	}

	T& push_back()
	{
		if (s >= c)
			reserve(c > 0 ? c * 2 : 1);
		return data[s++];
	}
	
	void pop_back()
	{
		s--;
	}

	void erase_idx(int i)
	{
		data[i] = data[s-1];
		pop_back();
	}

	void erase(const T &v)
	{
		for (int i = 0; i < s; i++)
		{
			if (data[i] == v)
			{
				data[i] = data[s-1];
				pop_back();
				break;
			}
		}
	}

	void erase_ordered(const T &v)
	{
		for (int i = 0; i < s; i++)
		{
			if (data[i] == v)
			{
				for (; i < s-1; i++)
					data[i] = data[i+1];
				pop_back();
				break;
			}
		}
	}

	int count(const T& v)
	{
		int n = 0;
		for (int i = 0; i < s; i++)
		{
			if (data[i] == v)
				n++;
		}
		return n;
	}

	bool contains(const T& v) {return count(v);}

	T* begin(){return data;}
	T* end(){return data + s;}
	
	
	T &first()
	{
		return data[0];
	}

	T &last()
	{
		return data[s-1];
	}

	// math operators
	void operator+=(const Array<T> &a)
	{
		assert(a.s == s);
		for (int i = 0; i < s; i++)
			data[i] += a.data[i];
	}

	void operator-=(const Array<T> &a)
	{
		assert(a.s == s);
		for (int i = 0; i < s; i++)
			data[i] -= a.data[i];
	}

	void operator=(const T &a)
	{
		for (int i = 0; i < s; i++)
			data[i] = a;
	}
};
