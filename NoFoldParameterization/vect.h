/** \file
\brief Tree View of scene
\author Josiah Manson

Modified: May 23, 2006
changed constructors to compile with gcc 4.0

Modified: July 07, 2005
added constructor to vect4<T>

Modified: May 04, 2009
rewrote to be templatized by length
**/

#pragma once

#ifndef __vect_h_
#define __vect_h_

//#pragma warning (disable: 4127) // conditional expression is constant (that is the point. dead code elimination)

#include <assert.h>
#include <cmath>
#include <algorithm>

#undef max
#undef min

// n d vector
template <int N, class T>
struct vect
{
	static const int DIMENSION = N;
	T v[N];

	static vect<N,T> make(const T v0) {vect<N,T> a; a.set(v0); return a;}
	static vect<N,T> make(const T v0,const T v1) {vect<N,T> a; a.set(v0,v1); return a;}
	static vect<N,T> make(const T v0,const T v1,const T v2) {vect<N,T> a; a.set(v0,v1,v2); return a;}
	static vect<N,T> make(const T v0,const T v1,const T v2,const T v3) {vect<N,T> a; a.set(v0,v1,v2,v3); return a;}
	static vect<N,T> make(const T v0,const T v1,const T v2,const T v3,const T v4) {vect<N,T> a; a.set(v0,v1,v2,v3,v4); return a;}
	
	// constructors (true constructors removed because compilers sometimes don't remove the empty constructor in arrays -- this actually was a performance bottleneck once)
	template <class A>
	void setall(const A v0)
	{
		for (int i = 0; i < N; i++)
			v[i] = (T)v0;
	}
	
	template <class A>
	void set(const A v0)
	{
		assert(N == 1);
		v[0] = (T)v0;
	}

	template <class A, class B>
	void set(const A v0, const B v1)
	{
		assert(N == 2);
		v[0] = (T)v0;
		v[1] = (T)v1;
	}
	
	template <class A, class B, class C>
	void set(const A v0, const B v1, const C v2)
	{
		assert(N == 3);
		v[0] = (T)v0;
		v[1] = (T)v1;
		v[2] = (T)v2;
	}
	
	template <class A, class B, class C, class D>
	void set(const A v0, const B v1, const C v2, const D v3)
	{
		assert(N == 4);
		v[0] = (T)v0;
		v[1] = (T)v1;
		v[2] = (T)v2;
		v[3] = (T)v3;
	}

	template <class A, class B, class C, class D, class E>
	void set(const A v0, const B v1, const C v2, const D v3, const E v4)
	{
		assert(N == 5);
		v[0] = (T)v0;
		v[1] = (T)v1;
		v[2] = (T)v2;
		v[3] = (T)v3;
		v[4] = (T)v4;
	}

	template <class S>
	void set(const vect<N, S> &a) 
	{
		for (int i = 0; i < N; i++)
			v[i] = a.v[i];
	}


	void setitem(int i, T val)
	{
		assert(i >= 0 && i < N);
		v[i] = val;
	}

	T getitem(int i)
	{
		assert(i >= 0 && i < N);
		return v[i];
	}

	
	
	template <class A>
	void operator()(const A v0)
	{
		assert(N == 1);
		v[0] = (T)v0;
	}

	template <class A, class B>
	void operator()(const A v0, const B v1)
	{
		assert(N == 2);
		v[0] = (T)v0;
		v[1] = (T)v1;
	}
	
	template <class A, class B, class C>
	void operator()(const A v0, const B v1, const C v2)
	{
		assert(N == 3);
		v[0] = (T)v0;
		v[1] = (T)v1;
		v[2] = (T)v2;
	}
	
	template <class A, class B, class C, class D>
	void operator()(const A v0, const B v1, const C v2, const D v3)
	{
		assert(N == 4);
		v[0] = (T)v0;
		v[1] = (T)v1;
		v[2] = (T)v2;
		v[3] = (T)v3;
	}

	template <class A, class B, class C, class D, class E>
	void operator()(const A v0, const B v1, const C v2, const D v3, const E v4)
	{
		assert(N == 5);
		v[0] = (T)v0;
		v[1] = (T)v1;
		v[2] = (T)v2;
		v[3] = (T)v3;
		v[4] = (T)v4;
	}

	template <int M, class S>
	void operator=(const vect<M, S> &a) 
	{
		int n;
		if (N < M)
			n = N;
		else
			n = M;

		for (int i = 0; i < n; i++)
			v[i] = a.v[i];
	}

	template <class S>
	void operator=(const S a) 
	{
		for (int i = 0; i < N; i++)
			v[i] = (T)a;
	}

	// accessors
	T &operator[](const int i)
	{
		assert(i >= 0 && i < N);
		return v[i];
	}

	// scalar functions
	template <class S> 
	void add(const S a)
	{
		for (int i = 0; i < N; i++)
			v[i] += (T)a;
	}

	void subtract(const T a)
	{
		for (int i = 0; i < N; i++)
			v[i] -= a;
	}
	
	vect<N, T> multiply(const vect<N, T> &a)
	{
		vect<N, T> r;
		for (int i = 0; i < N; i++)
			r[i] = v[i] * a.v[i];
		return r;
	}

	void multiply(const T a)
	{
		for (int i = 0; i < N; i++)
			v[i] *= a;
	}
	
	void divide(const T a)
	{
		for (int i = 0; i < N; i++)
			v[i] /= a;
	}

	// vector functions
	template <class S>
	void add(const vect<N, S> &a)
	{
		for (int i = 0; i < N; i++)
			v[i] += a.v[i];
	}

	template <class S>
	void subtract(const vect<N, S> &a)
	{
		for (int i = 0; i < N; i++)
			v[i] -= a.v[i];
	}

	template <class S>
	T dot(const vect<N, S> &a) const
	{
		T r = 0;
		for (int i = 0; i < N; i++)
			r += v[i]*a.v[i];
		return r;
	}
	
	template <class S>
	vect<3, T> cross(const vect<N, S> &a) const
	{
		vect<3, T> r;

		assert(N == 3);
		
		r.v[0] = v[1] * a.v[2] - v[2] * a.v[1];
		r.v[1] = v[2] * a.v[0] - v[0] * a.v[2];
		r.v[2] = v[0] * a.v[1] - v[1] * a.v[0];

		return r;
	}



	// unary functions
	void negate()
	{
		for (int i = 0; i < N; i++)
			v[i] = -v[i];
	}

	T length() const
	{
		return (T)sqrt((double)length2());
	}
	T length_safe() const
	{
		T l2 = length2();
		if (l2 < 1e-6)
			return 0;
		else
			return (T)sqrt(l2);
	}
	T length2() const
	{
		return dot(*this);
	}
	void normalize()
	{
		const T f = (T)(1.0 / length());
		for (int i = 0; i < N; i++)
			v[i] *= f;
	}

	// vector operators
	template <class S> 
	vect<3, T> operator%(const vect<N, S> &a) const 
	{
		return cross(a);
	}

	template <class S> 
	vect<N, T> operator+(const vect<N, S> &a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] + a.v[i];
		return r;
	}
	
	template <class S> 
	vect<N, T> operator-(const vect<N, S> &a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] - a.v[i];
		return r;
	}
	
	template <class S> 
	vect<N, T> operator/(const vect<N, S> &a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] / a.v[i];
		return r;
	}
	
	template <int M, class S> 
	T operator*(const vect<M, S> &a) const 
	{
		return dot(a);
	}
	
	// scalar operators
	template <class S> 
	friend vect<N, T> operator*(const S a, vect<N, T>);

	template <class S> 
	vect<N, T> operator*(const S a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = (T)(v[i] * a);
		return r;
	}
	template <class S> 
	vect<N, T> operator/(const S a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] / a;
		return r;
	}
	vect<N, T> operator+(const T a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] + a;
		return r;
	}
	vect<N, T> operator-(const T a) const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] - a;
		return r;
	}

	// unary operators
	vect<N, T> operator-() const 
	{
		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = -v[i];
		return r;
	}
	
	T operator!() const 
	{
		return length();
	}
	
	vect<N, T> operator~() const 
	{
		double f = 1.0 / length();

		vect<N, T> r; 
		for (int i = 0; i < N; i++)
			r.v[i] = v[i] * f;
		return r;
	}

	// update operators
	template <int M, class S> 
	void operator+=(const vect<M, S> &a)
	{
		for (int i = 0; i < N; i++)
			v[i] += (T)a.v[i];
	}
	template <int M, class S> 
	void operator-=(const vect<M, S> &a)
	{
		for (int i = 0; i < N; i++)
			v[i] -= (T)a.v[i];
	}
	template <int M, class S> 
	void operator*=(const vect<M, S> &a)
	{
		for (int i = 0; i < N; i++)
			v[i] *= (T)a.v[i];
	}
	template <int M, class S> 
	void operator/=(const vect<M, S> &a)
	{
		for (int i = 0; i < N; i++)
			v[i] /= (T)a.v[i];
	}

	template <class S> 
	void operator+=(const S a)
	{
		for (int i = 0; i < N; i++)
			v[i] += (T)a;
	}
	template <class S> 
	void operator-=(const S a)
	{
		for (int i = 0; i < N; i++)
			v[i] -= (T)a;
	}
	template <class S> 
	void operator*=(const S a)
	{
		for (int i = 0; i < N; i++)
			v[i] *= (T)a;
	}
	template <class S> 
	void operator/=(const S a)
	{
		for (int i = 0; i < N; i++)
			v[i] /= (T)a;
	}

	// comparison functions
	bool operator==(const vect<N, T> &a) const 
	{
		for (int i = 0; i < N; i++)
			if (v[i] != a.v[i])
				return false;
		return true;
	}
	bool operator!=(const vect<N, T> &a) const 
	{
		return !(*this == a);
	}
	bool operator<(const vect<N, T> &a) const 
	{
		return lexicographical_compare(v, v+N, a.v, a.v+N);
	}

	// utility
	T max()
	{
		T m = v[0];
		for (int i = 1; i < N; i++)
			if (m < v[i])
				m = v[i];
		return m;
	}
	T min()
	{
		T m = v[0];
		for (int i = 1; i < N; i++)
			if (m > v[i])
				m = v[i];
		return m;
	}
};

template <int N, class S, class T> 
vect<N, T> operator*(const S a, vect<N, T>);

#ifndef WIN32
#define __int8 char
#define __int16 short
#define __int32 int
#define __int64 long long
#endif

typedef vect<1, float> vect1f;
typedef vect<1, double> vect1d;
typedef vect<1, __int8> vect1b;
typedef vect<1, __int16> vect1s;
typedef vect<1, __int32> vect1i;
typedef vect<1, __int64> vect1l;
typedef vect<1, unsigned __int8> vect1ub;
typedef vect<1, unsigned __int16> vect1us;
typedef vect<1, unsigned __int32> vect1ui;
typedef vect<1, unsigned __int64> vect1ul;

typedef vect<2, float> vect2f;
typedef vect<2, double> vect2d;
typedef vect<2, __int8> vect2b;
typedef vect<2, __int16> vect2s;
typedef vect<2, __int32> vect2i;
typedef vect<2, __int64> vect2l;
typedef vect<2, unsigned __int8> vect2ub;
typedef vect<2, unsigned __int16> vect2us;
typedef vect<2, unsigned __int32> vect2ui;
typedef vect<2, unsigned __int64> vect2ul;

typedef vect<3, float> vect3f;
typedef vect<3, double> vect3d;
typedef vect<3, __int8> vect3b;
typedef vect<3, __int16> vect3s;
typedef vect<3, __int32> vect3i;
typedef vect<3, __int64> vect3l;
typedef vect<3, unsigned __int8> vect3ub;
typedef vect<3, unsigned __int16> vect3us;
typedef vect<3, unsigned __int32> vect3ui;
typedef vect<3, unsigned __int64> vect3ul;

typedef vect<4, float> vect4f;
typedef vect<4, double> vect4d;
typedef vect<4, __int8> vect4b;
typedef vect<4, __int16> vect4s;
typedef vect<4, __int32> vect4i;
typedef vect<4, __int64> vect4l;
typedef vect<4, unsigned __int8> vect4ub;
typedef vect<4, unsigned __int16> vect4us;
typedef vect<4, unsigned __int32> vect4ui;
typedef vect<4, unsigned __int64> vect4ul;

typedef vect<5, float> vect5f;
typedef vect<5, double> vect5d;
typedef vect<5, __int8> vect5b;
typedef vect<5, __int16> vect5s;
typedef vect<5, __int32> vect5i;
typedef vect<5, __int64> vect5l;
typedef vect<5, unsigned __int8> vect5ub;
typedef vect<5, unsigned __int16> vect5us;
typedef vect<5, unsigned __int32> vect5ui;
typedef vect<5, unsigned __int64> vect5ul;

#endif // __vect_h_
