#include "iso_pre.h"
#include "vect.h"

iso_pre::iso_pre()
{
}

iso_pre::iso_pre(double L, double x, double y)
{
	vect2d d21;
	d21[0] = L;
	d21[1] = 0;

	vect2d d31;
	d31[0] = x;
	d31[1] = y;

	vect2d d32;
	d32[0] = x-L;
	d32[1] = y;

	//ld21 = d21.length2();
	ld21 = L*L;

	//ld31 = d31.length2();
	ld31 = x*x + y*y;

	//ld32 = d32.length2();
	ld32 = (x-L)*(x-L) + y*y;

	Pt = L*y / 2.0f;
	Pt2 = Pt*Pt;

//	p3dotp2 = d31.dot(d21);
	p3dotp2 = L*x;

	/*cot[0] = -1*x / y;
	cot[1] = (x-L)/y;
	cot[2] = -1*(((x-L)*x) + y*y) / (L*y);
*/
	cot[0] = x / y;
	cot[1] = (L-x)/y;
	cot[2] = (-1*(L-x)*x + y*y) / (L*y);

	area_uv = 0;


	// scott precomp... these just encode the symmetric hessian coefficients from the quadratic part of the energy * (Ly)^2, should just be related to cotan weights
	c00 = ((L - x)*(L - x) + y * y) * (L * y);
	c01 = (x * (L - x) - y * y) * (L * y);
	c02 = (x - L) * L * L * y;
	c11 = (x * x + y * y) * (L * y);
	c12 = -x * L * L * y;
	c22 = L * L * L * y;
	ly = L * y;

	this->L = L;
	this->x = x;
	this->y = y;
}