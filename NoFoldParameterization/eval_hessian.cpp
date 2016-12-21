#include <math.h>

void jason_hessian(double **hessian, double L, double xx, double y, double u0, double v0, double u1, double v1, double u2, double v2, int i, int ip1, int j)
{
	double a, b, c, d, A, x23y23, x13y23, x12y23,
		x23y13, x13y13, x12y13,
		x23y12, x13y12, x12y12,
		x12, x13, x23,
		y12, y13, y23;

	a = u1* y - u0* y  - L* v0 + xx* v0 - xx* v1 + L* v2;
	b = L* u0 - xx* u0 + xx* u1 - L* u2 - y* v0 + y* v1;
	c = u1* y - u0* y + L* v0 - xx* v0 + xx* v1 - L* v2;
	d = xx*(u0 - u1) + L*(u2-u0 ) + y*(v1-v0);
	A = (u2*(v1-v0) + u1*(v0 - v2) + u0*(v2-v1));
	x23y23 = (u1 - u2)* (v1 - v2);
	x13y23 = (u0 - u2)* (v1 - v2);
	x12y23 = (u0 - u1)* (v1 - v2);

	x23y13 = (u1 - u2)* (v0 - v2);
	x13y13 = (u0 - u2)* (v0 - v2);
	x12y13 = (u0 - u1)* (v0 - v2);

	x23y12 = (u1 - u2)* (v0 - v1);
	x13y12 = (u0 - u2)* (v0 - v1);
	x12y12 = (u0 - u1)* (v0 - v1);

	x12 = (u0 - u1)*(u0 - u1);
	x13 = (u0 - u2)*(u0 - u2);
	x23 = (u1 - u2)*(u1 - u2);
	y12 = (v0 - v1)*(v0 - v1);
	y13 = (v0 - v2)*(v0 - v2);
	y23 = (v1 - v2)*(v1 - v2);

	//Du0 i... 
	hessian[2 * i][2 * i] += (2 * A*(A*(pow(L - xx, 2) + pow(y, 2))*
		(pow(A, 2) + pow(L, 2)*pow(y, 2)) -
		pow(A, 2)*(b*(L - xx) + d*(-L + xx) - (a + c)*y)*(v1 - v2) +
		((b*(L - xx) + d*(-L + xx) - (a + c)*y)*
		(pow(A, 2) + pow(L, 2)*pow(y, 2))*(v1 - v2)) / 2. +
		pow(A, 2)*(b*(L - xx) + d*(-L + xx) - (a + c)*y)*(-v1 + v2))\
		- 3 * (-v1 + v2)*(A*(b*(L - xx) + d*(-L + xx) - (a + c)*y)*
		(pow(A, 2) + pow(L, 2)*pow(y, 2)) -
		pow(A, 2)*(pow(a, 2) + pow(b, 2) + pow(c, 2) +
		pow(d, 2))*(v1 - v2) -
		((pow(a, 2) + pow(b, 2) + pow(c, 2) + pow(d, 2))*
		(b*L*(-L + xx)*y + d*L*(-L + xx)*y + a*L*pow(y, 2) -
		c*L*pow(y, 2) - 2 * pow(A, 2)*v1 + 2 * pow(A, 2)*v2)) / 2.
		)) / (pow(A, 4)*L*y);
	hessian[2 * i][2 * i + 1] += -((L*y*(a*(d*(2 * pow(u1, 2) - 4 * u1*u2 + 2 * pow(u2, 2) +
		pow(v1 - v2, 2)) +
		(2 * c*(u1 - u2) + b*(v1 - v2))*(v1 - v2)) -
		c*(4 * A*(u1 - u2)*y +
		d*(2 * pow(u1, 2) - 4 * u1*u2 + 2 * pow(u2, 2) +
		pow(v1 - v2, 2)) - 3 * b*pow(v1 - v2, 2)) +
		pow(a, 2)*(u1 - u2)*(v1 - v2) +
		3 * pow(c, 2)*(u1 - u2)*(v1 - v2) +
		2 * (pow(d, 2)*(u1 - u2) - A*b*y + A*d*y)*(v1 - v2))) /
		pow(A, 4));
	hessian[2 * i][2 * ip1] += (2 * (pow(A, 4)*L*xx - pow(A, 4)*(pow(xx, 2) + pow(y, 2)) +
		pow(L, 4)*pow(y, 2)*
		(-2 * A*(u0 - u2) - 3 *
		(pow(u0, 2) - 2 * u0*u2 + pow(u2, 2) + pow(v0 - v2, 2))*
		(v1 - v2))*(v0 - v2) -
		pow(L, 2)*pow(y, 2)*(pow(xx, 2) + pow(y, 2))*(v0 - v1)*
		((v0 - v1)*(pow(u2, 2) + 3 * (v0 - v2)*(v1 - v2)) +
		pow(u0, 2)*(v1 - v2) + 2 * u1*u2*(v1 - v2) +
		pow(u1, 2)*(-v0 + v2) +
		2 * u0*(u1*(v0 - v1) + u2*(-v0 + v2))) +
		pow(L, 3)*xx*pow(y, 2)*
		((v0 - v1)*(pow(u2, 2)*(3 * v0 + v1 - 4 * v2) +
		6 * pow(v0 - v2, 2)*(v1 - v2)) -
		pow(u1, 2)*pow(v0 - v2, 2) +
		pow(u0, 2)*(2 * v0 - v1 - v2)*(v1 - v2) -
		2 * u1*u2*(v0 - v2)*(v0 - 2 * v1 + v2) +
		2 * u0*(2 * u1*(v0 - v1)*(v0 - v2) +
		u2*(-2 * pow(v0, 2) + pow(v1, 2) + 4 * v0*v2 -
		v2*(2 * v1 + v2)))))) / (pow(A, 4)*L*y);
	hessian[2 * i][2 * ip1 + 1] += (L*y*(-2 * c*d*u0*u1 + 2 * c*d*u0*u2 + 2 * c*d*u1*u2 -
		2 * c*d*pow(u2, 2) + 3 * pow(c, 2)*u0*v1 +
		2 * pow(d, 2)*u0*v1 - 3 * pow(c, 2)*u2*v1 -
		2 * pow(d, 2)*u2*v1 + 3 * b*c*v0*v1 - c*d*v0*v1 +
		pow(a, 2)*(A + (u0 - u2)*(v1 - v2)) - 3 * pow(c, 2)*u0*v2 -
		2 * pow(d, 2)*u0*v2 + 3 * pow(c, 2)*u2*v2 +
		2 * pow(d, 2)*u2*v2 - 3 * b*c*v0*v2 + c*d*v0*v2 - 3 * b*c*v1*v2 +
		c*d*v1*v2 + 3 * b*c*pow(v2, 2) - c*d*pow(v2, 2) +
		A*(pow(b, 2) + 3 * pow(c, 2) + 4 * c*(-u1 + u2)*y +
		2 * b*y*(-v1 + v2) + d*(d + 2 * y*v1 - 2 * y*v2)) +
		a*(2 * A*c + (v1 - v2)*(2 * c*u0 - 2 * c*u2 + b*v0 - b*v2) +
		d*(2 * u0*u1 - 2 * u0*u2 - 2 * u1*u2 + 2 * pow(u2, 2) + v0*v1 -
		v0*v2 - v1*v2 + pow(v2, 2))))) / pow(A, 4);
	hessian[2 * i][2 * j] += -(pow(c, 4)*L + 4 * pow(c, 2)*pow(d, 2)*L + 3 * pow(d, 4)*L +
		pow(a, 4)*(L - xx) + 3 * pow(b, 4)*(L - xx) - pow(c, 4)*xx -
		4 * pow(c, 2)*pow(d, 2)*xx - 3 * pow(d, 4)*xx +
		8 * pow(A, 3)*c*u1 - 8 * pow(A, 3)*c*u2 + 2 * pow(c, 3)*d*y +
		2 * c*pow(d, 3)*y - 2 * pow(a, 3)*(b + 3 * d)*y +
		6 * pow(b, 3)*(2 * d*(L - xx) + c*y) +
		2 * pow(b, 2)*(4 * pow(c, 2)*(L - xx) + 9 * pow(d, 2)*(L - xx) +
		5 * c*d*y) - 2 * a*(4 * pow(A, 3)*(u1 - u2) +
		(pow(b, 3) + 3 * pow(b, 2)*d +
		5 * b*(pow(c, 2) + pow(d, 2)) +
		3 * d*(pow(c, 2) + pow(d, 2)))*y) +
		2 * pow(a, 2)*(2 * pow(b, 2)*(L - xx) + 4 * pow(d, 2)*(L - xx) +
		pow(c, 2)*(-L + xx) + 5 * c*d*y + 3 * b*(2 * d*(L - xx) + c*y)) +
		8 * pow(A, 3)*d*v1 - 8 * pow(A, 3)*d*v2 +
		2 * b*(6 * pow(c, 2)*d*(L - xx) + 6 * pow(d, 3)*(L - xx) +
		3 * pow(c, 3)*y + 3 * c*pow(d, 2)*y + 4 * pow(A, 3)*(-v1 + v2)
		)) / (8.*pow(A, 4)*y);
	hessian[2 * i][2 * j + 1] += (pow(a, 4) + pow(b, 4) - 3 * pow(c, 4) - 4 * pow(c, 2)*pow(d, 2) -
		pow(d, 4) - 4 * pow(c, 2)*d*L*u1 + 4 * pow(c, 2)*d*L*u2 +
		6 * pow(c, 3)*L*v1 + 2 * c*pow(d, 2)*L*v1 +
		2 * pow(a, 3)*(c + L*(v1 - v2)) +
		2 * pow(b, 2)*c*(c + 3 * L*(v1 - v2)) +
		2 * pow(a, 2)*(pow(b, 2) + pow(c, 2) + 2 * d*L*(u1 - u2) +
		3 * c*L*(v1 - v2)) + 2 * a*
		(-pow(c, 3) - c*pow(d, 2) + pow(b, 2)*(c + L*(v1 - v2)) +
		5 * pow(c, 2)*L*(v1 - v2) + 2 * b*d*L*(v1 - v2) +
		3 * pow(d, 2)*L*(v1 - v2)) + 4 * b*c*d*L*(v1 - v2) -
		6 * pow(c, 3)*L*v2 - 2 * c*pow(d, 2)*L*v2) / (4.*pow(A, 4));

	/*hessian[2 * i][2 * i + 1] = hessian[2 * i + 1][2 * i];
	hessian[2 * i][2 * ip1] = hessian[2 * ip1][2 * i];
	hessian[2 * i][2 * ip1 + 1] = hessian[2 * ip1 + 1][2 * i];
	hessian[2 * i][2 * j] = hessian[2 * j][2 * i];
	hessian[2 * i][2 * j + 1] = hessian[2 * j + 1][2 * i];*/

	//Dv0
	hessian[2 * i + 1][2 * i + 1] += (2 * A*(2 * pow(A, 2)*(u1 - u2)*(c*(L - xx) + a*(-L + xx) - (b + d)*y) +
		((u1 - u2)*(a*(L - xx) + c*(-L + xx) + (b + d)*y)*
		(pow(A, 2) + pow(L, 2)*pow(y, 2))) / 2. +
		A*(pow(L - xx, 2) + pow(y, 2))*
		(pow(A, 2) + pow(L, 2)*pow(y, 2))) -
		3 * (u1 - u2)*(-((pow(a, 2) + pow(b, 2) + pow(c, 2) +
		pow(d, 2))*(u1 - u2)*
		(pow(A, 2) + pow(L, 2)*pow(y, 2))) +
		A*(c*(L - xx) + a*(-L + xx) - (b + d)*y)*
		(pow(A, 2) + pow(L, 2)*pow(y, 2)) +
		pow(A, 2)*(pow(a, 2)*(u1 - u2) + pow(c, 2)*(u1 - u2) +
		2 * (pow(d, 2)*(u1 - u2) - A*b*y + A*d*y) +
		a*(b - d)*(v1 - v2) + c*(b - d)*(v1 - v2)))) /
		(pow(A, 4)*L*y);
	hessian[2 * i + 1][2 * ip1] += (L*y*(-2 * c*d*u0*u1 + 2 * c*d*u0*u2 + 2 * c*d*u1*u2 -
		2 * c*d*pow(u2, 2) + 3 * pow(c, 2)*u0*v1 +
		2 * pow(d, 2)*u0*v1 - 3 * pow(c, 2)*u2*v1 -
		2 * pow(d, 2)*u2*v1 + 3 * b*c*v0*v1 - c*d*v0*v1 +
		2 * A*(pow(c, 2) + 2 * c*(-u1 + u2)*y - (b - d)*y*(v1 - v2)) +
		pow(a, 2)*(u0 - u2)*(v1 - v2) - 3 * pow(c, 2)*u0*v2 -
		2 * pow(d, 2)*u0*v2 + 3 * pow(c, 2)*u2*v2 +
		2 * pow(d, 2)*u2*v2 - 3 * b*c*v0*v2 + c*d*v0*v2 - 3 * b*c*v1*v2 +
		c*d*v1*v2 + 3 * b*c*pow(v2, 2) - c*d*pow(v2, 2) +
		a*(2 * A*c + (v1 - v2)*(2 * c*u0 - 2 * c*u2 + b*v0 - b*v2) +
		d*(2 * u0*u1 - 2 * u0*u2 - 2 * u1*u2 + 2 * pow(u2, 2) + v0*v1 -
		v0*v2 - v1*v2 + pow(v2, 2))))) / pow(A, 4);
	hessian[2 * i + 1][2 * ip1 + 1] += (2 * (pow(A, 4)*L*xx - pow(A, 4)*(pow(xx, 2) + pow(y, 2)) -
		pow(L, 2)*(u0 - u1)*pow(y, 2)*(pow(xx, 2) + pow(y, 2))*
		(3 * pow(u0, 2)*(u1 - u2) + 3 * pow(u1, 2)*u2 -
		u2*(v0 - v1)*(v0 + v1 - 2 * v2) +
		u0*(-3 * pow(u1, 2) + 3 * pow(u2, 2) + 2 * v0*v1 -
		pow(v1, 2) - 2 * v0*v2 + pow(v2, 2)) +
		u1*(-3 * pow(u2, 2) + (v0 - v2)*(v0 - 2 * v1 + v2))) +
		pow(L, 4)*(u0 - u2)*pow(y, 2)*
		(3 * pow(u0, 2)*(-u1 + u2) +
		u2*(3 * pow(u2, 2) + (v0 + 2 * v1 - 3 * v2)*(v0 - v2)) -
		u1*(3 * pow(u2, 2) + pow(v0 - v2, 2)) +
		2 * u0*(3 * (u1 - u2)*u2 + (v0 - v2)*(-v1 + v2))) +
		pow(L, 3)*xx*pow(y, 2)*
		(6 * pow(u0, 3)*(u1 - u2) -
		pow(u1, 2)*(6 * pow(u2, 2) + pow(v0 - v2, 2)) -
		pow(u0, 2)*(6 * (u1 - u2)*(u1 + 2 * u2) -
		(4 * v0 - v1 - 3 * v2)*(v1 - v2)) +
		2 * u1*u2*(3 * pow(u2, 2) + 2 * (v0 - v2)*(v1 - v2)) +
		pow(u2, 2)*(v0 - v1)*(v0 + v1 - 2 * v2) +
		2 * u0*(6 * pow(u1, 2)*u2 -
		u2*(3 * pow(u2, 2) + pow(v0, 2) + 2 * v0*v1 -
		pow(v1, 2) - 4 * v0*v2 + 2 * pow(v2, 2)) +
		u1*(-3 * pow(u2, 2) + (v0 - v2)*(v0 - 2 * v1 + v2)))))) /
		(pow(A, 4)*L*y);
	hessian[2 * i + 1][2 * j] += (pow(a, 3)*(c + L*(v1 - v2)) +
		pow(a, 2)*(pow(c, 2) + 2 * d*L*(u1 - u2) + 3 * c*L*(v1 - v2)) +
		a*(-pow(c, 3) - c*pow(d, 2) + pow(b, 2)*(c + L*(v1 - v2)) +
		5 * pow(c, 2)*L*(v1 - v2) + 2 * b*d*L*(v1 - v2) +
		3 * pow(d, 2)*L*(v1 - v2)) -
		c*(pow(c, 3) + c*d*(d + 2 * L*u1 - 2 * L*u2) -
		pow(b, 2)*(c + 3 * L*(v1 - v2)) +
		3 * pow(c, 2)*L*(-v1 + v2) + 2 * b*d*L*(-v1 + v2) +
		pow(d, 2)*L*(-v1 + v2))) / (2.*pow(A, 4));
	hessian[2 * i + 1][2 * j + 1] += (-3 * pow(c, 4)*L - 4 * pow(c, 2)*pow(d, 2)*L - pow(d, 4)*L -
		3 * pow(a, 4)*(L - xx) + 3 * pow(c, 4)*xx +
		4 * pow(c, 2)*pow(d, 2)*xx + pow(d, 4)*xx +
		pow(b, 4)*(-L + xx) - 8 * pow(A, 3)*c*u1 + 8 * pow(A, 3)*c*u2 -
		6 * pow(b, 3)*c*y + 2 * pow(c, 3)*d*y + 2 * c*pow(d, 3)*y -
		2 * pow(a, 3)*(6 * c*(L - xx) + (b - 3 * d)*y) -
		2 * pow(a, 2)*(2 * pow(b, 2)*(L - xx) + 9 * pow(c, 2)*(L - xx) +
		4 * pow(d, 2)*(L - xx) + 3 * b*c*y - 5 * c*d*y) +
		pow(b, 2)*(-8 * pow(c, 2)*(L - xx) + 2 * pow(d, 2)*(L - xx) +
		10 * c*d*y) - 2 * a*(6 * pow(c, 3)*(L - xx) +
		6 * c*pow(d, 2)*(L - xx) - 4 * pow(A, 3)*u1 +
		4 * pow(A, 3)*u2 + pow(b, 3)*y - 3 * pow(c, 2)*d*y -
		3 * pow(d, 3)*y + 5 * b*(pow(c, 2) + pow(d, 2))*y +
		pow(b, 2)*(6 * c*(L - xx) - 3 * d*y)) - 8 * pow(A, 3)*d*v1 +
		8 * pow(A, 3)*d*v2 - 2 * b*
		(3 * pow(c, 3)*y + 3 * c*pow(d, 2)*y + 4 * pow(A, 3)*(-v1 + v2)))
		/ (8.*pow(A, 4)*y);

	/*hessian[2 * i + 1][2 * ip1] = hessian[2 * ip1][2 * i + 1];
	hessian[2 * i + 1][2 * ip1 + 1] = hessian[2 * ip1 + 1][2 * i + 1];
	hessian[2 * i + 1][2 * j] = hessian[2 * j][2 * i + 1];
	hessian[2 * i + 1][2 * j + 1] = hessian[2 * j + 1][2 * i + 1];
*/
	//Du1 i1p...
	hessian[2 * ip1][2 * ip1] += (2 * pow(A, 4)*(pow(xx, 2) + pow(y, 2)) +
		6 * pow(L, 4)*pow(y, 2)*(pow(u0 - u2, 2) + pow(v0 - v2, 2))*
		pow(v0 - v2, 2) + 2 * pow(L, 2)*pow(y, 2)*
		(pow(xx, 2) + pow(y, 2))*(v0 - v1)*
		((v0 - v1)*(pow(u2, 2) + 3 * pow(v0 - v2, 2)) +
		pow(u0, 2)*(3 * v0 - v1 - 2 * v2) + 2 * u1*u2*(v0 - v2) +
		2 * u0*(u1*(-v0 + v2) + u2*(-2 * v0 + v1 + v2))) -
		4 * pow(L, 3)*xx*pow(y, 2)*(v0 - v2)*
		((v0 - v1)*(2 * pow(u2, 2) + 3 * pow(v0 - v2, 2)) +
		u1*u2*(v0 - v2) + pow(u0, 2)*(3 * v0 - 2 * v1 - v2) +
		u0*(u1*(-v0 + v2) + u2*(-5 * v0 + 4 * v1 + v2)))) /
		(pow(A, 4)*L*y);
	hessian[2 * ip1][2 * ip1 + 1] += -((L*y*(-2 * c*d*pow(u0, 2) + 4 * c*d*u0*u2 - 2 * c*d*pow(u2, 2) +
		3 * pow(c, 2)*u0*v0 + 2 * pow(d, 2)*u0*v0 -
		3 * pow(c, 2)*u2*v0 - 2 * pow(d, 2)*u2*v0 +
		3 * b*c*pow(v0, 2) - c*d*pow(v0, 2) +
		a*(d*(2 * pow(u0, 2) - 4 * u0*u2 + 2 * pow(u2, 2) +
		pow(v0 - v2, 2)) +
		(2 * c*(u0 - u2) + b*(v0 - v2))*(v0 - v2)) +
		pow(a, 2)*(u0 - u2)*(v0 - v2) - 3 * pow(c, 2)*u0*v2 -
		2 * pow(d, 2)*u0*v2 + 3 * pow(c, 2)*u2*v2 +
		2 * pow(d, 2)*u2*v2 - 6 * b*c*v0*v2 + 2 * c*d*v0*v2 +
		3 * b*c*pow(v2, 2) - c*d*pow(v2, 2) +
		A*(pow(b, 2) + 2 * a*c + 2 * pow(c, 2) - pow(d, 2) -
		4 * c*u1*y + 4 * c*u2*y + 2 * d*y*v1 - 2 * d*y*v2 +
		2 * b*y*(-v1 + v2)))) / pow(A, 4));
	hessian[2 * ip1][2 * j] += -(pow(a, 4)*xx + 3 * pow(b, 4)*xx + pow(c, 4)*xx +
		4 * pow(c, 2)*pow(d, 2)*xx + 3 * pow(d, 4)*xx -
		8 * pow(A, 3)*c*u0 + 8 * pow(A, 3)*c*u2 - 2 * pow(c, 3)*d*y -
		2 * c*pow(d, 3)*y + 2 * pow(a, 3)*(b + 3 * d)*y +
		6 * pow(b, 3)*(2 * d*xx - c*y) +
		2 * pow(b, 2)*(4 * pow(c, 2)*xx + 9 * pow(d, 2)*xx - 5 * c*d*y) +
		2 * pow(a, 2)*(2 * pow(b, 2)*xx - pow(c, 2)*xx + 6 * b*d*xx +
		4 * pow(d, 2)*xx - 3 * b*c*y - 5 * c*d*y) +
		2 * a*(4 * pow(A, 3)*(u0 - u2) +
		(pow(b, 3) + 3 * pow(b, 2)*d +
		5 * b*(pow(c, 2) + pow(d, 2)) +
		3 * d*(pow(c, 2) + pow(d, 2)))*y) - 8 * pow(A, 3)*d*v0 +
		8 * pow(A, 3)*d*v2 - 2 * b*
		(-6 * pow(c, 2)*d*xx - 6 * pow(d, 3)*xx + 3 * pow(c, 3)*y +
		3 * c*pow(d, 2)*y + 4 * pow(A, 3)*(-v0 + v2))) /
		(8.*pow(A, 4)*y);
	hessian[2 * ip1][2 * j + 1] += (-pow(b, 4) + pow(c, 2)*pow(d, 2) + pow(d, 4) +
		4 * pow(c, 2)*d*xx*u0 - 4 * pow(c, 2)*d*xx*u1 -
		6 * pow(c, 3)*xx*v0 - 2 * c*pow(d, 2)*xx*v0 -
		pow(a, 2)*(pow(b, 2) - 5 * pow(d, 2) + 4 * d*xx*(u0 - u1) +
		6 * c*xx*(v0 - v1)) - 2 * a*
		(pow(b, 2) + 5 * pow(c, 2) + 2 * b*d + 3 * pow(d, 2))*xx*(v0 - v1)
		+ 6 * pow(c, 3)*xx*v1 + 2 * c*pow(d, 2)*xx*v1 +
		2 * pow(a, 3)*xx*(-v0 + v1) + 4 * b*c*d*xx*(-v0 + v1) +
		pow(b, 2)*c*(-5 * c + 6 * xx*(-v0 + v1))) / (4.*pow(A, 4));

	/*hessian[2 * ip1][2 * ip1 + 1] = hessian[2 * ip1 + 1][2 * ip1];
	hessian[2 * ip1][2 * j] = hessian[2 * j][2 * ip1];
	hessian[2 * ip1][2 * j + 1] = hessian[2 * j + 1][2 * ip1];*/

	//Dv1
	hessian[2 * ip1 + 1][2 * ip1 + 1] += (2 * (pow(A, 4)*(pow(xx, 2) + pow(y, 2)) +
		pow(L, 2)*(u0 - u1)*pow(y, 2)*(pow(xx, 2) + pow(y, 2))*
		(3 * pow(u0, 3) - 3 * pow(u0, 2)*(u1 + 2 * u2) -
		u1*(3 * pow(u2, 2) + pow(v0 - v2, 2)) +
		u0*(3 * u2*(2 * u1 + u2) + (v0 - v2)*(3 * v0 - 2 * v1 - v2)) -
		2 * u2*(v0 - v1)*(v0 - v2)) -
		2 * pow(L, 3)*xx*(u0 - u2)*pow(y, 2)*
		(3 * pow(u0, 3) - 3 * pow(u0, 2)*(u1 + 2 * u2) +
		u0*(3 * u2*(2 * u1 + u2) + (3 * v0 - v1 - 2 * v2)*(v0 - v2)) -
		u1*(3 * pow(u2, 2) + 2 * pow(v0 - v2, 2)) -
		u2*(v0 - v1)*(v0 - v2)) +
		3 * pow(L, 4)*pow(u0 - u2, 2)*pow(y, 2)*
		(pow(u0 - u2, 2) + pow(v0 - v2, 2)))) / (pow(A, 4)*L*y);
	hessian[2 * ip1 + 1][2 * j] += (pow(a, 4) - c*(pow(c, 3) + c*d*(d - 4 * xx*u0 + 4 * xx*u1) +
		pow(b, 2)*(5 * c + 6 * xx*(v0 - v1)) +
		6 * pow(c, 2)*xx*(v0 - v1) + 4 * b*d*xx*(v0 - v1) +
		2 * pow(d, 2)*xx*(v0 - v1)) -
		2 * a*(pow(b, 2) + 5 * pow(c, 2) + 2 * b*d + 3 * pow(d, 2))*xx*
		(v0 - v1) + 2 * pow(a, 3)*xx*(-v0 + v1) +
		pow(a, 2)*(pow(b, 2) + 5 * pow(d, 2) + 4 * d*xx*(-u0 + u1) +
		6 * c*xx*(-v0 + v1))) / (4.*pow(A, 4));
	hessian[2 * ip1 + 1][2 * j + 1] += -(3 * pow(a, 4)*xx + pow(b, 4)*xx + 3 * pow(c, 4)*xx +
		4 * pow(c, 2)*pow(d, 2)*xx + pow(d, 4)*xx -
		8 * pow(A, 3)*c*u0 + 8 * pow(A, 3)*c*u2 - 6 * pow(b, 3)*c*y +
		2 * pow(c, 3)*d*y + 2 * c*pow(d, 3)*y +
		2 * pow(a, 3)*(6 * c*xx - b*y + 3 * d*y) +
		2 * pow(b, 2)*(4 * pow(c, 2)*xx - pow(d, 2)*xx + 5 * c*d*y) +
		2 * pow(a, 2)*(2 * pow(b, 2)*xx + 9 * pow(c, 2)*xx +
		4 * pow(d, 2)*xx - 3 * b*c*y + 5 * c*d*y) -
		2 * a*(-6 * pow(c, 3)*xx - 6 * c*pow(d, 2)*xx - 4 * pow(A, 3)*u0 +
		4 * pow(A, 3)*u2 + pow(b, 3)*y - 3 * pow(c, 2)*d*y -
		3 * pow(d, 3)*y + 5 * b*(pow(c, 2) + pow(d, 2))*y -
		3 * pow(b, 2)*(2 * c*xx + d*y)) - 8 * pow(A, 3)*d*v0 +
		8 * pow(A, 3)*d*v2 - 2 * b*
		(3 * pow(c, 3)*y + 3 * c*pow(d, 2)*y + 4 * pow(A, 3)*(-v0 + v2))
		) / (8.*pow(A, 4)*y);

	/*hessian[2 * ip1 + 1][2 * j] = hessian[2 * j][2 * ip1 + 1];
	hessian[2 * ip1 + 1][2 * j + 1] = hessian[2 * j + 1][2 * ip1 + 1];*/

	//Du2  j...
	hessian[2 * j][2 * j] += (pow(a, 4)*L + 3 * pow(b, 4)*L + pow(c, 4)*L + 12 * pow(b, 3)*d*L +
		4 * pow(c, 2)*pow(d, 2)*L + 3 * pow(d, 4)*L +
		2 * pow(a, 2)*(2 * pow(b, 2) - pow(c, 2) + 6 * b*d +
		4 * pow(d, 2))*L + 2 * pow(b, 2)*
		(4 * pow(c, 2) + 9 * pow(d, 2))*L - 8 * pow(A, 3)*c*u0 +
		8 * a*pow(A, 3)*(u0 - u1) + 8 * pow(A, 3)*c*u1 -
		8 * pow(A, 3)*d*v0 + 4 * b*
		(3 * pow(c, 2)*d*L + 3 * pow(d, 3)*L + 2 * pow(A, 3)*(v0 - v1))\
		+ 8 * pow(A, 3)*d*v1) / (8.*pow(A, 4)*y);
	hessian[2 * j][2 * j + 1] += (-2 * L*y*(3 * (u0 - u1)*(pow(xx, 2) + pow(y, 2))*
		(pow(u0 - u1, 2) + pow(v0 - v1, 2))*(v0 - v1) -
		2 * L*xx*(-(u2*pow(v0 - v1, 3)) +
		2 * pow(u1, 2)*u2*(-v0 + v1) +
		u0*(4 * u1*u2*(v0 - v1) +
		pow(u1, 2)*(5 * v0 - 2 * v1 - 3 * v2) +
		pow(v0 - v1, 2)*(3 * v0 - v1 - 2 * v2)) -
		2 * u1*pow(v0 - v1, 2)*(v0 - v2) +
		pow(u0, 3)*(3 * v0 - 2 * v1 - v2) + pow(u1, 3)*(-v0 + v2) +
		pow(u0, 2)*(2 * u2*(-v0 + v1) + u1*(-7 * v0 + 4 * v1 + 3 * v2)))\
		+ pow(L, 2)*(u0*((v0 - v1)*
		(pow(u2, 2) + (v0 - v2)*(3 * v0 - 2 * v1 - v2)) +
		2 * u1*u2*(3 * v0 - v1 - 2 * v2) + 2 * pow(u1, 2)*(v0 - v2))\
		- u1*(v0 - v1)*(pow(u2, 2) + pow(v0 - v2, 2)) +
		pow(u0, 3)*(3 * v0 - v1 - 2 * v2) -
		2 * u2*pow(v0 - v1, 2)*(v0 - v2) +
		2 * pow(u1, 2)*u2*(-v0 + v2) +
		pow(u0, 2)*(2 * u2*(-2 * v0 + v1 + v2) +
		u1*(-5 * v0 + v1 + 4 * v2))))) / pow(A, 4);

	//hessian[2 * j][2 * j + 1] = hessian[2 * j + 1][2 * j];

	//Dv2
	hessian[2 * j + 1][2 * j + 1] += (3 * pow(a, 4)*L + pow(b, 4)*L + 12 * pow(a, 3)*c*L +
		3 * pow(c, 4)*L + 4 * pow(c, 2)*pow(d, 2)*L + pow(d, 4)*L +
		2 * pow(a, 2)*(2 * pow(b, 2) + 9 * pow(c, 2) + 4 * pow(d, 2))*L +
		pow(b, 2)*(8 * pow(c, 2)*L - 2 * pow(d, 2)*L) -
		8 * pow(A, 3)*c*u0 + 4 * a*
		(3 * pow(b, 2)*c*L + 3 * pow(c, 3)*L + 3 * c*pow(d, 2)*L +
		2 * pow(A, 3)*(u0 - u1)) + 8 * pow(A, 3)*c*u1 -
		8 * pow(A, 3)*d*v0 + 8 * pow(A, 3)*b*(v0 - v1) +
		8 * pow(A, 3)*d*v1) / (8.*pow(A, 4)*y);
}