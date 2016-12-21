#pragma once

class iso_pre
{
public:
	double Pt, Pt2; //L * y / 8.0f
	double ld21, ld31, ld32;
	double p3dotp2;

	double cot[3];

	double area_uv;
	double p1,p2;
	// scott precomp
	double c00, c01, c02, c11, c12, c22;
	double gradlylyu0, gradlylyv0, gradlylyu1, gradlylyv1, gradlylyu2, gradlylyv2;
	double ly;
	double num; // temp for value of quadratic numerator
	double denom; // temp for value of denominator in evaluation, stored as recip
	double denomQuad, denomNonlin;
	double L, x, y;

	iso_pre();
	iso_pre(double L, double x, double y);

};