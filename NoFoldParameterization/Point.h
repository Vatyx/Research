#ifndef POINT_H
#define POINT_H

class Point
{
public:
	double x, y;

	Point ( void ) {}

	Point ( double nx, double ny ) : x ( nx ), y ( ny ) {}

	Point operator + ( const Point &p ) const
	{
		Point rvalue = *this;
		rvalue.x += p.x;
		rvalue.y += p.y;

		return rvalue;
	}

	Point operator - ( const Point &p ) const
	{
		Point rvalue = *this;
		rvalue.x -= p.x;
		rvalue.y -= p.y;

		return rvalue;
	}

	Point operator / ( const Point &p ) const
	{
		Point rvalue = *this;

		rvalue.x /= p.x;
		rvalue.y /= p.y;

		return rvalue;
	}

	Point operator * ( double s ) const
	{
		Point rvalue = *this;

		rvalue.x *= s;
		rvalue.y *= s;

		return rvalue;
	}
};

#endif