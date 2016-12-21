#include "vect.h"

vect4d operator*(const int a, const vect4d v){return v*a;}
vect3d operator*(const int a, const vect3d v){return v*a;}
vect2d operator*(const int a, const vect2d v){return v*a;}

vect4d operator*(const float a, const vect4d v){return v*a;}
vect3d operator*(const float a, const vect3d v){return v*a;}
vect2d operator*(const float a, const vect2d v){return v*a;}

vect4d operator*(const double a, const vect4d v){return v*a;}
vect3d operator*(const double a, const vect3d v){return v*a;}
vect2d operator*(const double a, const vect2d v){return v*a;}

vect4f operator*(const int a, const vect4f v){return v*a;}
vect3f operator*(const int a, const vect3f v){return v*a;}
vect2f operator*(const int a, const vect2f v){return v*a;}

vect4f operator*(const float a, const vect4f v){return v*a;}
vect3f operator*(const float a, const vect3f v){return v*a;}
vect2f operator*(const float a, const vect2f v){return v*a;}

vect4f operator*(const double a, const vect4f v){return v*a;}
vect3f operator*(const double a, const vect3f v){return v*a;}
vect2f operator*(const double a, const vect2f v){return v*a;}