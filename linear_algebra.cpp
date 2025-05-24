#include "linear_algebra.h"

vec<3> cross(vec<3> a, vec<3> b) {
	vec<3> v;
	v[0] = a[1] * b[2] - a[2] * b[1];
	v[1] = a[2] * b[0] - a[0] * b[2];
	v[2] = a[0] * b[1] - a[1] * b[0];
	return v;
}

matrix<3, 3> rotation_matrix(vec<3> pivot, float angle) {
	float s = sin(angle);
	float c = cos(angle);
	return matrix<3, 3>({ {
		{(pivot[0] * pivot[0]) * (1 - c) + c,							(pivot[0] * pivot[1]) * (1 - c) - pivot[2] * s,	(pivot[0] * pivot[2]) * (1 - c) + pivot[1] * s},
		{(pivot[0] * pivot[1]) * (1 - c) + pivot[2] * s,	(pivot[1] * pivot[1]) * (1 - c) + c,							(pivot[1] * pivot[2]) * (1 - c) - pivot[0] * s},
		{(pivot[0] * pivot[2]) * (1 - c) - pivot[1] * s,	(pivot[1] * pivot[2]) * (1 - c) + pivot[0] * s,	(pivot[2] * pivot[2]) * (1 - c) + c}
	} });
}