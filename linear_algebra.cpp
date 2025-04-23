#include "linear_algebra.h"

vec<3> cross(vec<3> a, vec<3> b) {
	return { a.components[1] * b.components[2] - a.components[2] * b.components[1], a.components[2] * b.components[0] - a.components[0] * b.components[2], a.components[0] * b.components[1] - a.components[1] * b.components[0] };
}

matrix<3, 3> rotation_matrix(vec<3> pivot, float angle) {
	float s = sin(angle);
	float c = cos(angle);
	return { {
		{(pivot.components[0] * pivot.components[0]) * (1 - c) + c,							(pivot.components[0] * pivot.components[1]) * (1 - c) - pivot.components[2] * s,	(pivot.components[0] * pivot.components[2]) * (1 - c) + pivot.components[1] * s},
		{(pivot.components[0] * pivot.components[1]) * (1 - c) + pivot.components[2] * s,	(pivot.components[1] * pivot.components[1]) * (1 - c) + c,							(pivot.components[1] * pivot.components[2]) * (1 - c) - pivot.components[0] * s},
		{(pivot.components[0] * pivot.components[2]) * (1 - c) - pivot.components[1] * s,	(pivot.components[1] * pivot.components[2]) * (1 - c) + pivot.components[0] * s,	(pivot.components[2] * pivot.components[2]) * (1 - c) + c}
	} };
}