#pragma once
#include <math.h>
#include <iostream>

template <unsigned int D>
struct vec {
	float components[D];
};

template <unsigned int R, unsigned int C>
struct matrix {
	float components[R][C];
};

template <unsigned int D>
struct int_vec {
	unsigned long long int components[D];
};

template <unsigned int D>
float dist_sq(vec<D> a, vec<D> b) {
	float dist = 0.0;
	for (int i = 0; i < D; i++) {
		float delta = a.components[i] - b.components[i];
		dist += delta * delta;
	}
	return dist;
}


template <unsigned int D>
vec<D> clamp_vector(vec<D> v, vec<D> min, vec<D> max) {
	for (int i = 0; i < D; i++) {
		v.components[i] -= min.components[i];
		v.components[i] /= max.components[i] - min.components[i];
	}
	return v;
}


template <unsigned int D>
float dot(vec<D> a, vec<D> b) {
	float val = 0.0;
	for (int i = 0; i < D; i++) {
		val += a.components[i] * b.components[i];
	}
	return val;
}

vec<3> cross(vec<3> a, vec<3> b);

template <unsigned int D>
vec<D> scale_vec(vec<D> v, float s) {
	vec<D> out;
	for (int i = 0; i < D; i++) {
		out.components[i] = v.components[i] * s;
	}
	return out;
}

template <unsigned int D>
vec<D> normalize(vec<D> v) {
	float length = sqrt(dot(v, v));
	return scale_vec(v, 1.0 / length);
}

matrix<3, 3> rotation_matrix(vec<3> pivot, float angle);

template <unsigned int R, unsigned int C>
vec<R> mat_mult_vec(matrix<R, C> m, vec<C> v) {
	vec<R> out;
	for (int i = 0; i < R; i++) {
		float sum = 0.0;
		for (int j = 0; j < C; j++) {
			sum += m.components[i][j] * v.components[j];
		}
		out.components[i] = sum;
	}
	return out;
}

template <unsigned int D>
float unit_vec_arc_length(vec<D> a, vec<D> b) {
	//return the arc length between two unit vectors
	return acos(dot<D>(a, b));
}


template <unsigned int D>
void print_vector(vec<D> v) {
	std::cout << "[";
	for (int c = 0; c < D; c++) {
		std::cout << v.components[c];
		if (c < D - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "]";
}
template <unsigned int D>
void print_int_vector(int_vec<D> v) {
	std::cout << "[";
	for (int c = 0; c < D; c++) {
		std::cout << v.components[c];
		if (c < D - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "]";
}