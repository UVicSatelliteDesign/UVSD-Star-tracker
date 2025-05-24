#pragma once
#include <math.h>
#include <iostream>

#include <iomanip>
#include <string>
#include <io.h>
#include <fcntl.h>

template <typename T, unsigned int D>
struct fixed_array {
	T entries[D];
};


template <typename T, unsigned int D>
class gen_vec {
private:
	fixed_array<T, D> components;
public:
	gen_vec(){
		for (int i = 0; i < D; i++) {
			components.entries[i] = T(0);
		}
	}
	gen_vec(T s) {
		for (int i = 0; i < D; i++) {
			components.entries[i] = s;
		}
	}
	gen_vec(fixed_array<T, D> c) {
		components = c;
	}

	T& operator [](unsigned int i) {
		return components.entries[i];
	}

	float operator* (gen_vec v) {
		float val = 0.0;
		for (int i = 0; i < D; i++) {
			val += (*this)[i] * v[i];
		}
		return val;
	}
};

template <unsigned int D>
using vec = gen_vec<float, D>;

template <unsigned int D>
using int_vec = gen_vec<int, D>;

template <unsigned int R, unsigned int C>
class matrix {
private:
	fixed_array<vec<C>, R> rows;
public:
	matrix() {

	}
	matrix(fixed_array<fixed_array<float, C>, R> c) {
		for (int i = 0; i < R; i++) {
			rows.entries[i] = c.entries[i];
		}
	}
	vec<C>& operator [](unsigned int i) {
		return rows.entries[i];
	}
};



template <unsigned int D>
float dist_sq(vec<D> a, vec<D> b) {
	float dist = 0.0;
	for (int i = 0; i < D; i++) {
		float delta = a[i] - b[i];
		dist += delta * delta;
	}
	return dist;
}


template <unsigned int D>
vec<D> clamp_vector(vec<D> v, vec<D> min, vec<D> max) {
	for (int i = 0; i < D; i++) {
		v[i] -= min[i];
		v[i] /= max[i] - min[i];
	}
	return v;
}

template <unsigned int D>
vec<D> sub_vector(vec<D> a, vec<D> b) {
	vec<D> out;
	for (int i = 0; i < D; i++) {
		out[i] = a[i] - b[i];
	}
	return out;
}



vec<3> cross(vec<3> a, vec<3> b);

template <unsigned int D>
vec<D> scale_vec(vec<D> v, float s) {
	vec<D> out;
	for (int i = 0; i < D; i++) {
		out[i] = v[i] * s;
	}
	return out;
}

template <unsigned int D>
float length(vec<D> v) {
	return sqrt(dot(v, v));
}
template <unsigned int D>
vec<D> normalize(vec<D> v) {
	float length = sqrt(v * v);
	return scale_vec(v, 1.0 / length);
}

matrix<3, 3> rotation_matrix(vec<3> pivot, float angle);

template <unsigned int R, unsigned int C>
vec<R> mat_mult_vec(matrix<R, C> m, vec<C> v) {
	vec<R> out;
	for (int i = 0; i < R; i++) {
		float sum = 0.0;
		for (int j = 0; j < C; j++) {
			sum += m[i][j] * v[j];
		}
		out[i] = sum;
	}
	return out;
}

template <unsigned int D>
float unit_vec_arc_length(vec<D> a, vec<D> b) {
	//return the arc length between two unit vectors
	return acos(a * b);
}

template <unsigned int D>
vec<D> solve_system_of_equations(matrix<D, D> M, vec<D> augment) {
	vec<D> solution;

	//forward elimination
	for (int i = 0; i < D; i++) {
		for (int j = i + 1; j < D; j++) {
			float scaling_factor = -1.0 * M[j][i] / M[i][i];
			for (int k = i + 1; k < D; k++) {
				M[j][k] += scaling_factor * M[i][k];
			}
			augment[j] += scaling_factor * augment[i];
		}
	}

	//back substitution:
	solution[D-1] = augment[D - 1] / M[D-1][D-1];
	for (int i = D - 2; i >= 0; i--) {
		float sum = augment[i];
		for (int j = i + 1; j < D; j++) {
			sum -= M[i][j] * solution[j];
		}
		solution[i] = sum / M[i][i];
	}

	return solution;
}

template <unsigned int D>
void print_vector(vec<D> v) {
	std::cout << "[";
	for (int c = 0; c < D; c++) {
		std::cout << v[c];
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
		std::cout << v[c];
		if (c < D - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "]";
}

template <unsigned int R, unsigned int C>
void print_matrix(matrix<R, C> M) {
	std::wstring output = L"";
	for (int i = 0; i < R; i++) {
		//determine opening and closing characters:
		std::wstring open = L"\u2502";
		std::wstring close = L"\u2502";
		if (i == 0) {
			open = L'\u250C';
			close = L'\u2510';
		}
		else if (i == R - 1) {
			open = L'\u2514';
			close = L'\u2518';
		}
		output += open;

		for (int j = 0; j < C; j++) {
			output += std::to_wstring(M[i][j]);
			if (j != C - 1) {
				output += L"\t";
			}
		}
		output += close + L'\n';
	}
	_setmode(_fileno(stdout), _O_U16TEXT);
	std::wcout << output;
	//set back to normal so that printing works correctly in other function
	_setmode(_fileno(stdout), _O_TEXT);
}