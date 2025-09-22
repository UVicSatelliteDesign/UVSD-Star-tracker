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
struct vec {
	T components[D];

	vec(){
		for (int i = 0; i < D; i++) {
			components.entries[i] = T(0);
		}
	}
	vec(T s) {
		for (int i = 0; i < D; i++) {
			components.entries[i] = s;
		}
	}
	vec(fixed_array<T, D> c) {
		components = c;
	}

	T& operator [](unsigned int i) {
		return components.entries[i];
	}

	float operator* (vec v) {
		float val = 0.0;
		for (int i = 0; i < D; i++) {
			val += (*this)[i] * v[i];
		}
		return val;
	}
};


template <typename T, unsigned int R, unsigned int C>
struct mat {
	vec<T, C> components[R];

	vec<T, C>& operator [](unsigned int i) {
		return components[i];
	}
};




template <typename T, unsigned int D>
float dist_sq(vec<T, D> a, vec<T, D> b) {
	float dist = 0.0;
	for (int i = 0; i < D; i++) {
		float delta = (float)a[i] - (float)b[i];
		dist += delta * delta;
	}
	return dist;
}



template <typename T, unsigned int D>
vec<T, D> sub_vecs(vec<T, D> a, vec<T, D> b) {

	vec<T, D> out;
	for (int i = 0; i < D; i++) {
		out[i] = a[i] - b[i];
	}
	return out;
}


template <typename T>
vec<T, 3> cross(vec<T, 3> a, vec<T, 3> b);

template <typename T, unsigned int D>
vec<T, D> scale_vec(vec<T, D> v, T s) {
	vec<D> out;
	for (int i = 0; i < D; i++) {
		out[i] = v[i] * s;
	}
	return out;
}

template <typename T, unsigned int D>
float length(vec<T, D> v) {
	return sqrt(dot(v, v));
}
template <typename T, unsigned int D>
vec<T, D> normalize(vec<T, D> v) {
	float length = sqrt(v * v);
	return scale_vec(v, 1.0 / length);
}

template <typename T>
mat<T, 3, 3> rotation_mat(vec<T, 3> pivot, T angle);

template <typename T, unsigned int R, unsigned int C>
vec<T, R> mat_mult_vec(mat<T, R, C> m, vec<T, C> v) {
	vec<T, R> out;
	for (int i = 0; i < R; i++) {
		float sum = 0.0;
		for (int j = 0; j < C; j++) {
			sum += m[i][j] * v[j];
		}
		out[i] = sum;
	}
	return out;
}

template <typename T, unsigned int D>
float unit_vec_arc_length(vec<T, D> a, vec<T, D> b) {
	//return the arc length between two unit vectors
	return acos(a * b);
}

template <typename T, unsigned int D>
vec<T, D> solve_system_of_equations(mat<T, D, D> M, vec<T, D> augment) {
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

template <typename T, unsigned int D>
void print_vector(vec<T, D> v) {
	std::cout << "[";
	for (int c = 0; c < D; c++) {
		std::cout << v[c];
		if (c < D - 1) {
			std::cout << ", ";
		}
	}
	std::cout << "]";
}

template <typename T, unsigned int R, unsigned int C>
void print_matrix(mat<T, R, C> M) {
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