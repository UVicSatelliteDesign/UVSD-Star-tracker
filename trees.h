#pragma once
#include <intrin.h>
#include <algorithm>
#include "linear_algebra.h"
#include <iomanip>
#include <io.h>
#include <fcntl.h>

template <typename T>
struct key_z_index_pair {
	T* object;
	unsigned long long int z_index;
};


template <typename T, unsigned int D>
struct binary_node {
	//the key must have a z-index element
	key_z_index_pair<T> key;
	int_vec<D> components;
	unsigned long long int split_position;
	int split_direction;//the direction can be infered by the remainder after dividing by 3 0 for x, 1 for y and 2 for z, -1 for undivided. This can also be used to check if there are child nodes
	binary_node* left;
	binary_node* right;
	binary_node* parent;
	binary_node* sibling;
};

constexpr unsigned int count_leading_zero_compile_time(unsigned int Z) {
	unsigned int leading_zeros = 0;
	for (int i = 31; i >= 0; i--) {
		if (((Z >> i) & 1) != 0) {
			break;
		}
		leading_zeros++;
	}
	return leading_zeros;
}

template <unsigned int D>
unsigned long long int spread_bits(unsigned long long int x) {
	//determine the number of steps to take by finding the base 2 logarithm of the number of bits per component and round up
	constexpr unsigned int total_bits_used = 64 - 64 % D;
	constexpr unsigned int bit_count = total_bits_used / D;


	//the built in instrinsic function __lzcnt is sadly not recognized as a constant expression and therefore cannot be evaluted at compile time
	//that is why the custom implementation count_leading_zero_compile_time is used. It is much less efficient than __lzcnt, but only needs to be 
	//evaluated once per template initialiaztion as compile time so it does not matter.
	constexpr unsigned int step_count = 32 - count_leading_zero_compile_time(bit_count);
	unsigned long long int filter = 0;
	for (int b = 0; b < bit_count; b++) {
		filter += (unsigned long long int)(1) << b;
	}
	unsigned long long int masks[step_count];


	//generate bit masks:
	for (int i = 0; i < step_count; i++) {
		//Bits per cluster in this mask
		unsigned int cluster_size = 1 << (step_count - 1 - i);


		//Bits per group, including 0s
		unsigned int group_size = D * cluster_size;


		//set bits
		masks[i] = 0;
		unsigned int bits_used = 0;
		for (int b = 0; b < 64; b++) {
			//determine if the bit is a zero or a one:
			unsigned int bit = (b % group_size) < cluster_size;


			//check if all the bits have been used up:
			bit = bit && (bits_used < bit_count);


			//Increment the bit counter
			bits_used += bit;


			//Set the bit in its place
			masks[i] += (unsigned long long int)(bit) << b;
		}
		/*
		std::cout << "Mask " << i << ":\t";
		print_binary_int64(masks[i]);
		std::cout << "\n";
		*/
	}




	x = x & filter;


	for (int i = 0; i < step_count; i++) {
		x = (x | (x << ((D - 1) << (step_count - 1 - i)))) & masks[i];
	}
	/*
		start code :	0000000000000000000000000000000000000000000111111111111111111111 = x
		Shifted bits:	0000000000011111111111111111111100000000000000000000000000000000 = S = x << 32
		Combined bits:	0000000000011111111111111111111100000000000111111111111111111111 = C = S | x
		Final bits:		0000000000011111000000000000000000000000000000001111111111111111 = F = C & mask


		Start code:		0000000000011111000000000000000000000000000000001111111111111111
		Shifted bits:	0000000000000000000000000000000011111111111111110000000000000000 << 16
		Combined bits:	0000000000011111000000000000000011111111111111111111111111111111
		Final bits:		0000000000011111000000000000000011111111000000000000000011111111


		Start code:		0000000000011111000000000000000011111111000000000000000011111111
		Shifted bits:	0001111100000000000000001111111100000000000000001111111100000000 << 8
		Combined bits:	0001111100011111000000001111111111111111000000001111111111111111
		final bits:		0001000000001111000000001111000000001111000000001111000000001111


		Start code:		0001000000001111000000001111000000001111000000001111000000001111
		Shifted bits:	0000000011110000000011110000000011110000000011110000000011110000 << 4
		Combined bits:	0001000011111111000011111111000011111111000011111111000011111111
		Final bits:		0001000011000011000011000011000011000011000011000011000011000011


		Start code:		0001000011000011000011000011000011000011000011000011000011000011
		Shifted bits:	0100001100001100001100001100001100001100001100001100001100001100
		Combined bits:	0101001111001111001111001111001111001111001111001111001111001111
		Final bits:		0001001001001001001001001001001001001001001001001001001001001001


	*/
	return x;
}


template <unsigned int D>
unsigned long long int cluster_bits(unsigned long long int Z) {




	//determine the number of steps to take by finding the base 2 logarithm of the number of bits per component and round up
	constexpr unsigned int total_bits_used = 64 - 64 % D;
	constexpr unsigned int bit_count = total_bits_used / D;


	//the built in instrinsic function __lzcnt is sadly not recognized as a constant expression and therefore cannot be evaluted at compile time
	//that is why the custom implementation count_leading_zero_compile_time is used. It is much less efficient than __lzcnt, but only needs to be 
	//evaluated once per template initialiaztion as compile time so it does not matter.
	constexpr unsigned int step_count = 32 - count_leading_zero_compile_time(bit_count);
	unsigned long long int filter = 0;
	for (int b = 0; b < 64; b++) {
		filter += (unsigned long long int)((b % D) < 1) << b;
	}
	unsigned long long int masks[step_count];


	//generate bit masks:
	for (int i = 0; i < step_count; i++) {
		//Bits per cluster in this mask
		unsigned int cluster_size = 2 << i;


		//Bits per group, including 0s
		unsigned int group_size = D * cluster_size;


		//set bits
		masks[i] = 0;
		unsigned int bits_used = 0;
		for (int b = 0; b < 64; b++) {
			//determine if the bit is a zero or a one:
			unsigned int bit = (b % group_size) < cluster_size;


			//check if all the bits have been used up:
			bit = bit && (bits_used < bit_count);


			//Increment the bit counter
			bits_used += bit;


			//Set the bit in its place
			masks[i] += (unsigned long long int)(bit) << b;
		}
		/*
		std::cout << "Mask " << i << ":\t";
		print_binary_int64(masks[i]);
		std::cout << "\n";
		*/
	}


	//filter bits
	Z = Z & filter;


	unsigned int shift_factor = D - 1;
	for (int i = 0; i < step_count; i++) {
		Z = (Z | (Z >> shift_factor)) & masks[i];
		shift_factor *= 2;
	}


	/*
		D = 3 again:
		Start code:		0001001001001001001001001001001001001001001001001001001001001001	21 groups
		Shifted bits:	0000010010010010010010010010010010010010010010010010010010010010 >> 2
		Combined bits:	0001011011011011011011011011011011011011011011011011011011011011
		Final bits:		0001000011000011000011000011000011000011000011000011000011000011


		Start code:		0001000011000011000011000011000011000011000011000011000011000011	10 + 1 groups
		Shifted bits:	0000000100001100001100001100001100001100001100001100001100001100 >> 4
		Combined bits:	0001000111001111001111001100001111001111001111001111001111001111
		Final bits:		0001000000001111000000001111000000001111000000001111000000001111


		Start code:		0001000000001111000000001111000000001111000000001111000000001111	5 + 1 groups
		Shifted bits:	0000000000010000000011110000000011110000000011110000000011110000 >> 8
		Combined bits:	0001000000011111000011111111000011111111000011111111000011111111
		Final bits:		0000000000011111000000000000000011111111000000000000000011111111


		Start code:		0000000000011111000000000000000011111111000000000000000011111111	2 + 1 + 1 groups
		Shifted bits:	0000000000000000000000000001111100000000000000001111111100000000
		Combined bits:	0000000000000000000000000001111100000000000000001111111100000000 >> 16
		Final bits:		0000000000011111000000000000000000000000000000001111111111111111


		Start code:		0000000000011111000000000000000000000000000000001111111111111111
		Shifted bits:	0000000000000000000000000000000000000000000111110000000000000000
		Combined bits:	0000000000011111000000000000000000000000000111111111111111111111
		Final bits:		0000000000000000000000000000000000000000000111111111111111111111




		Example with dimension = 3:


		Start code:		0001001001001001001001001001001001001001001001001001001001001001 = Z
		Target bits:	0000001000001000001000001000001000001000001000001000001000001000 = T = Z & mask
		Shifted bits:	0000000010000010000010000010000010000010000010000010000010000010 = S = T >> dimension - 1 = 2
		Combined bits:	0001001011001011001011001011001011001011001011001011001011001011
		Final bits:		0001000011000011000011000011000011000011000011000011000011000011


		Start code:		0001000011000011000011000011000011000011000011000011000011000011
		Target bits:	0001000011000000000011000000000011000000000011000000000011000000
		Shifted bits:	0000000000001100000000001100000000001100000000001100000000001100
		Combined bits:	0001000011001111000011001111000011001111000011001111000011001111
		Final bits:		0001000000001111000000001111000000001111000000001111000000001111


		 Start code:	0001000000001111000000001111000000001111000000001111000000001111
		 Target bits:	0001000000000000000000001111000000000000000000001111000000000000
		 Shifted bits:	0000000000010000000000000000000011110000000000000000000011110000
		 Combined bits:	0001000000011111000000001111000011111111000000001111000011111111
		 Final bits:	0000000000011111000000000000000011111111000000000000000011111111


		 Start code:	0000000000011111000000000000000011111111000000000000000011111111
		 Target bits:	0000000000000000000000000000000011111111000000000000000000000000
		 Shifted bits:	0000000000000000000000000000000000000000000000001111111100000000
		 Final bits:	0000000000011111000000000000000000000000000000001111111111111111


		 Start code:	0000000000011111000000000000000000000000000000001111111111111111
		 Target bits:	0000000000011111000000000000000000000000000000000000000000000000
		 Shifted bits:	0000000000000000000000000000000000000000000111110000000000000000
		 Combined bits:	0000000000011111000000000000000000000000000111111111111111111111
		 Final bits:	0000000000000000000000000000000000000000000111111111111111111111


________________________________________________________________________________________________________________________________________________________


		Example with dimension = 4:


		Start code: 	0001000100010001000100010001000100010001000100010001000100010001 = Z
		Target bits:	0001000000010000000100000001000000010000000100000001000000010000 = T = Z & mask (mask depends on dimension number)
		Shifted bits:	0000001000000010000000100000001000000010000000100000001000000010 = S = T >> dimension - 1 = 3
		Combined bits:	0001001100010011000100110001001100010011000000110001001100010011 = C = S | Z
		Final bits:		0000001100000011000000110000001100000011000000110000001100000011 = F = C & mask;
													|
													V
		Start code:		0000001100000011000000110000001100000011000000110000001100000011 = Z' = F
		Target bits:	0000001100000000000000110000000000000011000000000000001100000000 = T' = Z' & mask
		Shifted bits:	0000000000001100000000000000110000000000000011000000000000001100 = S' = T' >> 2 * (dimension - 1) = 6
		Combined bits:	0000001100001111000000110000111100000011000011110000001100001111 = C' = S' | Z'
		Final bits:		0000000000001111000000000000111100000000000011110000000000001111 = F' = C' & mask
													|
													V
		Start code:		0000000000001111000000000000111100000000000011110000000000001111 = Z'' = F'
		Target bits:	0000000000001111000000000000000000000000000011110000000000000000 = T'' = Z'' & mask;
		Shifted bits:	0000000000000000000000001111000000000000000000000000000011110000 = S'' = T'' >> 4 * (dimension - 1) = 12
		Combined bits:	0000000000001111000000001111111100000000000011110000000011111111 = C'' = S'' | Z''
		Final bits:		0000000000000000000000001111111100000000000000000000000011111111 = F'' = C'' & mask
													|
													V
		Start code:		0000000000000000000000001111111100000000000000000000000011111111 = Z''' = F''
		Target bits:	0000000000000000000000001111111100000000000000000000000000000000 = T''' = Z''' & mask
		Shifted bits:	0000000000000000000000000000000000000000000000001111111100000000 = S''' = T''' >> 8 * (dimension -1) = 24
		Combined bits:	0000000000000000000000001111111100000000000000001111111111111111 = C''' = S''' | Z'''
		Final bits:		0000000000000000000000000000000000000000000000001111111111111111 = F''' = C''' & mask


________________________________________________________________________________________________________________________________________________________


		Example with dimension = 5:


		Start code:		0000000010000100001000010000100001000010000100001000010000100001 12 groupd
		Shifted bits:	0000000000001000010000100001000010000100001000010000100001000010
		Combined bits:	0000000010001100011000110001100011000110001100011000110001100011
		Final bits:		0000000000001100000000110000000011000000001100000000110000000011


		Start code:		0000000000001100000000110000000011000000001100000000110000000011 6 groups
		Final bits:		0000000000000000000011000000001100000000110000000011000000001100
		Combined bits:	0000000000001100000011110000001111000000111100000011110000001111
		Final bits:		0000000000000000000011110000000000000000111100000000000000001111


		Start code:		0000000000000000000011110000000000000000111100000000000000001111 3 groups
		Shifted bits:	0000000000000000000000000000000000001111000000000000000011110000
		Combined bits:	0000000000000000000011110000000000001111000000000000000011111111
		Final bits:		0000000000000000000011110000000000000000000000000000000011111111


		Start code:		0000000000000000000011110000000000000000000000000000000011111111
		Final bits:		0000000000000000000000000000000000000000000000000000111100000000
		Combined bits:	0000000000000000000011110000000000000000000000000000111111111111
		Final bits:		0000000000000000000000000000000000000000000000000000111111111111


________________________________________________________________________________________________________________________________________________________


		Example with dimension = 6:


		Start code:		0000000001000001000001000001000001000001000001000001000001000001 = Z
		Target bits:	0000000001000000000001000000000001000000000001000000000001000000 = T = Z & mask
		Shifted bits:	0000000000000010000000000010000000000010000000000010000000000010 = S = T >> dimension - 1 = 5
		Combined bits:	0000000001000011000001000011000001000011000001000011000001000011 = C = S | Z
		Final bits:		0000000000000011000000000011000000000011000000000011000000000011 = F = C & mask;
													|
													V
		Start code:		0000000000000011000000000011000000000011000000000011000000000011 = Z' = F
		Target bits:	0000000000000000000000000011000000000000000000000011000000000000 = T' = Z' & mask //this mask is a lil different
		Shifted bits:	0000000000000000000000000000000000001100000000000000000000001100 = S' = T' >> 2 * (dimension - 1) = 10
		Combined bits:	0000000000000011000000000011000000001111000000000011000000001111 = C' = S' | Z'
		Final bits:		0000000000000011000000000000000000001111000000000000000000001111 = F' = C' & mask
													|
													V
		Start code:		0000000000000011000000000000000000001111000000000000000000001111 = Z'' = F'
		Target bits:	0000000000000000000000000000000000001111000000000000000000000000 = T'' = Z'' & mask
		Shifted bits:	0000000000000000000000000000000000000000000000000000000011110000 = S'' = T'' >> 4 * (dimension - 1) = 20
		Combined bits:	0000000000000011000000000000000000001111000000000000000011111111 = C'' = S'' | Z''
		Final bits:		0000000000000011000000000000000000000000000000000000000011111111 = F'' = C'' & mask
													|
								Auxillary step: move straggling bit to the left
													|
													V
		Start code:		0000000000000011000000000000000000000000000000000000000011111111 = Z''' = F''
		Target bits:	0000000000000011000000000000000000000000000000000000000000000000 = T''' = Z''' & mask
		Shifted bits:	0000000000000000000000000000000000000000000000000000001100000000 = S''' = T'''' >> 8 * (dimension - 1) = 40
		Combined bits:	0000000000000011000000000000000000000000000000000000001111111111 = C''' = S''' | Z'''
		Final bits:		0000000000000000000000000000000000000000000000000000001111111111 = F''' = C''' & mask
	*/
	return Z;
}

template <unsigned int D>
unsigned long long int clamped_vec_to_z_index(vec<D> v) {
	//assumes that all compoents of v are in the range 0.0 to 1.0
	unsigned long long int Z = 0;
	for (int i = 0; i < D; i++) {
		float x = v.components[i] * (1 << ((64 - (64 % D)) / D));
		Z += spread_bits<D>((unsigned long long int)(x)) << i;
	}
	return Z;
}
template <unsigned int D>
unsigned long long int vec_to_z_index(vec<D> v, vec<D> min, vec<D> max) {
	unsigned long long int Z = 0;
	for (int i = 0; i < D; i++) {
		Z += spread_bits<D>((unsigned long long int)(((v.components[i] - min.components[i]) / (max.components[i] - min.components[i])) * ((1 << ((64 - (64 % D)) / D)) - 1))) << i;
	}
	return Z;
}
template <unsigned int D>
int_vec<D> int_components_from_vec(vec<D> v, vec<D> min, vec<D> max) {
	int_vec<D> V;
	for (int i = 0; i < D; i++) {
		V.components[i] = ((v.components[i] - min.components[i]) / (max.components[i] - min.components[i])) * ((1 << ((64 - (64 % D)) / D)) - 1);
	}
	return V;
}
template <unsigned int D>
vec<D> z_index_to_vec(unsigned long long int Z, vec<D> min, vec<D> max) {
	vec<D> v;
	for (int i = 0; i < D; i++) {
		v.components[i] = cluster_bits<D>(Z) * (max.components[i] - min.components[i]) + min.components[i];
	}
	return v;
}
template <unsigned int D>
unsigned long long int extract_int_component(unsigned long long int Z, unsigned int k) {
	return cluster_bits<D>(Z >> k);//cluster_bits(Z >> 2 - k);
}


template <typename T, unsigned int D>
int find_split(key_z_index_pair<T>* sorted_pairs, unsigned int first, unsigned int last, int* split_direction, unsigned long long int* split_pos) {

	//count the number of bits the start and end values have in common:
	//the first non-zero bit (going from left to right) is the first bit that differs between the two values. The split point will be the index of the last value to have a zero in that spot
	unsigned long long int commonality = sorted_pairs[first].z_index ^ sorted_pairs[last].z_index;
	unsigned int bit_of_interest = _lzcnt_u64(commonality);
	*split_direction = (63 - bit_of_interest) % D;//0 for x, 1 for y, 2 for z
	unsigned long long int bitmask = 1ull << (63 - bit_of_interest);
	//use binary search to find the where this 0 becomes a 1:


	unsigned int index = first;
	//for the algorithm to work correctly, the step size needs to be a power of two. Using the leading zero count to efficiently approximate the base 2 logarithm of the range
	unsigned int step_size = 1 << (32 - __lzcnt(last - first));
	step_size /= 2;


	/*
	std::cout << "start: " << first << " end: " << last << ".\tMask:" << "\t";
	print_binary_int64(bitmask);
	std::cout << "\n";
	*/
	while (step_size >= 1) {
		/*
		std::cout << "index: " << index << " step size: " << step_size << "\t\t";
		print_binary_int64(z_indices[index].z_index);
		std::cout << "  ";
		print_binary_int64(z_indices[index + step_size].z_index);
		std::cout << "  ";
		print_binary_int64(filtered);
		std::cout << "\n";
		*/
		if (index + step_size < last) {
			unsigned long long int filtered = sorted_pairs[index + step_size].z_index & bitmask;
			if (filtered == 0ull) {
				index += step_size;
			}
		}
		step_size /= 2;
	}
	unsigned long long int split_plane = sorted_pairs[last].z_index & (0xffffffffffffffffull << (63 - bit_of_interest));
	*split_pos = extract_int_component<D>(split_plane, *split_direction);
	/*
	std::cout << "Split plane: ";
	print_binary_int64(split_plane);
	std::cout << " pos: " << *split_pos << "\n";
	std::cout << "Final:" << index << "\n\n";
	*/
	return index;
}
template <typename T, unsigned int D>
binary_node<T, D>* divide_cell(key_z_index_pair<T>* pairs, unsigned int first, unsigned int last) {
	unsigned long long int split_pos;
	int split_direction;
	int split_index = find_split<T, D>(pairs, first, last, &split_direction, &split_pos);
	unsigned long long int z_index = pairs[split_index].z_index;


	binary_node<T, D>* left = NULL;
	binary_node<T, D>* right = NULL;
	binary_node<T, D>* node = new binary_node<T, D>({ pairs[split_index],
		{},
		split_pos,
		-1,
		NULL, NULL, NULL, NULL });
	for (int i = 0; i < D; i++) {
		node->components.components[i] = extract_int_component<D>(z_index, i);
	}
	if (last - first > 0) {
		left = divide_cell<T, D>(pairs, first, split_index);
		right = divide_cell<T, D>(pairs, split_index + 1, last);
		left->parent = node;
		right->parent = node;
		left->sibling = right;
		right->sibling = left;
		node->split_direction = split_direction;
	}
	node->left = left;
	node->right = right;
	return node;
}


template <typename T>
bool key_z_index_pair_compare(key_z_index_pair<T> a, key_z_index_pair<T> b) {
	return a.z_index < b.z_index;
}
template <typename T, unsigned int D>
binary_node<T, D>* generate_binary_tree(T* objects, unsigned int object_count, unsigned long long int (*z_index_extractor)(T)) {
	//generate pairs:
	key_z_index_pair<T>* pairs = new key_z_index_pair<T>[object_count];

	//initialize pairs:
	for (int i = 0; i < object_count; i++) {
		pairs[i].object = objects + i;
		pairs[i].z_index = z_index_extractor(objects[i]);
	}

	//sort pairs by z_index;
	std::sort(pairs, pairs + object_count, key_z_index_pair_compare<T>);


	//create the tree:
	binary_node<T, D>* root = divide_cell<T, D>(pairs, 0, object_count - 1);


	delete[] pairs;


	return root;
}


template <unsigned int D>
unsigned long long int dist_squared_between_z_indices(unsigned long long int a, unsigned long long int b) {
	unsigned long long int s = 0;


	for (int i = 0; i < D; i++) {
		unsigned long long int c1 = extract_int_component<D>(a, i);
		unsigned long long int c2 = extract_int_component<D>(b, i);;
		unsigned long long int delta;
		if (c1 < c2) {
			delta = c2 - c1;
		}
		else {
			delta = c1 - c2;
		}
		s += delta * delta;
	}


	return s;
}



template <typename T, unsigned int D>
T** brute_force_neighbors(T* objects, unsigned int object_count, binary_node<T, D>* reference, unsigned int neighbor_count, unsigned long long int (*z_index_extractor)(T)) {
	unsigned long long int* distances = new unsigned long long int[neighbor_count];
	T** close_objects = new T * [neighbor_count];
	for (int i = 0; i < neighbor_count; i++) {
		distances[i] = 0xffffffffffffffff;
	}
	for (int i = 0; i < object_count; i++) {
		if (objects + i != reference->key.object) {
			unsigned long long int dist = dist_squared_between_z_indices<D>(z_index_extractor(objects[i]), reference->key.z_index);//dist_sq(reference->key.object->dir, obj[i].dir);


			if (dist < distances[neighbor_count - 1]) {
				distances[neighbor_count - 1] = dist;
				close_objects[neighbor_count - 1] = objects + i;
				int current_k = neighbor_count - 2;
				while (current_k >= 0) {
					if (dist < distances[current_k]) {
						distances[current_k + 1] = distances[current_k];
						close_objects[current_k + 1] = close_objects[current_k];
						distances[current_k] = dist;
						close_objects[current_k] = objects + i;
					}
					current_k--;
				}
			}
		}
	}
	return close_objects;
}


template <typename T, unsigned int D>
long long int split_plane_dist(binary_node<T, D>* node, int_vec<D> reference) {
	long long int split_plane = node->split_position;//*(reinterpret_cast<float*>(node->position) + node->split_direction);
	long long int point_pos = reference.components[node->split_direction];//*(reinterpret_cast<float*>(&reference->key->dir) + node->split_direction);
	long long int delta = point_pos - split_plane;
	return delta;
}
template <typename T, unsigned int D>
unsigned long long int dist_to_node(binary_node<T, D>* a, int_vec<D> v) {
	unsigned long long int distance_squared = 0;
	for (int i = 0; i < D; i++) {
		unsigned long long int c1 = a->components.components[i];
		unsigned long long int c2 = v.components[i];
		unsigned long long int delta;
		if (c1 < c2) {
			delta = c2 - c1;
		}
		else {
			delta = c1 - c2;
		}
		distance_squared += delta * delta;
	}
	return distance_squared;
}




template <typename T, unsigned int D>
void inspect_node(binary_node<T, D>* node, binary_node<T, D>** neighbors, unsigned long long int* distances, int_vec<D> reference, unsigned int k) {
	/*
	std::cout << "Visited child node: ";
	print_vector<D>(*node->key.object);
	std::cout << "\t" << node->split_direction << "\n";
	*/
	if (node->split_direction != -1) {
		//find distance to split plane:
		long long int dist = split_plane_dist(node, reference);


		if (dist * dist < distances[k - 1]) {
			//if the distance from the reference node to the split plane is closer than the current estimate for the furthest point
			//then than both children could be potential candidate
			inspect_node(node->left, neighbors, distances, reference, k);
			inspect_node(node->right, neighbors, distances, reference, k);


		}
		else {
			//The space occupied by this node does intersect the current furthest point of interest bubble,
			//but since the distance between this node's split plane and the reference node is too far, only one of the
			//two children is a viable candidate
			if (dist > 0) {
				inspect_node(node->right, neighbors, distances, reference, k);
			}
			else {
				inspect_node(node->left, neighbors, distances, reference, k);
			}
		}
	}
	else {


		unsigned long long int distance_squared = dist_to_node(node, reference);


		if (distance_squared < distances[k - 1]) {
			distances[k - 1] = distance_squared;
			neighbors[k - 1] = node;
			int current_k = k - 2;
			while (current_k >= 0) {
				if (distance_squared < distances[current_k]) {
					distances[current_k + 1] = distances[current_k];
					neighbors[current_k + 1] = neighbors[current_k];
					distances[current_k] = distance_squared;
					neighbors[current_k] = node;
				}
				current_k--;
			}
		}
	}
}
template <typename T, unsigned int D>
binary_node<T, D>** find_k_nearest_neighbors(binary_node<T, D>* host, int_vec<D> reference, unsigned int k, int log) {
	binary_node<T, D>** neighbors = new binary_node<T, D> *[k];


	//maintain ascending order. Also these are sqaured distances
	unsigned long long int* distances = new unsigned long long int[k];


	//Set the distances to the largest possible values:
	for (int i = 0; i < k; i++) {
		distances[i] = 0xffffffffffffffff;
	}


	/*	start at the reference node. The reason for starting here is twofold :
	*	1:	Starting at the reference node eliminates the possibility of counting
	*		itself as the closest point, as it would have a distance of 0, which is unbeatable
	*
	*	2:	By starting at the reference node, the closest neightbors are investigated first,
	*		increasing the likely hood of finding most of the close neighbors early on
	*		which allows for a greater portion of the tree to be pruned during traversal
	*/
	binary_node<T, D>* current = host;
	while (current->sibling) {
		//check the distance between this node and the plane separating it from its sibling:
		unsigned long long int dist = split_plane_dist<T>(current->parent, reference);//You can trust that this node has a parant because it has a sibling


		if (dist * dist < distances[k - 1]) {
			/*
			std::cout << "Visited parent node: ";
			print_vector<D>(*current->key.object);
			std::cout << "\t" << current->split_direction << "\n";
			*/
			inspect_node<T, D>(current->sibling, neighbors, distances, reference, k);


		}




		current = current->parent;
	}
	std::sort(distances, distances + k);
	if (log) {
		for (int i = 0; i < k; i++) {
			std::cout << distances[i] << "\n";
		}
	}
	return neighbors;
}
template <typename T, unsigned int D>
binary_node<T, D>* find_cell(binary_node<T, D>* root, vec<D> v, vec<D> min, vec<D> max) {
	//compute the integer components of the vector:
	int_vec<D> int_components;
	for (int i = 0; i < D; i++) {
		int_components.components[i] = (unsigned long long int)(((v.components[i] - min.components[i]) / (max.components[i] - min.components[i])) * (1 << ((64 - (64 % D)) / D)));
	}


	//check root node
	binary_node<T, D>* current = root;
	while (current->split_direction != -1) {
		//check which side of the split plane this node is one:
		unsigned long long int split_plane = current->split_position;
		if (int_components.components[current->split_direction] < split_plane) {
			current = current->left;
		}
		else {
			current = current->right;
		}
	}
	return current;
}
template <typename T, unsigned int k>
void rake_leaf_nodes(binary_node<T, k>* node, binary_node<T, k>** leaf_node_array, unsigned int* leaf_index) {
	if (node->split_direction == -1) {
		leaf_node_array[(*leaf_index)++] = node;
	}
	else {
		rake_leaf_nodes(node->left, leaf_node_array, leaf_index);
		rake_leaf_nodes(node->right, leaf_node_array, leaf_index);
	}
}

template <typename T, unsigned int k>
binary_node<T, k>** collect_leaf_nodes(binary_node<T, k>* root, unsigned int expected_node_count) {
	binary_node<T, k>** leaf_nodes = new binary_node<T, k>*[expected_node_count];
	unsigned int node_index = 0;
	rake_leaf_nodes(root, leaf_nodes, &node_index);
	return leaf_nodes;
}


template <typename T, unsigned int k>
void delete_branch(binary_node<T, k>* node) {
	if (node->split_direction != -1) {
		delete_branch(node->left);
		delete_branch(node->right);
	}
	delete node;
}



template <typename T, unsigned  D>
void print_tree(binary_node<T, D>* node, void (*print_object)(T) = NULL) {
	//climb up the tree
	binary_node<T, D>* current = node;
	std::wstring output = L"";
	if (current->parent) {
		//right sibling
		if (current->sibling->key.z_index < current->key.z_index) {
			output = L"\u2514\u2500\u2500\u2500\u2500" + output;
		}
		else {
			//left sibling
			output = L"\u251C\u2500\u2500\u2500\u2500" + output;
		}
		current = current->parent;
	}
	while (current) {
		//check if this node is a left or right child


		if (current->sibling) {
			if (current->sibling->key.z_index < current->key.z_index) {
				//std::cout << "   ";
				output = L"     " + output;
			}
			else {
				//std::cout << "|  ";
				output = L"\u2502    " + output;
			}
		}
		current = current->parent;
	}

	//set to U16 to properly dispaly box drawing characters
	_setmode(_fileno(stdout), _O_U16TEXT);
	std::wcout << output;
	//set back to normal so that printing works correctly in other function
	_setmode(_fileno(stdout), _O_TEXT);


	//std::cout << node->key->name;//print_star(*node->star);
	if (print_object) {
		print_object(*node->key.object);
	}


	std::cout << " dir: " << node->split_direction << " pos: " << node->split_position << "\t";
	print_int_vector<D>(node->components);
	std::cout << "\n";
	//check if there are child nodes
	if (node->split_direction != -1) {
		print_tree(node->left, print_object);
		print_tree(node->right, print_object);
	}
}
void print_binary_int64(unsigned long long int N);