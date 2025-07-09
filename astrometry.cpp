#include "astrometry.h"
#include "database.h"
#include <random>
#include "linear_algebra.h"
#include "trees.h"
#include <string>
#define PI 3.14159265359


vec<3> point_on_sphere(float theta, float phi) {
	vec<3> v;
	v[0] = cos(theta) * cos(phi);
	v[1] = sin(phi);
	v[2] = sin(theta) * cos(phi);


	return v;
}
float random_float(int seed) {
	//in the future it might be worth while to use a global seed and a single mersenne twiser instead of generating one per number


	std::mt19937 gen(seed);

	std::uniform_int_distribution<unsigned int> distrib(0, 0xffffffff);


	return float(distrib(gen)) / float(0xffffffff);
}

star* generate_random_stars(unsigned int star_count) {
	unsigned int seed = 0;


	star* stars = new star[star_count];
	for (int i = 0; i < star_count; i++) {
		float ra = random_float(seed++) * 2.0 * PI;
		float dec = asin(random_float(seed++) * 2.0 - 1.0);

		stars[i].dir = point_on_sphere(ra, dec);
	}
	return stars;
}
int compare_stars_simple(star a, star b) {
	int same = true;
	for (int i = 0; i < 3; i++) {
		if (a.dir[i] != b.dir[i]) {
			same = false;
			break;
		}
	}
	return same;
}
void find_star_neighbors(binary_node<star, 3>** stars, unsigned int star_count) {
	feature* stellar_web = new feature[star_count];
	for (int i = 0; i < star_count; i++) {
		binary_node<star, 3>** close_stars = find_k_nearest_neighbors<star, 3>(stars[i], stars[i]->components, 3, false);
		float distances[6];
		distances[0] = unit_vec_arc_length(stars[i]->key.object->dir, close_stars[0]->key.object->dir);
		distances[1] = unit_vec_arc_length(stars[i]->key.object->dir, close_stars[1]->key.object->dir);
		distances[2] = unit_vec_arc_length(stars[i]->key.object->dir, close_stars[2]->key.object->dir);
		distances[3] = unit_vec_arc_length(close_stars[0]->key.object->dir, close_stars[1]->key.object->dir);
		distances[4] = unit_vec_arc_length(close_stars[1]->key.object->dir, close_stars[2]->key.object->dir);
		distances[5] = unit_vec_arc_length(close_stars[2]->key.object->dir, close_stars[0]->key.object->dir);
		float longest = 0.0;
		for (int j = 0; j < 6; j++) {
			if (distances[j] > longest) {
				longest = distances[j];
			}
		}


		for (int c = 0; c < 6; c++) {
			stars[i]->key.object->primary.distances[c] = distances[c];
		}
		delete[] close_stars;
	}
}
feature* find_star_neighbors_redundancy(binary_node<star, 3>** stars, unsigned int star_count) {
	feature* stellar_web = new feature[star_count * 4];
	int s = 0;
	for (int i = 0; i < star_count; i++) {
		binary_node<star, 3>** close_stars = find_k_nearest_neighbors<star, 3>(stars[i], stars[i]->components, 6, false);


		int neighbors[4][3] = {
			{0, 1, 2},
			{1, 2, 3},
			{0, 2, 3},
			{0, 1, 3}
			/*
			{0, 3, 4},
			{1, 3, 4},
			{2, 3, 4},
			{3, 4, 5}
			*/
		};
		for (int j = 0; j < 4; j++) {
			float distances[6];
			distances[0] = unit_vec_arc_length(stars[i]->key.object->dir, close_stars[neighbors[j][0]]->key.object->dir);
			distances[1] = unit_vec_arc_length(stars[i]->key.object->dir, close_stars[neighbors[j][1]]->key.object->dir);
			distances[2] = unit_vec_arc_length(stars[i]->key.object->dir, close_stars[neighbors[j][2]]->key.object->dir);
			distances[3] = unit_vec_arc_length(close_stars[neighbors[j][0]]->key.object->dir, close_stars[neighbors[j][1]]->key.object->dir);
			distances[4] = unit_vec_arc_length(close_stars[neighbors[j][1]]->key.object->dir, close_stars[neighbors[j][2]]->key.object->dir);
			distances[5] = unit_vec_arc_length(close_stars[neighbors[j][2]]->key.object->dir, close_stars[neighbors[j][0]]->key.object->dir);
			float longest = 0.0;
			for (int k = 0; k < 6; k++) {
				if (distances[k] > longest) {
					longest = distances[k];
				}
			}


			for (int c = 0; c < 6; c++) {
				stellar_web[s].distances[c] = distances[c];
			}
			s++;
		}
		delete[] close_stars;
	}
	return stellar_web;
}
unsigned long long int z_index_from_star(star s) {
	return  vec_to_z_index(s.dir, vec<3>(-1.0f), vec<3>(1.0f));
}



star** find_matches(binary_node<star, 6>* root, feature reference, unsigned int quad_count, int match_count) {


	int_vec<6> int_components = int_components_from_vec<6>(reference.distances, vec<6>(0.0), vec<6>(1.0));


	//construct star quad tree:
	//print_tree(root, print_star_quad);


	//find the cell which the reference quad falls into
	binary_node<star, 6>* host_cell = find_cell(root, reference.distances, vec<6>(0.0), vec<6>(1.0));


	//the star quad corresponding to the host cell is not nessessarily the most similar so nearest neighbor search is used:
	binary_node<star, 6>** nearby_quad_nodes = find_k_nearest_neighbors<star, 6>(host_cell, int_components, match_count - 1, false);

	star** possible_stars = new star*[match_count];
	for (int i = 0; i < match_count - 1; i++) {
		possible_stars[i] = nearby_quad_nodes[i]->key.object;
	}
	possible_stars[match_count - 1] = host_cell->key.object;

	return possible_stars;
}

feature* generate_star_quads_from_star_centroids(vec<2>* centroids, unsigned int star_count) {
	feature* star_quads = new feature[star_count];
	for (int i = 0; i < star_count; i++) {
		float close_distances[3] = { INFINITY, INFINITY, INFINITY };
		vec<2>* close_stars[3];
		for (int j = 0; j < star_count; j++) {
			if (i != j) {
				float dist = dist_sq(centroids[i], centroids[j]);
				if (dist < close_distances[2]) {
					close_distances[2] = dist;
					close_stars[2] = centroids + j;
					if (close_distances[2] < close_distances[1]) {
						close_distances[2] = close_distances[1];
						close_stars[2] = close_stars[1];
						close_distances[1] = dist;
						close_stars[1] = centroids + j;


						if (close_distances[1] < close_distances[0]) {
							close_distances[1] = close_distances[0];
							close_stars[1] = close_stars[0];
							close_distances[0] = dist;
							close_stars[0] = centroids + j;
						}
					}
				}
			}
		}

		float distances[6];
		distances[0] = sqrt(dist_sq(centroids[i], *close_stars[0]));
		distances[1] = sqrt(dist_sq(centroids[i], *close_stars[1]));
		distances[2] = sqrt(dist_sq(centroids[i], *close_stars[2]));
		distances[3] = sqrt(dist_sq(*close_stars[0], *close_stars[1]));
		distances[4] = sqrt(dist_sq(*close_stars[1], *close_stars[2]));
		distances[5] = sqrt(dist_sq(*close_stars[2], *close_stars[0]));


		float longest = 0.0;
		for (int j = 0; j < 6; j++) {
			if (distances[j] > longest) {
				longest = distances[j];
			}
		}
		for (int c = 0; c < 6; c++) {
			star_quads[i].distances[c] = distances[c] / longest;
		}
	}
	return star_quads;
}

vec<2>* generate_synthetic_image_data(vec<3> normal, float ccw_rotation, float fov, star* stars, unsigned int star_count, unsigned int* visible_star_count, star*** visible_stars_out, float noise, int random_seed) {
	vec<3> reference_up(0.0);
	reference_up[1] = 1.0;

	vec<3> right = normalize(cross(normal, reference_up));
	//apply twist:
	matrix<3, 3> ccw_twist_matrix = rotation_matrix(normal, ccw_rotation);
	right = normalize(mat_mult_vec(ccw_twist_matrix, right));


	vec<3> up = normalize(cross(right, normal));


	matrix<3, 3> pitch_up_matrix = rotation_matrix(right, fov / 2);
	matrix<3, 3> pitch_down_matrix = rotation_matrix(right, -fov / 2);
	matrix<3, 3> yaw_left_matrix = rotation_matrix(up, fov / 2);
	matrix<3, 3> yaw_right_matrix = rotation_matrix(up, -fov / 2);


	//determine the 4 cutting plane normals:
	vec<3> top_plane = normalize(mat_mult_vec(pitch_up_matrix, up));
	vec<3> bottom_plane = normalize(mat_mult_vec(pitch_down_matrix, scale_vec(up, -1.0)));
	vec<3> left_plane = normalize(mat_mult_vec(yaw_left_matrix, scale_vec(right, -1.0)));
	vec<3> right_plane = normalize(mat_mult_vec(yaw_right_matrix, right));


	/*
	std::cout << "nornal vector: ";
	print_vector(normal);
	std::cout << "\tup vector: ";
	print_vector(up);
	std::cout << "\tright vector: ";
	print_vector(right);


	std::cout << "\ntop plane: ";
	print_vector(top_plane);
	std::cout << "\tbottom plane: ";
	print_vector(bottom_plane);
	std::cout << "\tleft plane: ";
	print_vector(left_plane);
	std::cout << "\tright plane: ";
	print_vector(right_plane);
	std::cout << "\n";
	*/
	//create array of stars in the field of view of the camera
	star** visible_stars = new star * [star_count];

	//create array of output vectors:
	vec<2>* star_centroids = new vec<2>[star_count];

	//Keep track of how many visible stars are in the image
	*visible_star_count = 0;
	int seed = 0;

	float normalization_factor = tan(fov / 2);

	//random 
	std::default_random_engine generator;
	generator.seed(random_seed);
	std::normal_distribution<float> guassian(0.0, noise);

	//increment over the aray of stars and pick out the ones which are within the field of view
	for (unsigned int i = 0; i < star_count; i++) {
		vec<3> star_norm = stars[i].dir;


		//project star direction onto normal axis and get distance:
		float z = stars[i].dir * normal;


		//get horizontal coordinate:
		float x = right * stars[i].dir / (z * normalization_factor) + guassian(generator);


		//get vertical coordinate:
		float y = up * stars[i].dir / (z * normalization_factor) + guassian(generator);


		if (z > 0.0 && abs(x) < 1.0 && abs(y) < 1.0) {
			visible_stars[*visible_star_count] = stars + i;
			star_centroids[*visible_star_count] = vec<2>({ x , y });
			(*visible_star_count)++;
		}
	}
	*visible_stars_out = visible_stars;


	return star_centroids;
}

unsigned long long int z_index_from_star_quad(star s) {
	return vec_to_z_index(s.primary.distances, vec<6>(0.0), vec<6>(1.0));
}

std::string RA_to_string(float RA) {
	//get hour:
	float normalized = RA / (2.0 * PI);
	float hour = 24.0 * normalized;
	float minute = 60.0 * (hour - floor(hour));
	float second = 60.0 * (minute - floor(minute));


	return std::to_string(int(hour)) + "h " + std::to_string(int(minute)) + "m " + std::to_string(int(second)) + "s";
}
std::string DEC_to_string(float DEC) {
	float deg = 360 * DEC / (2.0 * PI);
	float minute = 60.0 * (deg - floor(deg));
	float second = 60.0 * (minute - floor(minute));
	return std::to_string(int(deg)) + " " + std::to_string(int(minute)) + "m " + std::to_string(int(second)) + "s";
}

void synthesize_database(unsigned int star_count, const char* path) {
	star* stars = generate_random_stars(star_count);
	binary_node<star, 3>* root = generate_binary_tree<star, 3>(stars, star_count, z_index_from_star);
	binary_node<star, 3>** leaf_nodes = collect_leaf_nodes(root, star_count);
	find_star_neighbors(leaf_nodes, star_count);
	store_database<star>(path, stars, star_count);
}

void test_orientation_from_centroids(vec<2>* centroids, int visible_stars, float fov, float identification_theshold, binary_node<star, 6>* root, vec<3>* forward, vec<3>* right, vec<3>* up, star** test_stars, int* matches, int* top_three_matches) {
	int variety = 3;
	//identify edge stars
	
	//extract possible star quads
	feature* star_quads = generate_star_quads_from_star_centroids(centroids, visible_stars);
	star*** star_candidates = new star * *[visible_stars];
	for (int i = 0; i < visible_stars; i++) {
		star** matches = find_matches(root, star_quads[i], visible_stars, variety);
		star_candidates[i] = matches;
	}
	float** scores = new float* [visible_stars];
	for (int i = 0; i < visible_stars; i++) {
		scores[i] = new float[variety];
		for (int j = 0; j < variety; j++) {
			scores[i][j] = 0.0;
		}
	}
	for (int i = 0; i < visible_stars; i++) {
		for (int j = i + 1; j < visible_stars; j++) {
			//find the ratio between the distance between the two stars and the distance between the first star and its nearest neghbor:
			float r = sqrt(dist_sq(centroids[i], centroids[j]));
			float observed_separation = sqrt(dist_sq(centroids[i], centroids[j])) * fov * 0.5;

			for (int a = 0; a < variety; a++) {
				for (int b = 0; b < variety; b++) {
					//check if it matches the values of the real data base
					float real_r = unit_vec_arc_length(star_candidates[i][a]->dir, star_candidates[j][b]->dir);
					float real_separation = unit_vec_arc_length(star_candidates[i][a]->dir, star_candidates[j][b]->dir);
					//float diff = real_r - r;
					float diff = real_separation - observed_separation;
					float factor = identification_theshold / (identification_theshold + 50.0 * diff * diff);
					scores[i][a] += factor;
					scores[j][b] += factor;
					/*
					if (abs(diff) < identification_theshold) {
						scores[i][a] += 1.0;
						scores[j][b] += 1.0;
					}
					*/
					
				}
			}
		}
	}
	//star_candidates, 17
	float highest_scores[] = { 0.0, 0.0, 0.0 };
	//first index is the centroid, the second is the variant index
	int best_stars[3][2] = { {0, 0}, {0, 0}, {0, 0} };
	for (int i = 0; i < visible_stars; i++) {
		//find the variation with the highest score:
		float high_score = 0.0;
		int local_best_star = 0;
		for (int j = 0; j < variety; j++) {
			if (scores[i][j] > high_score) {
				high_score = scores[i][j];
				local_best_star = j;
			}
		}
		if (high_score > highest_scores[2]) {
			highest_scores[2] = high_score;
			best_stars[2][0] = i;
			best_stars[2][1] = local_best_star;
			if (high_score > highest_scores[1]) {
				highest_scores[2] = highest_scores[1];
				highest_scores[1] = high_score;
				best_stars[2][0] = best_stars[1][0];
				best_stars[2][1] = best_stars[1][1];
				best_stars[1][0] = i;
				best_stars[1][1] = local_best_star;
				if (high_score > highest_scores[0]) {
					highest_scores[1] = highest_scores[0];
					highest_scores[0] = high_score;
					best_stars[1][0] = best_stars[0][0];
					best_stars[1][1] = best_stars[0][1];
					best_stars[0][0] = i;
					best_stars[0][1] = local_best_star;
				}
			}
		}
		//if (test_stars[i] == star_candidates[i][local_best_star]) {
		if (compare_stars_simple(*test_stars[i], *star_candidates[i][local_best_star])){
			*matches += 1;
		}
	}

	//check top 3:
	for (int i = 0; i < 3; i++) {
		if (compare_stars_simple(*test_stars[best_stars[i][0]], *star_candidates[best_stars[i][0]][best_stars[i][1]])) {
			*top_three_matches += 1;
		}
	}

	float epsilon = sin(fov / 2);
	vec<3> alpha = star_candidates[best_stars[0][0]][best_stars[0][1]]->dir;
	vec<3> beta = star_candidates[best_stars[1][0]][best_stars[1][1]]->dir;
	vec<3> gamma = star_candidates[best_stars[2][0]][best_stars[2][1]]->dir;
	matrix<3, 3> coefficents({ {
		{alpha[0], alpha[1], alpha[2]},
		{beta[0], beta[1], beta[2]},
		{gamma[0], gamma[1], gamma[2]}
	} });
	vec<3> a1({ epsilon * centroids[best_stars[0][0]][0], epsilon * centroids[best_stars[1][0]][0],epsilon * centroids[best_stars[2][0]][0] });
	vec<3> a2({ epsilon * centroids[best_stars[0][0]][1], epsilon * centroids[best_stars[1][0]][1],epsilon * centroids[best_stars[2][0]][1] });

	vec<3> x = solve_system_of_equations(coefficents, a1);
	vec<3> y = solve_system_of_equations(coefficents, a2);
	*right = x;
	*up = y;
	*forward = cross(y, x);
	return;
}
void orientation_from_centroids(vec<2>* centroids, int visible_stars, float fov, binary_node<star, 6>* root, vec<3>* forward, vec<3>* right, vec<3>* up) {
	int variety = 3;
	float threshold = 0.1;
	//identify edge stars

	//extract possible star quads
	feature* star_quads = generate_star_quads_from_star_centroids(centroids, visible_stars);
	star*** star_candidates = new star * *[visible_stars];
	for (int i = 0; i < visible_stars; i++) {
		star** matches = find_matches(root, star_quads[i], visible_stars, variety);
		star_candidates[i] = matches;
	}
	float** scores = new float* [visible_stars];
	for (int i = 0; i < visible_stars; i++) {
		scores[i] = new float[variety];
		for (int j = 0; j < variety; j++) {
			scores[i][j] = 0.0;
		}
	}
	for (int i = 0; i < visible_stars; i++) {
		for (int j = 0; j < visible_stars; j++) {
			if (i != j) {
				//find the ratio between the distance between the two stars and the distance between the first star and its nearest neghbor:
				float r = sqrt(dist_sq(centroids[i], centroids[j]));

				for (int a = 0; a < variety; a++) {
					for (int b = 0; b < variety; b++) {
						//check if it matches the values of the real data base
						float real_r = unit_vec_arc_length(star_candidates[i][a]->dir, star_candidates[j][b]->dir);


						if (abs(real_r - r) < threshold) {
							scores[i][a] += 1.0;
						}
					}
				}
			}
		}
	}

	float highest_scores[] = { 0.0, 0.0, 0.0 };
	//first index is the centroid, the second is variation
	int best_stars[3][2] = { {0, 0}, {0, 0}, {0, 0} };
	for (int i = 0; i < visible_stars; i++) {
		//find the variation with the highest score:
		float high_score = 0.0;
		int local_best_star = 0;
		for (int j = 0; j < variety; j++) {
			if (scores[i][j] > high_score) {
				high_score = scores[i][j];
				local_best_star = j;
			}
			if (scores[i][j] > highest_scores[2]) {
				highest_scores[2] = scores[i][j];
				best_stars[2][0] = i;
				best_stars[2][1] = j;
				if (scores[i][j] > highest_scores[1]) {
					highest_scores[2] = highest_scores[1];
					highest_scores[1] = scores[i][j];
					best_stars[2][0] = best_stars[1][0];
					best_stars[2][1] = best_stars[1][1];
					best_stars[1][0] = i;
					best_stars[1][1] = j;
					if (scores[i][j] > highest_scores[0]) {
						highest_scores[1] = highest_scores[0];
						highest_scores[0] = scores[i][j];
						best_stars[1][0] = best_stars[0][0];
						best_stars[1][1] = best_stars[0][1];
						best_stars[0][0] = i;
						best_stars[0][1] = j;
					}
				}
			}
		}
	}

	float epsilon = sin(fov / 2);
	vec<3> alpha = star_candidates[best_stars[0][0]][best_stars[0][1]]->dir;
	vec<3> beta = star_candidates[best_stars[1][0]][best_stars[1][1]]->dir;
	vec<3> gamma = star_candidates[best_stars[2][0]][best_stars[2][1]]->dir;
	matrix<3, 3> coefficents({ {
		{alpha[0], alpha[1], alpha[2]},
		{beta[0], beta[1], beta[2]},
		{gamma[0], gamma[1], gamma[2]}
	} });
	vec<3> a1({ epsilon * centroids[best_stars[0][0]][0], epsilon * centroids[best_stars[1][0]][0],epsilon * centroids[best_stars[2][0]][0] });
	vec<3> a2({ epsilon * centroids[best_stars[0][0]][1], epsilon * centroids[best_stars[1][0]][1],epsilon * centroids[best_stars[2][0]][1] });

	vec<3> x = solve_system_of_equations(coefficents, a1);
	vec<3> y = solve_system_of_equations(coefficents, a2);
	*right = x;
	*up = y;
	*forward = cross(y, x);
	return;
}

sky_atlas::sky_atlas(database<star> data, unsigned int subdivisions) {
	//create array to store maps
	int tiles_per_face = (subdivisions + 1) * (subdivisions + 1);
	maps = new sky_map[6 * tiles_per_face];
	tile_count = 6 * tiles_per_face;
	//std::cout << "var tiles = [";

	//iterate over the maps and set the star counts to 0 and set their directions:
	for (int i = 0; i < 6 * tiles_per_face; i++) {
		maps[i].star_count = 0;

		//determine u v cordinates, face, axis and direction
		int face = i / tiles_per_face;
		int u = (i - face * tiles_per_face) / (subdivisions + 1);
		int v = i - face * tiles_per_face - u * (subdivisions + 1);
		int axis = face % 3;
		float direction = (face < 3) * 2.0 - 1.0;

		//the direction of the map is the average of the directions of the 4 corners
		vec<3> corner_directions[4];
		for (int j = 0; j < 4; j++) {
			corner_directions[j][axis] = direction;
			corner_directions[j][(axis + 1) % 3] = 2.0 * (float(u + (j + (j / 2) % 2) % 2) / (subdivisions + 1)) - 1.0;
			corner_directions[j][(axis + 2) % 3] = 2.0 * (float(v + (j / 2) % 2) / (subdivisions + 1)) - 1.0;
		}
		vec<3> center_dir;
		center_dir[axis] = direction;
		center_dir[(axis + 1) % 3] = 2.0 * (float(u + 0.5) / (subdivisions + 1)) - 1.0;
		center_dir[(axis + 2) % 3] = 2.0 * (float(v + 0.5) / (subdivisions + 1)) - 1.0;
		center_dir = normalize(center_dir);
		//std::cout << "{vertices: [";
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 3; k++) {
				//std::cout << corner_directions[j][k] << ",";
			}
		}
		//std::cout << "],";
		//normalize corner vectors:
		for (int j = 0; j < 4; j++) {
			corner_directions[j] = normalize(corner_directions[j]);
		}
		float min_dot_product = 1.0;
		for (int j = 0; j < 4; j++) {
			for (int k = j; k < 4; k++) {
				if (j != k) {
					float dot_product = corner_directions[j] * corner_directions[k];
					if (dot_product < min_dot_product) {
						min_dot_product = dot_product;
					}
				}
			}
		}

		vec<3> vector_sum = { { 0, 0, 0 } };
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 4; k++) {
				vector_sum[j] += corner_directions[k][j];
			}
			vector_sum[j] /= 4.0;
		}
		maps[i].dir = center_dir;//vector_sum;
		maps[i].bounding_angle = min_dot_product;
		/*
		std::cout << "norm_vertices: [";
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 3; k++) {
				std::cout << corner_directions[j][k] << ",";
			}
		}
		std::cout << "], Dir: [";
		for (int k = 0; k < 3; k++) {
			std::cout << vector_sum[k] << ",";
		}
		std::cout << "], angle:" << maps[i].bounding_angle << "},";
		*/
	}
	//std::cout << "];\nexport {tiles};";

	//Iterate over the stars and assign each one to a tile in the star atlas
	unsigned int* indices = new unsigned int[data.object_count];
	for (int i = 0; i < data.object_count; i++) {
		//find dot profucts between the direction vector and each axis:
		float axis_similarities[] = {
			data.objects[i].dir * vec<3>({1.0, 0.0, 0.0}),
			data.objects[i].dir * vec<3>({0.0, 1.0, 0.0}),
			data.objects[i].dir * vec<3>({0.0, 0.0, 1.0})
		};

		//determine the closest axis;
		unsigned int closest_axis = 0;
		float most_similar = 0.0;
		for (int j = 0; j < 3; j++) {
			if (abs(axis_similarities[j]) > abs(most_similar)) {
				most_similar = axis_similarities[j];
				closest_axis = j;
			}
		}
		closest_axis += 3 * (axis_similarities[closest_axis] < 0.0);
		//project direction vector onto cube face
		unsigned int u = (subdivisions + 1) * (((axis_similarities[(closest_axis + 1) % 3] / abs(axis_similarities[closest_axis % 3])) + 1.0) / 2.0);//(subdivisions + 1) * (1.0 + axis_similarities[(closest_axis + 1) % 3]) / 2.0;
		unsigned int v = (subdivisions + 1) * (((axis_similarities[(closest_axis + 2) % 3] / abs(axis_similarities[closest_axis % 3])) + 1.0) / 2.0);//(subdivisions + 1) * (1.0 + axis_similarities[(closest_axis + 2) % 3]) / 2.0;
		unsigned int face = closest_axis;// +3 * (axis_similarities[closest_axis] < 0.0);
		unsigned int index = face * tiles_per_face + u * (subdivisions + 1) + v;

		indices[i] = index;
		maps[index].star_count++;
	}

	for (int i = 0; i < 6 * tiles_per_face; i++) {
		maps[i].stars = new star[maps[i].star_count];
		maps[i].star_count = 0;
	}
	for (int i = 0; i < data.object_count; i++) {
		maps[indices[i]].stars[maps[indices[i]].star_count++] = data.objects[i];
	}
	delete[] indices;
}
star* sky_atlas::get_visible_stars(vec<3> dir, float angle, int* star_count) {
	float dot_threshold = cos(angle);
	bool* maps_needed = new bool[tile_count];
	unsigned int stars_visible = 0;
	for (int i = 0; i < tile_count; i++) {
		if (maps[i].dir * dir > dot_threshold) {
			maps_needed[i] = true;
			stars_visible += maps[i].star_count;
		}
		else {
			maps_needed[i] = false;
		}
	}
	*star_count = stars_visible;
	star* stars = new star[stars_visible];
	unsigned int star_index = 0;
	for (int i = 0; i < tile_count; i++) {
		if (maps_needed[i]) {
			for (int j = 0; j < maps[i].star_count; j++) {
				stars[star_index++] = maps[i].stars[j];
			}
		}
	}
	return stars;
}
void sky_atlas::get_orientation(vec<2>* centroids, int visible_stars, float fov, float safety_fov, float identification_theshold, vec<3> estimated_dir, vec<3>* forward, vec<3>* right, vec<3>* up, star** true_stars, int* matches, int *top_three_matches, int* search_size) {
	//get possible stars
	star* possible_stars = get_visible_stars(estimated_dir, safety_fov, search_size);
	binary_node<star, 6>* quad_root = generate_binary_tree<star, 6>(possible_stars, *search_size, z_index_from_star_quad);
	test_orientation_from_centroids(centroids, visible_stars, fov, identification_theshold, quad_root, forward, right, up, true_stars, matches, top_three_matches);
	//orientation_from_centroids(centroids, visible_stars, fov, quad_root, forward, right, up);
	delete[] possible_stars;
	delete_branch(quad_root);
	return;
}