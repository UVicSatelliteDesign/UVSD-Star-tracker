#include "astrometry.h"
#include "database.h"
#include <random>
#include "linear_algebra.h"
#include "trees.h"
#include <string>
#define PI 3.14159265359


vec<3> point_on_sphere(float theta, float phi) {
	float x = cos(theta) * cos(phi);
	float y = sin(phi);
	float z = sin(theta) * cos(phi);


	return { x, y, z };
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
		std::string name = std::string("Star #") + std::to_string(i);
		int c = 0;
		while (c < name.length()) {
			stars[i].name[c] = name.at(c);
			c++;
		}
		stars[i].name[c] = '\0';
		stars[i].right_ascension = ra;
		stars[i].magnitude = 0.0;
		stars[i].declination = dec;
		stars[i].dir = point_on_sphere(ra, dec);
	}
	return stars;
}


int find_true_index(binary_node<star, 3>** stars, unsigned int star_count, binary_node<star, 3>* node) {
	int i = 0;
	while (node != stars[i] && i < star_count) {
		i++;
	}
	return i;
}
void find_star_neighbors(binary_node<star, 3>** stars, unsigned int star_count, star* start) {
	star_quad* stellar_web = new star_quad[star_count];
	for (int i = 0; i < star_count; i++) {
		binary_node<star, 3>** close_stars = find_k_nearest_neighbors<star, 3>(stars[i], stars[i]->components, 3, false);
		vec<3> d = stars[i]->key.object->dir;
		int c1 = find_true_index(stars, star_count, close_stars[0]);
		int c2 = find_true_index(stars, star_count, close_stars[1]);
		int c3 = find_true_index(stars, star_count, close_stars[2]);
		std::cout << "{dir: [" << d.components[0] << "," << d.components[1] << "," << d.components[2] << "], close: [" << c1 << "," << c2 << "," << c3 << "]},";
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
			stars[i]->key.object->primary.distances.components[c] = distances[c] / longest;
		}
		stars[i]->key.object->primary.longest_arc = longest;
		delete[] close_stars;
	}
}
star_quad* find_star_neighbors_redundancy(binary_node<star, 3>** stars, unsigned int star_count) {
	star_quad* stellar_web = new star_quad[star_count * 4];
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
				stellar_web[s].distances.components[c] = distances[c] / longest;
				stellar_web[s].longest_arc = longest;
			}
			s++;
		}
		delete[] close_stars;
	}
	return stellar_web;
}
unsigned long long int z_index_from_star(star s) {
	return  vec_to_z_index(s.dir, { -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 });
}



star** find_matches(binary_node<star, 6>* root, star_quad reference, unsigned int quad_count, int match_count) {


	int_vec<6> int_components = int_components_from_vec<6>(reference.distances, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });


	//construct star quad tree:
	//print_tree(root, print_star_quad);


	//find the cell which the reference quad falls into
	binary_node<star, 6>* host_cell = find_cell(root, reference.distances, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });


	//the star quad corresponding to the host cell is not nessessarily the most similar so nearest neighbor search is used:
	binary_node<star, 6>** nearby_quad_nodes = find_k_nearest_neighbors<star, 6>(host_cell, int_components, match_count - 1, false);

	star** possible_stars = new star*[match_count];
	for (int i = 0; i < match_count - 1; i++) {
		possible_stars[i] = nearby_quad_nodes[i]->key.object;
	}
	possible_stars[match_count - 1] = host_cell->key.object;

	return possible_stars;
}

star_quad* generate_star_quads_from_star_centroids(vec<2>* centroids, unsigned int star_count) {
	star_quad* star_quads = new star_quad[star_count];
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
			star_quads[i].distances.components[c] = distances[c] / longest;
		}
		star_quads[i].longest_arc = longest;
	}
	return star_quads;
}

vec<2>* generate_synthetic_image_data(vec<3> normal, float ccw_rotation, float fov, star* stars, unsigned int star_count, unsigned int* visible_star_count, star*** visible_stars_out, float noise) {
	vec<3> reference_up = { 0.0, 1.0, 0.0 };


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


	//Keep track of how many visible stars are in the image
	*visible_star_count = 0;


	//increment over the aray of stars and pick out the ones which are within the field of view
	for (unsigned int i = 0; i < star_count; i++) {
		vec<3> star_norm = stars[i].dir;


		float top_check = dot(star_norm, top_plane);
		float bottom_check = dot(star_norm, bottom_plane);
		float left_check = dot(star_norm, left_plane);
		float right_check = dot(star_norm, right_plane);


		if (top_check < 0.0 && bottom_check < 0.0 && left_check < 0.0 && right_check < 0.0) {
			visible_stars[*visible_star_count] = stars + i;
			(*visible_star_count)++;
		}
	}
	*visible_stars_out = visible_stars;
	//create array of output vectors:
	vec<2>* star_centroids = new vec<2>[*visible_star_count];

	float normalization_factor = sin(fov/2);
	int seed = 0;
	//project stars onto sensor space:
	for (unsigned int i = 0; i < *visible_star_count; i++) {
		//project star direction onto normal axis and get distance:
		float z = dot(visible_stars[i]->dir, normal);


		//get horizontal coordinate:
		float x = dot(right, visible_stars[i]->dir) / normalization_factor;


		//get vertical coordinate:
		float y = dot(up, visible_stars[i]->dir) / normalization_factor;

		star_centroids[i] = { x + random_float(seed++) * noise * 2, y + random_float(seed++) * noise * 2 };
	}


	return star_centroids;
}

unsigned long long int z_index_from_star_quad(star s) {
	return vec_to_z_index(s.primary.distances, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
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
void print_star(star s) {
	std::cout << s.name << ":\tR.A: " << RA_to_string(s.right_ascension) << "\tDEC: " << DEC_to_string(s.declination) << "\t x: " << s.dir.components[0] << "\ty: " << s.dir.components[1] << "\tz: " << s.dir.components[2] << "\tZ index: ";
}
void print_stars(star* stars, unsigned int star_count) {
	for (int i = 0; i < star_count; i++) {
		print_star(stars[i]);
		std::cout << "\n";
	}
}


void synthesize_database(unsigned int star_count, const char* path) {
	star* stars = generate_random_stars(star_count);
	binary_node<star, 3>* root = generate_binary_tree<star, 3>(stars, star_count, z_index_from_star);
	binary_node<star, 3>** leaf_nodes = collect_leaf_nodes(root, star_count);
	find_star_neighbors(leaf_nodes, star_count, stars);
	store_database<star>(path, stars, star_count);
}

void test_orientation_from_centroids(vec<2>* centroids, int visible_stars, float fov, binary_node<star, 6>* root, vec<3>* forward, vec<3>* right, vec<3>* up, star** test_stars, int* matches, int* top_three_matches) {
	int variety = 3;
	float threshold = 0.1;
	//identify edge stars

	//extract possible star quads
	star_quad* star_quads = generate_star_quads_from_star_centroids(centroids, visible_stars);
	star*** star_candidates = new star * *[visible_stars];
	for (int i = 0; i < visible_stars; i++) {
		star** matches = find_matches(root, star_quads[i], visible_stars, variety);
		star_candidates[i] = matches;
	}
	float** scores = new float*[visible_stars];
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
				float r = sqrt(dist_sq(centroids[i], centroids[j])) / star_quads[i].longest_arc;

				for (int a = 0; a < variety; a++) {
					for (int b = 0; b < variety; b++) {
						//check if it matches the values of the real data base
						float real_r = unit_vec_arc_length(star_candidates[i][a]->dir, star_candidates[j][b]->dir) / star_candidates[i][a]->primary.longest_arc;


						if (abs(real_r - r) < threshold) {
							scores[i][a] += 1.0;
						}
					}
				}
			}
		}
	}
	//star_candidates, 17
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
		if (test_stars[i] == star_candidates[i][local_best_star]) {
			*matches += 1;
		}
	}

	//check top 3:
	for (int i = 0; i < 3; i++) {
		if (test_stars[best_stars[i][0]] == star_candidates[best_stars[i][0]][best_stars[i][1]]) {
			*top_three_matches += 1;
		}
	}
	float epsilon = sin(fov / 2);
	vec<3> alpha = star_candidates[best_stars[0][0]][best_stars[0][1]]->dir;
	vec<3> beta = star_candidates[best_stars[1][0]][best_stars[1][1]]->dir;
	vec<3> gamma = star_candidates[best_stars[2][0]][best_stars[2][1]]->dir;
	matrix<3, 3> coefficents({ {
		{alpha.components[0], alpha.components[1], alpha.components[2]},
		{beta.components[0], beta.components[1], beta.components[2]},
		{gamma.components[0], gamma.components[1], gamma.components[2]}
	} });
	vec<3> a1 = { epsilon * centroids[best_stars[0][0]].components[0], epsilon * centroids[best_stars[1][0]].components[0],epsilon * centroids[best_stars[2][0]].components[0] };
	vec<3> a2 = { epsilon * centroids[best_stars[0][0]].components[1], epsilon * centroids[best_stars[1][0]].components[1],epsilon * centroids[best_stars[2][0]].components[1] };

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
	star_quad* star_quads = generate_star_quads_from_star_centroids(centroids, visible_stars);
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
				float r = sqrt(dist_sq(centroids[i], centroids[j])) / star_quads[i].longest_arc;

				for (int a = 0; a < variety; a++) {
					for (int b = 0; b < variety; b++) {
						//check if it matches the values of the real data base
						float real_r = unit_vec_arc_length(star_candidates[i][a]->dir, star_candidates[j][b]->dir) / star_candidates[i][a]->primary.longest_arc;


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
		{alpha.components[0], alpha.components[1], alpha.components[2]},
		{beta.components[0], beta.components[1], beta.components[2]},
		{gamma.components[0], gamma.components[1], gamma.components[2]}
	} });
	vec<3> a1 = { epsilon * centroids[best_stars[0][0]].components[0], epsilon * centroids[best_stars[1][0]].components[0],epsilon * centroids[best_stars[2][0]].components[0] };
	vec<3> a2 = { epsilon * centroids[best_stars[0][0]].components[1], epsilon * centroids[best_stars[1][0]].components[1],epsilon * centroids[best_stars[2][0]].components[1] };

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
			corner_directions[j].components[axis] = direction;
			corner_directions[j].components[(axis + 1) % 3] = 2.0 * (float(u+(j+ (j / 2) % 2)%2) / (subdivisions + 1)) - 1.0;
			corner_directions[j].components[(axis + 2) % 3] = 2.0 * (float(v+(j/2)%2) / (subdivisions + 1)) - 1.0;
		}
		//std::cout << "{vertices: [";
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 3; k++) {
				//std::cout << corner_directions[j].components[k] << ",";
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
					float dot_product = dot(corner_directions[j], corner_directions[k]);
					if (dot_product < min_dot_product) {
						min_dot_product = dot_product;
					}
				}
			}
		}

		vec<3> vector_sum = {{ 0, 0, 0 }};
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 4; k++) {
				vector_sum.components[j] += corner_directions[k].components[j];
			}
			vector_sum.components[j] /= 4.0;
		}
		maps[i].dir = vector_sum;
		maps[i].bounding_angle = min_dot_product;
		/*
		std::cout << "norm_vertices: [";
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 3; k++) {
				std::cout << corner_directions[j].components[k] << ",";
			}
		}
		std::cout << "], Dir: [";
		for (int k = 0; k < 3; k++) {
			std::cout << vector_sum.components[k] << ",";
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
			dot(data.objects[i].dir, {1.0, 0.0, 0.0}),
			dot(data.objects[i].dir, {0.0, 1.0, 0.0}),
			dot(data.objects[i].dir, {0.0, 0.0, 1.0})
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

		//project direction vector onto cube face
		unsigned int u = (subdivisions + 1) * (((axis_similarities[(closest_axis + 1) % 3] / axis_similarities[closest_axis]) + 1.0) / 2.0);//(subdivisions + 1) * (1.0 + axis_similarities[(closest_axis + 1) % 3]) / 2.0;
		unsigned int v = (subdivisions + 1) * (((axis_similarities[(closest_axis + 2) % 3] / axis_similarities[closest_axis]) + 1.0) / 2.0);//(subdivisions + 1) * (1.0 + axis_similarities[(closest_axis + 2) % 3]) / 2.0;
		unsigned int face = closest_axis + 3 * (most_similar < 0.0);
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
		if (dot(maps[i].dir, dir) > dot_threshold) {
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
void sky_atlas::get_orientation(vec<2>* centroids, int visible_stars, float fov, float safety_fov, vec<3> estimated_dir, vec<3>* forward, vec<3>* right, vec<3>* up, star** true_stars) {
	//get possible stars
	int star_count = 0;
	int matches = 0;
	int top_three_matches = 0;
	star* possible_stars = get_visible_stars(estimated_dir, safety_fov, &star_count);
	binary_node<star, 6>* quad_root = generate_binary_tree<star, 6>(possible_stars, star_count, z_index_from_star_quad);
	test_orientation_from_centroids(centroids, visible_stars, fov, quad_root, forward, right, up, true_stars, &matches, &top_three_matches);
	//orientation_from_centroids(centroids, visible_stars, fov, quad_root, forward, right, up);
	return;
}