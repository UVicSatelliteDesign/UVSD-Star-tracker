#pragma once
#include <string>
#include "linear_algebra.h"
#include "trees.h"

struct star_quad {
	vec<6> distances;
	float longest_arc;
};
struct star {
	char name[32];
	float right_ascension;//between 0 and 2PI
	float declination;//between -pi and +pi
	float magnitude;//just the plain old astronomical magnitude -1.0 to ~7.0
	vec<3> dir;//the normal direction of the star from the perspective of the earth and the satellite 
	star_quad primary;
};
vec<3> point_on_sphere(float theta, float phi);
unsigned long long int z_index_from_star(star s);
star* generate_random_stars(unsigned int star_count);
void synthesize_database(unsigned int star_count, const char* path);
void find_star_neighbors(binary_node<star, 3>** stars, unsigned int star_count);
star** find_matches(binary_node<star, 6>* root, star_quad reference, unsigned int quad_count, int match_count);
star_quad* generate_star_quads_from_star_centroids(vec<2>* centroids, unsigned int star_count);
vec<2>* generate_synthetic_image_data(float right_ascension, float declineation, float ccw_rotation, float fov, star* stars, unsigned int star_count, unsigned int* visible_star_count, star*** visible_stars_out, float noise);
void orientation_from_centroids(vec<2>* centroids, binary_node<star, 6>* root, vec<3>* forward, vec<3>* right, vec<3>* up);
float random_float(int seed);
std::string RA_to_string(float RA);
std::string DEC_to_string(float DEC);
unsigned long long int z_index_from_star_quad(star s);
