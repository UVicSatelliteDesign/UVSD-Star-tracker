#include <iostream>
#include <fstream>
#include <intrin.h>
#include <random>
#include <string>
#include <math.h>


#include "database.h"
#include "astrometry.h"
#include "linear_algebra.h"
#include "trees.h"
#define PI 3.14159265359

void test_star_quad_identification_succes_rate_vs_noise(database<star> data, unsigned int star_limit, unsigned int samples, float start_noise, float end_noise, unsigned int steps, const char* folder) {
	if (star_limit > data.object_count) {
		star_limit = data.object_count;
	}
	//build tree:
	binary_node<star, 6>* quad_root = generate_binary_tree<star, 6>(data.objects, star_limit, z_index_from_star_quad);
	unsigned int seed = 0;
	FILE* destination;
	std::string path = std::string(folder) + std::to_string(star_limit) + "_stars_" + std::to_string(samples) + "_samples_" + std::to_string(start_noise) + "_to_" + std::to_string(end_noise) + "_noise.txt";
	fopen_s(&destination, path.c_str(), "w");
	if (!destination) {
		return;
	}
	float noise = start_noise;
	float noise_step = (end_noise - start_noise) / (steps - 1);
	for (int i = 0; i < steps; i++) {
		unsigned int successes = 0;
		for (int s = 0; s < samples; s++) {
			unsigned int star_index = (unsigned int)(random_float(seed++) * (star_limit - 1));
			star_quad reference = data.objects[star_index].primary;
			for (int c = 0; c < 6; c++) {
				reference.distances.components[c] += noise * random_float(seed++);
			}
			star** guess = find_matches(quad_root, reference, star_limit, 1);
			if (*guess == data.objects + star_index) {
				successes++;
			}
		}
		//fprintf(destination, "%f,%f\n", noise, float(successes)/ float(samples));
		fprintf(destination, "%f\n", float(successes)/ float(samples));
		noise += noise_step;
	}
	fclose(destination);
}

//k is the degree if the edge stars, k = 0 corresponds to checking for stars where the edge is closer than the closest neighbor, k = 2 correponds the edge being closer than the thrid neighbor
template <unsigned int k>
void test_edge_star_proportion_vs_stars_fov(database<star> data, unsigned int samples, float start_fov, float end_fov, unsigned int steps, const char* folder) {
	FILE* destination;
	std::string path = std::string(folder) + std::to_string(data.object_count) + "_stars_" + std::to_string(samples) + "_samples_edge_stars_vs_visible_stars.txt";
	fopen_s(&destination, path.c_str(), "w");
	if (!destination) {
		return;
	}
	float fov = start_fov;
	float fov_step = (end_fov - start_fov) / (steps - 1);
	unsigned int seed = 0;
	for (int i = 0; i < steps; i++) {
		unsigned int total_stars = 0;
		unsigned int total_edge_stars = 0;
		for (int s = 0; s < samples; s++) {
			star** visible_stars;
			unsigned int visible_star_count = 0;
			vec<2>* centroids = generate_synthetic_image_data(random_float(seed++) * 6.28, asin(random_float(seed++) * 2.0 - 1.0), random_float(seed++) * 6.28, fov, data.objects, data.object_count, &visible_star_count, &visible_stars, 0.0);
			if (visible_star_count > 3) {
				star_quad* image_quads = generate_star_quads_from_star_centroids(centroids, visible_star_count);
				for (int j = 0; j < visible_star_count; j++) {
					//check distance to the edge
					float dist = std::min(std::min(1.0f - centroids[j].components[0], 1.0f + centroids[j].components[0]), std::min(1.0f - centroids[j].components[1], 1.0f + centroids[j].components[1]));
					float target_star = image_quads[j].distances.components[k] * image_quads[j].longest_arc;
					if (dist < target_star) {
						total_edge_stars++;
					}
				}
				delete[] image_quads;
			}
			else {
				total_edge_stars += visible_star_count;
			}
			total_stars += visible_star_count;
			delete[] visible_stars;
			delete[] centroids;
		}
		fprintf(destination, "%f,%f\n", float(total_stars) / samples, float(total_edge_stars) / float(samples));
		fov += fov_step;
	}
	fclose(destination);
}


void test_average_distance(database<star> data, unsigned int samples, const char* folder) {
	if (samples > data.object_count) {
		samples = data.object_count;
	}
	FILE* destination;
	std::string path = std::string(folder) + std::to_string(samples) + "_star_distance_histogram.txt";
	fopen_s(&destination, path.c_str(), "w");
	if (!destination) {
		return;
	}

	float* first_distances = new float[samples];
	float* second_distances = new float[samples];
	float* third_distances = new float[samples];
	float distance_sum[3] = { 0.0, 0.0, 0.0};
	for (int i = 0; i < samples; i++) {
		float distance1 = data.objects[i].primary.distances.components[0] * data.objects[i].primary.longest_arc;
		float distance2 = data.objects[i].primary.distances.components[1] * data.objects[i].primary.longest_arc;
		float distance3 = data.objects[i].primary.distances.components[2] * data.objects[i].primary.longest_arc;
		distance_sum[0] += distance1;
		distance_sum[1] += distance2;
		distance_sum[2] += distance3;
		if (distance1 != distance1) {
			printf("here is da problem\n");
		}
		first_distances[i] = distance1;
		second_distances[i] = distance2;
		third_distances[i] = distance3;
		//fprintf_s(destination, "%f,%f,%f\n", distance1, distance2, distance3);
	}
	float min = 0.0;
	float max = 0.0;
	for (int i = 0; i < samples; i++) {
		if (first_distances[i] > max) {
			max = first_distances[i];
		}
		if (second_distances[i] > max) {
			max = second_distances[i];
		}
		if (third_distances[i] > max) {
			max = third_distances[i];
		}
	}
	float range = max - min;
	int bin_count = 100;
	float* first_bins = new float[bin_count];
	float* second_bins = new float[bin_count];
	float* third_bins = new float[bin_count];
	for (int i = 0; i < bin_count; i++) {
		first_bins[i] = 0.0;
		second_bins[i] = 0.0;
		third_bins[i] = 0.0;
	}
	for (int i = 0; i < samples; i++) {
		//get bin index:
		int index1 = bin_count * (first_distances[i] - min) / (range);
		int index2 = bin_count * (second_distances[i] - min) / (range);
		int index3 = bin_count * (third_distances[i] - min) / (range);
		first_bins[index1] += first_distances[i] / distance_sum[0];
		second_bins[index2] += second_distances[i] / distance_sum[1];
		third_bins[index3] += third_distances[i] / distance_sum[2];
	}
	fprintf_s(destination, "min:%f, max:%f\nfirst,second,third\n", min, max);
	for (int i = 0; i < bin_count; i++) {
		float bin_val = range * (float(i) / bin_count);
		fprintf_s(destination, "%f,%f,%f\n", first_bins[i], second_bins[i], third_bins[i]);
	}
	std::cout << distance_sum[0] / samples << ", " << distance_sum[1] / samples << ", " << distance_sum[2] / samples << "\n";
	fclose(destination);
}
void test_identification(database<star> data) {
	unsigned int seed = 0;
	unsigned int visible_star_count = 0;
	float fov = 0.4;
	star** visible_stars;
	vec<2>* centroids = generate_synthetic_image_data(random_float(seed++) * 6.28, asin(random_float(seed++) * 2.0 - 1.0), random_float(seed++) * 6.28, fov, data.objects, data.object_count, &visible_star_count, &visible_stars, 0.0);
	binary_node<star, 6>* quad_root = generate_binary_tree<star, 6>(data.objects, data.object_count, z_index_from_star_quad);

	orientation_from_centroids(centroids, visible_star_count, quad_root, NULL, NULL, NULL, visible_stars);
}

int main() {
	srand(34);
	std::setprecision(4);
	std::cout << std::fixed;

	/*
	unsigned int star_count = 2000;
	star* stars = generate_random_stars(star_count);
	binary_node<star, 3>* root = generate_binary_tree<star, 3>(stars, star_count, z_index_from_star);//divide_cell<star, 3>(stars, 0, star_count - 1);
	//std::cout << "\n\nStar tree:\n";
	//print_tree(root, print_star);


	binary_node<star, 3>** leaf_nodes = collect_leaf_nodes<star, 3>(root, star_count);
	find_star_neighbors(leaf_nodes, star_count);
	binary_node<star, 6>* quad_root = generate_binary_tree<star, 6>(stars, star_count, z_index_from_star_quad);
	int* match_count = new int[500];
	float* successs = new float[500];
	for (int i = 0; i < 500; i++) {
		match_count[i] = 0;
		successs[i] = 0.0;
	}
	int seed = 521;
	for (int k = 50; k < 60; k++) {
		int total_stars_detected = 0;
		int total_successful_identifications = 0;
		int best_star_successes = 0;
		int samples = 10;
		int valid_samples = 0;///sample where there are enough visible stars to form a quad
		for (int s = 0; s < samples; s++) {
			star** test_stars;
			unsigned int visible_stars;
			vec<2>* centroids = generate_synthetic_image_data(random_float(seed++) * 6.28, asin(random_float(seed++) * 2.0 - 1.0), random_float(seed++) * 6.28, float(k) * 0.006 + 0.065, stars, star_count, &visible_stars, &test_stars, 0.01);
			if (visible_stars >= 4) {
				valid_samples++;
				star_quad* image_quads = generate_star_quads_from_star_centroids(centroids, visible_stars);


				star** extracted_quads = new star*[visible_stars];
				for (int i = 0; i < visible_stars; i++) {
					extracted_quads[i] = find_match(quad_root, image_quads[i], star_count, 0.1f);
				}


				int matches = 0;
				float* scores = new float[visible_stars];
				float threshold = 0.1;
				for (int i = 0; i < visible_stars; i++) {
					scores[i] = 0.0;
					for (int j = 0; j < visible_stars; j++) {
						if (i != j) {
							//find the ratio between the distance between the two stars and the distance between the first star and its nearest neghbor:
							float r = sqrt(dist_sq(centroids[i], centroids[j])) / image_quads[i].longest_arc;


							//check if it matches the values of the real data base
							float real_r = unit_vec_arc_length(extracted_quads[i]->dir, extracted_quads[j]->dir) / extracted_quads[i]->primary.longest_arc;


							if (abs(real_r - r) < threshold) {
								scores[i] += 1.0;
							}
						}
					}
					std::cout << "score " << i << ":\t" << scores[i] << ": ";
					if (test_stars[i] == extracted_quads[i]) {
						std::cout << " Match";
						matches++;
					}
					std::cout << "\n";
				}


				float high_score = 0;
				float second_high_score = 0;
				int best_star = 0;
				int next_best_star;
				for (int i = 0; i < visible_stars; i++) {
					if (scores[i] > high_score) {
						high_score = scores[i];
						best_star = i;
					}
					else if (scores[i] > second_high_score) {
						second_high_score = scores[i];
						next_best_star = i;
					}
				}






				match_count[visible_stars] += 1;
				successs[visible_stars] += float(matches) / visible_stars;


				total_successful_identifications += matches;
				std::cout << "Run #" << s << "\t Stars in image: " << visible_stars << "\t Matches: " << matches << "\t highest score: " << high_score;


				if (test_stars[best_star] == extracted_quads[best_star]) {
					std::cout << " ***";
					best_star_successes++;
				}
				std::cout << "\n";
				delete[] image_quads;
				delete[] extracted_quads;
			}
			total_stars_detected += visible_stars;
			delete[] test_stars;


			delete[] centroids;
		}
		std::cout << "Valid samples: " << valid_samples << "\tIndividial star identification rate: " << float(total_successful_identifications) / float(total_stars_detected) << "\taverage stars detected: " << float(total_stars_detected) / float(samples) << "\tBest star identified correctly: " << float(best_star_successes) / float(samples) << "\tBest star identification rate: " << float(best_star_successes) / float(valid_samples) << "\t average matches: " << float(total_successful_identifications) / samples << "\n";
	}
	*/
	/*
	for (int i = 0; i < 200; i++) {
		std::cout << successs[i] / match_count[i] << "\n";
	}
	*/
	const char* test_folder = "C:/Users/logac/Desktop/UVSD Star tracker/tests/";
	//synthesize_database(150000, "C:/Users/logac/Desktop/UVSD Star tracker/database_150000.star");
	database<star> data = load_database<star>("C:/Users/logac/Desktop/UVSD Star tracker/database_16000.star");
	
	test_identification(data);
	//test_average_distance(data, 150000, test_folder);
	//test_edge_star_proportion_vs_stars_fov(data, 1000, 0.1, 1.5, 250, "C:/Users/logac/Desktop/UVSD Star tracker/tests/");
	//export_synthetic_centroids(data, 0.5, 0.3, 0.0, 1.5, test_folder);
	return 0;
}