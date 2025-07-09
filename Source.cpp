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
			feature reference = data.objects[star_index].primary;
			for (int c = 0; c < 6; c++) {
				reference.distances[c] += noise * random_float(seed++);
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
			vec<3> normal = point_on_sphere(random_float(seed++) * 6.28, asin(random_float(seed++) * 2.0 - 1.0));
			vec<2>* centroids = generate_synthetic_image_data(normal, random_float(seed++) * 6.28, fov, data.objects, data.object_count, &visible_star_count, &visible_stars, 0.0, seed++);
			if (visible_star_count > 3) {
				feature* image_quads = generate_star_quads_from_star_centroids(centroids, visible_star_count);
				for (int j = 0; j < visible_star_count; j++) {
					//check distance to the edge
					float dist = std::min(std::min(1.0f - centroids[j][0], 1.0f + centroids[j][0]), std::min(1.0f - centroids[j][1], 1.0f + centroids[j][1]));
					float target_star = image_quads[j].distances[k];
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
		float distance1 = data.objects[i].primary.distances[0];
		float distance2 = data.objects[i].primary.distances[1];
		float distance3 = data.objects[i].primary.distances[2];
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
void test_identification(database<star> data, int samples) {
	unsigned int seed = 0;
	unsigned int visible_star_count = 0;
	float fov = 0.12;

	std::cout << "Orientation determination test. Fov:" << fov << "\tSamples:" << samples << "\n";

	int total_visible_stars = 0;
	int total_matches = 0;
	int total_successes = 0;
	binary_node<star, 6>* quad_root = generate_binary_tree<star, 6>(data.objects, data.object_count, z_index_from_star_quad);
	for (int s = 0; s < samples; s++) {
		star** visible_stars;
		vec<3> normal = point_on_sphere(random_float(seed++) * 6.28, asin(random_float(seed++) * 2.0 - 1.0));

		vec<2>* centroids = generate_synthetic_image_data(normal, random_float(seed++) * 6.28, fov, data.objects, data.object_count, &visible_star_count, &visible_stars, 0.0, seed++);

		vec<3> right;
		vec<3> up;
		vec<3> forward;
		int matches = 0;
		int top_three_matches = 0;
		if (visible_star_count >= 3) {
			test_orientation_from_centroids(centroids, visible_star_count, fov, 0.01, quad_root, &forward, &right, &up, visible_stars, &matches, &top_three_matches);
		}
		if (top_three_matches == 3) {
			total_successes++;
		}

		total_matches += matches;
		total_visible_stars += visible_star_count;

		delete[] centroids;
		delete[] visible_stars;
	}
	std::cout << "Average visible stars: " << float(total_visible_stars) / samples << "\tAverage matches: " << float(total_matches) / samples << "\tSuccess rate: " << 100.0 * float(total_successes) / samples << "%\t";
}

void test_tiled_identification(database<star> data, unsigned int subdivisions, unsigned int samples, float fov, float proportional_position_noise, float identification_theshold) {
	unsigned int seed = 0;
	sky_atlas test_atlas(data, subdivisions);

	float tile_angle = acos(1.0 / sqrt(1.0 + 2.0 / ((subdivisions + 1) * (subdivisions + 1))));

	int total_visible_stars = 0;
	int total_matches = 0;
	int total_successes = 0;
	int total_stars_searched = 0;
	for (int s = 0; s < samples; s++) {

		

		star** visible_stars;
		unsigned int visible_star_count = 0;
		vec<3> normal = point_on_sphere(random_float(seed++) * 6.28, asin(random_float(seed++) * 2.0 - 1.0));
		vec<2>* centroids = generate_synthetic_image_data(normal, random_float(seed++) * 6.28, fov, data.objects, data.object_count, &visible_star_count, &visible_stars, proportional_position_noise, seed++);
		int matches = 0;
		vec<3> right;
		vec<3> up;
		vec<3> forward;
		int top_three_matches = 0;
		int search_size = 0;

		if (visible_star_count > 3) {
			test_atlas.get_orientation(centroids, visible_star_count, fov, fov + tile_angle, identification_theshold, normal, &forward, &right, &up, visible_stars, &matches, &top_three_matches, &search_size);
			vec<3> error = sub_vector(forward, normal);
			total_stars_searched += search_size;

			//print_vector(forward);
			//std::cout << s << "\n";

			//std::cout << "error: " << length(error) << "\n";
		}
		else {
			//std::cout << "Insufficient stars \n";
		}
		if (top_three_matches == 3) {
			total_successes++;
		}
		else {
			//std::cout <<"Sample " << s << " failed \n";
		}
		total_matches += matches;
		total_visible_stars += visible_star_count;


		delete[] centroids;
		delete[] visible_stars;
	}
	std::cout << float(total_visible_stars) / samples << ", " << proportional_position_noise << ", " << 100.0 * float(total_successes) / samples << ",\n";

	//std::cout << "Fov: " << fov << " Average visible stars: " << float(total_visible_stars) / samples << "\tAverage matches: " << float(total_matches) / samples << "\tSuccess rate: " << 100.0 * float(total_successes) / samples << "%\n";

}


int main() {
	srand(34);
	std::setprecision(4);
	std::cout << std::fixed;

	const char* test_folder = "C:/Users/logac/Desktop/UVSD Star tracker/tests/";
	synthesize_database(16000, "C:/Users/Logan/Desktop/UVSD-Star-tracker/database_16000.star");
	//database<star> data = load_database<star>("./database_16000.star");


	//test_tiled_identification(data, 10, 10, 0.13, 0.01, 0.0001);
	/*
	//test_identification(data, 6000);
	int grid_resolution = 40;
	for (int i = 0; i < grid_resolution; i++) {
		for (int j = 1; j < grid_resolution + 1; j++) {
			test_tiled_identification(data, 10, 1000, 0.17 * float(i) / (grid_resolution - 1), 0.02 * float(j) / (grid_resolution - 1), 0.0001);
		}
	}
	*/
	//test_average_distance(data, 150000, test_folder);
	//test_edge_star_proportion_vs_stars_fov(data, 1000, 0.1, 1.5, 250, "C:/Users/logac/Desktop/UVSD Star tracker/tests/");
	//export_synthetic_centroids(data, 0.5, 0.3, 0.0, 1.5, test_folder);
	return 0;
}