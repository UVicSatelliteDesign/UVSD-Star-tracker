#pragma once
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "sensor_sim.h"
#include <random>
#include <iostream>
#include <algorithm>
//#include "stb_image_write.h"


unsigned short** simulate_image(int width, int height, int adc_bits) {
	//intialize flux array
	float** flux = new float* [width];
	for (int i = 0; i < width; i++) {
		flux[i] = new float[height];
		for (int j = 0; j < height; j++) {
			flux[i][j] = 0.0;
		}
	}

	float threshold_flux_density = 0.001;

	int time_steps = 25;
	fvec2** centroids = new fvec2*[time_steps];
	for (int i = 0; i < time_steps; i++) {
		//project stars onto image plane

		//only keep ones 
	}

	float sensor_size = 0.2;


	std::mt19937 generator(0);

	float offscreen_factor = 1.4;
	std::uniform_real_distribution<float> uniform(-offscreen_factor * sensor_size /2, offscreen_factor * sensor_size / 2);
	std::uniform_real_distribution<float> brightness_gen(0.1, 1.0);


	float aspect_ratio = float(width) / height;
	for (int k = 0; k < 1200; k++) {
		float c_x0 = uniform(generator);
		float c_y0 = uniform(generator);
		float bounding_dist = 0.005;
		int steps = 45;

		float magnitude = brightness_gen(generator);
		magnitude *= magnitude / steps;
		if (k < 40) {
			magnitude *= 4.0;
		}


		float step_size = sensor_size / (18 * steps);
		for (int t = 0; t < steps; t++) {

			float center_x = c_x0 + step_size * t;
			float center_y = c_y0 + 0.5 * step_size* t;

			float max_x = center_x + bounding_dist;
			float min_x = center_x - bounding_dist;
			float max_y = center_y + bounding_dist;
			float min_y = center_y - bounding_dist;

			max_x = std::clamp(max_x, -aspect_ratio * sensor_size / 2, aspect_ratio * sensor_size / 2);
			min_x = std::clamp(min_x, -aspect_ratio * sensor_size / 2, aspect_ratio * sensor_size / 2);
			max_y = std::clamp(max_y, -sensor_size / 2, sensor_size / 2);
			min_y = std::clamp(min_y, -sensor_size / 2, sensor_size / 2);


			int top = int(height * (max_y / (sensor_size)+0.5));
			int bottom = int(height * (min_y / (sensor_size)+0.5));
			int left = int(width * (min_x / (sensor_size * aspect_ratio) + 0.5));
			int right = int(width * (max_x / (sensor_size * aspect_ratio) + 0.5));

			float x_step = (max_x - min_x) / (right - left);
			float y_step = (max_y - min_y) / (top - bottom);

			//top = std::clamp(top, 0, height);
			//bottom = std::clamp(bottom, 0, height);
			//left = std::clamp(left, 0, width);
			//right = std::clamp(right, 0, width);

			float y = min_y;
			for (int i = bottom; i < top; i++) {
				float x = min_x;
				for (int j = left; j < right; j++) {
					float vecs_dist_sq = (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y);
					flux[i][j] += magnitude * 2500.0 * exp(-0.5 * vecs_dist_sq / (0.0003 * 0.0003));
					x += x_step;
				}
				y += y_step;
			}
		}
	}

	int max_val = 1 << adc_bits - 1;




	unsigned short** pixels = new unsigned short* [width];
	for (int i = 0; i < width; i++) {
		pixels[i] = new unsigned short[height];
		for (int j = 0; j < height; j++) {
			if (flux[i][j] > 0) {
				std::poisson_distribution<> poisson(flux[i][j]);
				pixels[i][j] = unsigned short(std::min(max_val, poisson(generator)));
			}
			else {
				pixels[i][j] = 0;
			}
		}
	}

	//don't forget that exposure time will affect the thermal noise by averaging it out.
	float thermal_noise_mean = 3.5;
	float thermal_noise_deviation = 0.3;
	std::lognormal_distribution<float> distribution(thermal_noise_mean, thermal_noise_deviation);


	unsigned char* image_data = new unsigned char[width * height];
	unsigned char* p = image_data;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			//std::cout << distribution(generator) / max_val << "\n";
			*(p++) = unsigned char(255 * std::min(float(pixels[i][j]) / max_val + distribution(generator) / 255.0f, 1.0f));
		}
	}	stbi_write_png("./test_img.png", width, height, 1, image_data, 0);

	return pixels;
}
float compute_apparent_star_radius(star_camera camera, float magnitude) {
	float e_threshold = 1;
	float f_0 = 2.518e-8;
	float flux = camera.transmittance * camera.aperture * camera.aperture * std::_Pi * f_0 * log10(-0.4 * magnitude);
	return 2 * camera.planar_spread_parameter * camera.planar_spread_parameter * log(camera.quantum_efficeincy * flux  / (e_threshold * 2 * std::_Pi * camera.planar_spread_parameter * camera.planar_spread_parameter * camera.photon_energy));
}
void compute_guassian_template(star_field_generator* gen) {
	//compute the angular size of one pixel
	float pixel_angle = gen->camera.fov / gen->camera.height;

	//convert the standard deviation of the point spread from angular to 
	float sigma = gen->camera.angular_spread_parameter / pixel_angle;

	//choose an image size the at least covers the desired region
	unsigned int image_size = ceil(4 * sigma);

	//make sure the image size is an odd number to ensure the the brightest part of the distribution is preserved
	image_size += (1 - image_size % 2);

	//generate the distribution
	gen->guassian_template = generate_gaussian_image(image_size, image_size * pixel_angle, gen->camera.angular_spread_parameter);
}
bitmap<float> generate_gaussian_image(unsigned int image_size, float scale, float standard_deviation) {
	bitmap<float> image = init_bitmap(image_size, image_size, 1, 0.0f);
	float normalization_factor = 1.0 / (2.0 * std::_Pi * standard_deviation * standard_deviation);
	float exponent_coefficient = -0.5 / (standard_deviation * standard_deviation);
	float step = scale / image_size;
	float x = -scale * 0.5 + 0.5 * step;
	for (int i = 0; i < image_size; i++) {
		float y = -scale * 0.5 + 0.5 * step;
		for (int j = 0; j < image_size; j++) {
			float intensity = normalization_factor * exp(exponent_coefficient * (x * x + y * y));
			image_set_at(&image, i, j, intensity);
			y += step;
		}
		x += step;
	}	
	return image;
}
fvec2* visible_centroids_from_star_field_generator(star_field_generator gen, mat<float, 3, 3> orientation) {
	//this should be replaced with either a function which computed frustom geometry, or by some sort of property of the star tracker camera struct
	mat<float, 3, 3> pitch_up_matrix = rotation_mat(orientation[0], gen.camera.fov / 2);
	mat<float, 3, 3> pitch_down_matrix = rotation_mat(orientation[0], -gen.camera.fov / 2);
	mat<float, 3, 3> yaw_left_matrix = rotation_mat(orientation[1], gen.camera.fov / 2);
	mat<float, 3, 3> yaw_right_matrix = rotation_mat(orientation[1], -gen.camera.fov / 2);


	//determine the 4 cutting plane normals:
	fvec3 top_plane = normalize(mat_mult_vec(pitch_up_matrix, orientation[1]));
	fvec3 bottom_plane = normalize(mat_mult_vec(pitch_down_matrix, scale_vec(orientation[1], -1.0f)));
	fvec3 left_plane = normalize(mat_mult_vec(yaw_left_matrix, scale_vec(orientation[0], -1.0f)));
	fvec3 right_plane = normalize(mat_mult_vec(yaw_right_matrix, orientation[0]));


	//create array of output vectors:
	fvec2* centroids = new fvec2[gen.star_count];

	
	for (unsigned int i = 0; i < gen.star_count; i++) {
		//project the star's direction onto the orientation basis
		fvec3 position = transpose(orientation) * gen.stars[i].direction;


		if (position[2] > 0.0 && abs(position[0]) < gen.camera.aspect_ratio && abs(position[1]) < 1.0) {
			centroids[i][0] = position[0] / position[2];
			centroids[i][1] = position[1] / position[2];
		}
		else {
			centroids[i][0] = NAN;
			centroids[i][1] = NAN;
		}

	}

	return centroids;
}
star_field generate_star_field(star_field_generator gen, float start_time, float end_time) {
	//number of standard deviations between each neighboring point in the exposure
	const float gaussian_spacing = 0.25;
	
	//compute first and last orientation matrices
	mat<float, 3, 3> start_matrix = gen.path(start_time);
	mat<float, 3, 3> end_matrix = gen.path(end_time);

	//apply these matrices to the four corners of the frame.
	//these points will have the greatest apparant motion
	float y_offset = tan(gen.camera.fov / 2);
	float x_offset = gen.camera.aspect_ratio * y_offset;
	fvec3 top_left = init_vec({ 1.0f, -x_offset, y_offset });
	fvec3 top_right = init_vec({ 1.0f, x_offset, y_offset });
	fvec3 bottom_left = init_vec({ 1.0f, -x_offset, -y_offset });
	fvec3 bottom_right = init_vec({ 1.0f, x_offset, -y_offset });
	fvec3 start_top_left = start_matrix * top_left;
	fvec3 start_top_right = start_matrix * top_right;
	fvec3 start_bottom_left = start_matrix * bottom_left;
	fvec3 start_bottom_right = start_matrix * bottom_right;
	fvec3 end_top_left = end_matrix * top_left;
	fvec3 end_top_right = end_matrix * top_right;
	fvec3 end_bottom_left = end_matrix * bottom_left;
	fvec3 end_bottom_right = end_matrix * bottom_right;

	//Determine which points moved the most as an estimate on how far stars can move during the frame
	float rotation_angles[] = {
		vec_arc_length(start_top_left, end_top_left),
		vec_arc_length(start_top_right, end_top_right),
		vec_arc_length(start_bottom_left, end_bottom_left),
		vec_arc_length(start_bottom_right, end_bottom_right)
	};
	float largest_angle = 0.0f;
	for (int i = 0; i < 4; i++) {
		if (rotation_angles[i] > largest_angle) {
			largest_angle = rotation_angles[i];
		}
	}

	//determine the time step necesarry to acheive the desired smoothness
	int steps = ceil(largest_angle / (gaussian_spacing * gen.camera.angular_spread_parameter));
	float time_step = (end_time - start_time) / (steps - 1);
	
	//Initialize flux bitmap
	bitmap flux_bitmap = init_bitmap(gen.camera.width, gen.camera.height, 1, 0.0f);

	star_field field;
	mat<float, 3, 3> average_orientation = gen.path(0.5f * (start_time + end_time));
	float time = start_time;



	for (int i = 0; i < steps; i++) {
		mat<float, 3, 3> orientation = gen.path(time);

		fvec2* centroids = visible_centroids_from_star_field_generator(gen, orientation);

		//think: if a star goes off the screen during the exposure, should the centroid correspond to it's center of mass within the image, or it's true center mass including when it was off screen
		//if the latter, you need to first determin which stars are visible at all. You could argue that due to the nature of the gaussian distribution, all stars will be visible in all frames.
		//but this is unrealistic, so you need to determin a range / perimeter from the edge of the screen that depend on star brightness.
		for (int j = 0; j < gen.star_count; j++) {
			//if (!isnan(centroids[i][0])) {
				print_vector(centroids[i]);
				printf("\n");
			//}
		}


		time += time_step;
	}

	//Sample the flux bitmap to get a an image output;
	bitmap brightness_bitmap = init_bitmap(gen.camera.width, gen.camera.height, 1, (short)0);

	field.star_count = 0;
	field.bitmap = brightness_bitmap;
	
	
	return field;
}
void write_bitmap_to_png(const char* path, bitmap<unsigned char> img) {
	stbi_flip_vertically_on_write(true);
	stbi_write_png(path, img.width, img.height, img.channels, img.data, 0);
}