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
	vec<2>** centroids = new vec<2>*[time_steps];
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
					float dist_sq = (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y);
					flux[i][j] += magnitude * 2500.0 * exp(-0.5 * dist_sq / (0.0003 * 0.0003));
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
float compute_apparent_star_radius(star_camera camera, star_field_star star) {
	return 0.0;
}
void compute_guassian_template(star_field_generator gen) {
	float brightest_magnitude = -1.0;
}
void write_image_to_png(const char* path, image<unsigned char> img) {
	stbi_flip_vertically_on_write(true);
	stbi_write_png(path, img.width, img.height, img.channels, img.data, 0);
}