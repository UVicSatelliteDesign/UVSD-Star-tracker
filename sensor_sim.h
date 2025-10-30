#pragma once
#include "linear_algebra.h"
#include <algorithm>
#include "stb_image_write.h"
#include <random>

struct star_field_centroid {
	unsigned int id;
	fvec2 centroid;
};
struct star_field_star {
	fvec3 direction;
	float magnitude;
	float bounding_radius;//this is just used to speed up the drawing process so that the entire gaussian template isn't copied every time
};

template <typename T>
struct bitmap {
	unsigned int width;
	unsigned int height;
	unsigned int channels;
	T* data;
};
struct star_field {
	bitmap<short> bitmap;
	unsigned int star_count;
	star_field_centroid* centroids;//the centroids of all stars that were visible at any point during the frame
	mat<float, 3, 3> start_orientation;
	mat<float, 3, 3> end_orientation;
	mat<float, 3, 3> average_orientation;//time averaged orientation
};
struct star_camera {
	float aperture;//lens diameter. meters
	float fov;//vertical field of view in radians
	float aspect_ratio;//ratio of the width to the height of the image
	float focal_length;//meters
	float quantum_efficeincy;//portion of photons which are converted to electons
	float photon_energy;//average energy of photon detected by this sensor
	float angular_spread_parameter;
	float planar_spread_parameter;//one standard deviation of a star on the sensor. Measured in meters
	float transmittance;//the portion of light which passes through the optics and is received at the sensor. Ideally 100%
	float pixel_area;//The physical size of a single pixel on the sensor. Which may be smaller than the area of the sensor divided by the number of pixels, hence why it is a separate parameter. square meters


	float sensor_size;//height of the sensor in meters
	unsigned int width;//image width in pixels
	unsigned int height;//image height in pixels
	unsigned int full_well;//number of electrons required to saturate a pixel

	//this should probably be part of a separate structure describing the capture parameters, as it's not an attribute of the optical system itself.

};
struct star_field_generator {
	star_camera camera;
	unsigned int star_count;
	star_field_star* stars;
	bool positional_noise;
	float frame_interval;//amound of time between frames (realistically longer than the exposure time)
	float exposure_time;//amount of time collecting light
	std::default_random_engine generator;
	std::normal_distribution<float> distribution;
	mat<float, 3, 3>(*path)(float);
	bitmap<float> guassian_template;
};
star_camera default_star_camera();

template <typename T>
bitmap<T> init_bitmap(unsigned int width, unsigned height, unsigned int channels, T initial_value) {
	bitmap<T> img;
	img.width = width;
	img.height = height;
	img.channels = channels;
	unsigned int size = width * height * channels;
	img.data = new T[size];
	for (int i = 0; i < size; i++) {
		img.data[i] = initial_value;
	}
	return img;
}

template <typename T>
bitmap<T> init_checkerboard(unsigned int width, unsigned int height, unsigned int channels, T high_value, unsigned int period) {
	bitmap<T> img;
	img.width = width;
	img.height = height;
	img.channels = channels;
	unsigned int size = width * height * channels;
	img.data = new T[size];
	int i = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			for (int c = 0; c < channels; c++) {


				img.data[i++] = (T)(high_value * ((((x % (2 * period)) >= period) + ((y % (2 * period)) >= period)) % 2));
			}
		}
	}
	return img;
}

template <typename T>
void image_set_at(bitmap<T>* img, unsigned int x, unsigned int y, T val) {
	img->data[img->channels * (img->width * y + x)] = val;
}

template <typename T>
T image_get_at(bitmap<T> img, unsigned int x, unsigned int y, bool extend) {
	if (!vec_in_span<unsigned int, 2>({{ x,y } }, { { { 0u, 0u } }, { { img.width - 1, img.height - 1 } } })) {
		if (extend) {
			x = std::clamp<int>(x, 0, img.width - 1);
			y = std::clamp<int>(y, 0, img.height - 1);
			return img.data[img.channels * (img.width * y + x)];
		}
		return T(0);
	}
	return img.data[img.channels * (img.width * y + x)];
}

template <typename T>
void image_add_at(bitmap<T>* img, unsigned int x, unsigned int y, T val) {
	img->data[img->channels * (img->width * y + x)] += val;
}

template <typename T>
void super_impose_image(bitmap<T> source, bitmap<T>* destination, fvec2 offset) {
	//the offset is the position of the bottom left corner of the source image relative to the bottom left corner of the destintion image
	frect source_area = { {{0.0f, 0.0f}},{{(float)source.width, (float)source.height}} };
	frect destination_area = { {{0.0f, 0.0f}},{{(float)destination->width, (float)destination->height}} };//rect(vec2(0.0), vec2({ float(destination.width), float(destination.height) }));
	frect destination_overlap = expand_span_to_grid(span_overlap(source_area + offset, destination_area), { {0.0f,0.0f} }, { {1.0f, 1.0f} });
	frect source_overlap = expand_span_to_grid(destination_overlap - offset, { {0.0f,0.0f} }, { {1.0f, 1.0f} });
	/*
                 
		+------------------------------+								+------------------------------+	+------------------+-----------+	
		|           Source             |								|           Source             |	| Destination      |           |
		|                              |								|                              |	|   Overlap        |           |
		|           +------------------+-----------+					|           +------------------+	|                  |           |
		|           |     Overlap      |           |					|           | Source Overlap   |	|                  |           |
		|           |                  |           |					|           |                  |	+------------------+           |
		|           |                  |           |					|           |                  |	|                              |
		|           |                  |           |					|           |                  |	|                              |
		+-----------+------------------+           |	------------->	0-----------+------------------+	|                              |
	     ^ offset   |                              |														|                              |
	       \        |                              |														|          Destination         |
	         \      |                              |														0------------------------------+
	           \    |                              |														
	             \  |          Destination         |														
	               \+------------------------------+														
	*/

	//compute distribution from offset:
	/*		
		+----------+----------+
		|          |          |
		|          |          |      
		|     +----|----+     |      i   = (1-a)*(1-b)
		|(1-b)| iii| iv |     |      ii  = a*(1-b)	
		+----------+----------+      iii = (1-a)*b
		|   b | i  | ii |     |		 iv  = a*b
		|     +----|----+     |
		|        a |  (1-a)   |
		|          |          |
		+----------+----------+
	
	
	*/

	fvec2 modular_offset = mod_vec(offset, 1.0f);
	float a = abs(modular_offset[0]);
	float b = abs(modular_offset[1]);

	float pixel_weights[2][2] = 
	{ 
		{a * b,				(1.0f - a) * b	},
		{a * (1.0f - b),	(1.0f - a) * (1.0f - b)		}
	};

	int source_x = source_overlap.min[0];
	int source_y = source_overlap.min[1];
	for (int destination_x = destination_overlap.min[0]; destination_x < destination_overlap.max[0]; destination_x++) {
		source_y = source_overlap.min[0];
		for (int destination_y = destination_overlap.min[1]; destination_y < destination_overlap.max[1]; destination_y++) {

			float val = (
				pixel_weights[0][0] * image_get_at(source, source_x, source_y, false)+
				pixel_weights[1][0] * image_get_at(source, source_x + 1, source_y, false)+
				pixel_weights[0][1] * image_get_at(source, source_x, source_y + 1, false)+
				pixel_weights[1][1] * image_get_at(source, source_x + 1, source_y + 1, false)
			);
			image_add_at(destination, destination_x, destination_y, val);

			source_y++;
		}
		source_x++;
	}
};

template <typename T, typename U>
bitmap<T> bitmap_cast(bitmap<U> img, U max_in, T max_out) {
	bitmap<T> out;
	out.channels = img.channels;
	out.width = img.width;
	out.height = img.height;
	unsigned int size = img.width * img.height * img.channels;
	out.data = new T[size];
	for (int i = 0; i < size; i++) {
		out.data[i] = (T)(std::clamp(float(max_out) * float(img.data[i]) / float(max_in), 0.0f, float(max_out)));
	}
	return out;
}

void write_bitmap_to_png(const char* path, bitmap<unsigned char> img);

float compute_apparent_star_radius(star_camera camera, star_field_star star);
bitmap<float> generate_gaussian_image(unsigned int image_size, float scale, float standard_deviation);
star_field generate_star_field(star_field_generator gen, float start_time, float end_time);
float get_star_bounding_radius(star_field_star s);
unsigned short** simulate_image(int width, int height, int adc_bits);
