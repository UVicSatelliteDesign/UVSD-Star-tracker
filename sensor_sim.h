#pragma once
#include "linear_algebra.h"
#include <algorithm>
#include "stb_image_write.h"

struct star_field_centroid {
	unsigned int id;
	vec<2> centroid;
};
struct star_field_star {
	vec<3> direction;
	float magnitude;
	float bounding_radius;
};

template <typename T>
struct image {
	unsigned int width;
	unsigned int height;
	unsigned int channels;
	T* data;
};
struct star_field {
	image<short> image;
	unsigned int star_count;
	star_field_centroid* centroids;
	matrix<3, 3> orientation;
};
struct star_camera {
	float aperture;//lens diameter
	float fov;
	float aspect_ratio;
	float focal_length;
	float quantum_efficeincy;
	float angular_spread_parameter;
	float transmittance;
	float pixel_area;
	float exposure_time;//amount of time collecting light
	float frame_interval;//amound of time between frames (realistically longer than the exposure time)
	float sensor_size;//height of the sensor
	unsigned int width;//image width in pixels
	unsigned int height;//image height in pixels
	unsigned int full_well;//number of electrons required to saturate a pixel
};
struct star_field_generator {
	star_camera camera;
	star_field_star* stars;
	matrix<3, 3>(*path)(float);
	float* guassian_template;
};
star_camera default_star_camera();

template <typename T>
image<T> initialize_image(unsigned int width, unsigned height, unsigned int channels, T initial_value) {
	image<T> img;
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
image<T> initialize_checkerboard(unsigned int width, unsigned int height, unsigned int channels, T high_value, unsigned int period) {
	image<T> img;
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
void image_set_at(image<T> img, unsigned int x, unsigned int y, T val) {
	img.data[img.channels * (img.width * y + x)] = val;
}
template <typename T>
T image_get_at(image<T> img, unsigned int x, unsigned int y, bool extend) {
	if (!vec_in_span(uvec2({ x,y }), urect(uvec2({ 0, 0 }), uvec2({ img.width, img.height })))) {
		if (extend) {
			x = std::clamp<int>(x, 0, img.width - 1);
			y = std::clamp<int>(y, 0, img.height - 1);
			return img.data[img.channels * (img.width * y + x)];
		}
		return T(0);
	}
	x = std::clamp<int>(x, 0, img.width - 1);
	y = std::clamp<int>(y, 0, img.height - 1);
	return img.data[img.channels * (img.width * y + x)];
}
template <typename T>
void image_add_at(image<T> img, unsigned int x, unsigned int y, T val) {
	img.data[img.channels * (img.width * y + x)] += val;
}

template <typename T>
void super_impose_image(image<T> source, image<T> destination, vec2 offset) {
	//the offset is the position of the bottom left corner of the source image relative to the bottom left corner of the destintion image
	rect source_area = offset_span(rect(vec2(0.0), vec2({ float(source.width), float(source.height) })), offset);
	rect destination_area = rect(vec2(0.0), vec2({ float(destination.width), float(destination.height) }));
	rect overlap = clamp_rect(source_area, destination_area);
	rect working_area = expand_span_to_grid(overlap, vec2(0.0), 1.0f);//irect(offset_span(overlap, negate_vector(offset)));
	/*
                 
		+------------------------------+
		|           Source             |
		|                              |
		|           +------------------+-----------+
		|           |     Overlap      |           |
		|           |                  |           |
		|           |                  |           |
		|           |                  |           |
		+-----------+------------------+           |
	     ^ offset   |                              |
	       \        |                              |
	         \      |                              |
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
		|   b | iii| iv |     |      ii  = a*(1-b)	
		+----------+----------+      iii = (1-a)*b
		|(1-b)| i  | ii |     |		 iv  = a*b
		|     +----|----+     |
		|    (1-a) |  a       |
		|          |          |
		+----------+----------+
	
	
	*/
	vec2 modular_offset = mod_vector(offset, 1.0f);
	float a = abs(modular_offset[0]);
	float b = abs(modular_offset[1]);
	float pixel_weights[2][2] = 
	{ 
		{(1.0 - a) * (1.0 - b),		a * (1.0 - b)	}, 
		{(1.0 - a) * b,				a * b		}
	};

	ivec2 source_offsets[2][2] =
	{
		{ivec2({(int)floor(-offset[0]), (int)floor(-offset[1])}),		ivec2({(int)floor(-offset[0]) + 1, (int)floor(-offset[1])})},
		{ivec2({(int)floor(-offset[0]), (int)floor(-offset[1]) + 1}),	ivec2({(int)floor(-offset[0]) + 1, (int)floor(-offset[1]) + 1})}
	};

	for (int x = working_area.min[0]; x < working_area.max[0]; x++) {
		for (int y = working_area.min[1]; y < working_area.max[1]; y++) {
			float val = (
				pixel_weights[0][0] * image_get_at(source, x + source_offsets[0][0][0], y + source_offsets[0][0][1], false) +
				pixel_weights[0][1] * image_get_at(source, x + source_offsets[0][1][0], y + source_offsets[0][1][1], false) +
				pixel_weights[1][0] * image_get_at(source, x + source_offsets[1][0][0], y + source_offsets[1][0][1], false) +
				pixel_weights[1][1] * image_get_at(source, x + source_offsets[1][1][0], y + source_offsets[1][1][1], false)
			);
			image_add_at(destination, x, y, (T)val);
		}
	}
};

template <typename T, typename U>
image<T> image_cast(image<U> img, U max_in, T max_out) {
	image<T> out;
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

void write_image_to_png(const char* path, image<unsigned char> img);

float compute_apparent_star_radius(star_camera camera, star_field_star star);
void compute_guassian_template(star_field_generator gen);
star_field generate_star_field(star_field_generator gen, float time);

unsigned short** simulate_image(int width, int height, int adc_bits);
