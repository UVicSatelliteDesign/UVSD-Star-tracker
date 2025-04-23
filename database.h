#pragma once
#include <iostream>

template <typename T>
struct database {
	unsigned int object_count;
	T* objects;
};


template <typename T>
database<T> load_database(const char* path) {
	FILE* source;
	fopen_s(&source, path, "rb");
	if (!source) {
		return {0, NULL};
	}
	
	//get number of elements in the file:
	unsigned int object_count;
	fread(&object_count, sizeof(unsigned int), 1, source);
	T* objects = new T[object_count];
	for (int i = 0; i < object_count; i++) {
		fread(&objects[i], sizeof(T), 1, source);
	}
	return {object_count, objects};
}

template <typename T>
void store_database(const char* path, T* objects, unsigned int object_count) {
	FILE* destination;
	fopen_s(&destination, path, "wb");
	if (!destination) {
		return;
	}

	fwrite(&object_count, sizeof(unsigned int), 1, destination);
	for (int i = 0; i < object_count; i++) {
		fwrite(objects + i, sizeof(T), 1, destination);
	}
	fclose(destination);
}