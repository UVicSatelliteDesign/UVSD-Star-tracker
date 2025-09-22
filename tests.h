#pragma once
#include "astrometry.h";

//create a synthetic star file. Then read is and check if the stars match the original data
void test_database_read_and_write();

void test_k_nearest_neighbor_search_accuracy(unsigned int k, unsigned int samples, unsigned int point_count);
void test_k_nearest_neighbor_search_performance(unsigned int k, unsigned int samples, unsigned int point_count);
void 