/*
 * PSO
 * main.cpp
 *
 *  Created on: Dec 1, 2012
 *      Author: aj.alves@zerokol.com
 */
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;

// BUG Eclipse
// #define CLOCKS_PER_SEC 1000000

/* Parâmetros do algoritmo */
#define POPULATION_SIZE 8
#define MAX_ITERATIONS 500
/* Parâmetros do problema */
#define PARAMS_SIZE 2
#define UPPER_BOUND 1
#define LOWER_BOUND 0
#define C1 1.8
#define C2 1.8
#define MAX_SPEED 2

double inertia_bound[2] = { 0.4, 0.9 };
double w = 0.0;
double bounds_matrix[PARAMS_SIZE][2] = { { -20, 40 }, { -30, 50 } };

double global_best_solution[PARAMS_SIZE];
double global_best_fitness;

struct particle {
	double solution_array[PARAMS_SIZE];
	double fitness;
	double best_solution_array[PARAMS_SIZE];
	double best_fitness;
	double speed_array[PARAMS_SIZE];
};

particle swarm[POPULATION_SIZE];

/* Funções auxiliares */
string number_to_String(double n);
double calculate_time(clock_t start, clock_t end);
void init();
double get_random_number();
double get_inertia_weight();
void calculate_speed_at(int index);
void calculate_displacement_at(int index);
void calculate_fitness_at(int index);
void get_local_best_at(int index);
void get_global_best();

int main(int argc, char *argv[]) {
	clock_t time_start = clock();
	int iteration = 0;
	init();
	while (iteration < MAX_ITERATIONS and global_best_fitness > 0.0001) {
		w = get_inertia_weight();
		for (int k = 0; k < POPULATION_SIZE; k++) {
			calculate_speed_at(k);
			calculate_displacement_at(k);
			calculate_fitness_at(k);
			get_local_best_at(k);
		}
		get_global_best();
		string temp = "Iteração(" + number_to_String(iteration) + ")-> ";
		for (int i = 0; i < PARAMS_SIZE; i++) {
			temp += "X(" + number_to_String(i + 1) + ")="
					+ number_to_String(global_best_solution[i]) + " ";
		}
		cout << temp << endl;
		iteration++;
	}
	cout << "f(x, y) = " << global_best_fitness << endl;
	cout << "Tempo de exec (PSO): " << calculate_time(time_start, clock())
			<< " ms" << endl;
	return 0;
}

string number_to_String(double n) {
	stringstream out;
	out << n;
	return out.str();
}

double calculate_time(clock_t start, clock_t end) {
	return 1000.0 * ((double) (end - start) / (double) CLOCKS_PER_SEC);
}

void init() {
	for (int i = 0; i < POPULATION_SIZE; i++) {
		// calcular aleatoriamente uma solução
		for (int j = 0; j < PARAMS_SIZE; j++) {
			swarm[i].solution_array[j] = get_random_number()
					* (bounds_matrix[j][UPPER_BOUND]
							- bounds_matrix[j][LOWER_BOUND])
					+ bounds_matrix[j][LOWER_BOUND];
		}
		// calcular aleatoriamente uma velocidade
		for (int j = 0; j < PARAMS_SIZE; j++) {
			swarm[i].speed_array[j] = MAX_SPEED * get_random_number();
		}
		// configurar a melhor solução como o atual
		for (int j = 0; j < PARAMS_SIZE; j++) {
			swarm[i].best_solution_array[j] = swarm[i].solution_array[j];
		}
		// calcular o fitness dessa particula
		calculate_fitness_at(i);
		// configurar o melhor fitness como o atual
		swarm[i].best_fitness = swarm[i].fitness;
	}
	// dentre as soluções aleatórias, pegar a melhor
	int flag = 0;
	for (int i = 1; i < POPULATION_SIZE; i++) {
		if (swarm[i].fitness < swarm[flag].fitness) {
			swarm[flag].fitness = swarm[i].fitness;
			flag = i;
		}
	}
	for (int j = 0; j < PARAMS_SIZE; j++) {
		global_best_solution[j] = swarm[flag].solution_array[PARAMS_SIZE];
	}
	global_best_fitness = swarm[flag].fitness;
}

double get_random_number() {
	return ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
}

double get_inertia_weight() {
	return (inertia_bound[UPPER_BOUND] - inertia_bound[LOWER_BOUND])
			* get_random_number() + inertia_bound[LOWER_BOUND];
}

void calculate_speed_at(int index) {
	for (int i = 0; i < PARAMS_SIZE; i++) {
		swarm[index].speed_array[i] = w * swarm[index].speed_array[i]
				+ C1 * get_random_number()
						* (swarm[index].best_solution_array[i]
								- swarm[index].solution_array[i])
				+ C2 * get_random_number()
						* (global_best_solution[i]
								- swarm[index].solution_array[i]);
		if (swarm[index].speed_array[i] > MAX_SPEED) {
			swarm[index].speed_array[i] = MAX_SPEED;
		}
	}
}

void calculate_displacement_at(int index) {
	for (int i = 0; i < PARAMS_SIZE; i++) {
		swarm[index].solution_array[i] += swarm[index].speed_array[i];
		if (swarm[index].solution_array[i] > bounds_matrix[i][UPPER_BOUND]) {
			swarm[index].solution_array[i] = 2 * bounds_matrix[i][UPPER_BOUND]
					- swarm[index].solution_array[i];
		}
		if (swarm[index].solution_array[i] < bounds_matrix[i][LOWER_BOUND]) {
			swarm[index].solution_array[i] = 2 * bounds_matrix[i][LOWER_BOUND]
					- swarm[index].solution_array[i];
		}
	}
}

void calculate_fitness_at(int index) {
	// MIN f(x,y) = 100 * (y - x^2)^2 + (x - 1)^2
	//swarm[index].fitness = 100
	//		* (swarm[index].solution_array[1]
	//				- swarm[index].solution_array[0]
	//						* swarm[index].solution_array[0])
	//		* (swarm[index].solution_array[1]
	//				- swarm[index].solution_array[0]
	//						* swarm[index].solution_array[0])
	//		+ (swarm[index].solution_array[0] - 1)
	//				* (swarm[index].solution_array[0] - 1);
	// MIN f(x,y) = 100*(y-x^2)^2+(1 -x)^2
	swarm[index].fitness = 100
			* pow(
					swarm[index].solution_array[1]
							- pow(swarm[index].solution_array[0], 2), 2)
			+ pow(1 - swarm[index].solution_array[0], 2);
}

void get_local_best_at(int index) {
	if (swarm[index].best_fitness > swarm[index].fitness) {
		for (int j = 0; j < PARAMS_SIZE; j++) {
			swarm[index].best_solution_array[j] =
					swarm[index].solution_array[j];
		}
		swarm[index].best_fitness = swarm[index].fitness;
	}
}

void get_global_best() {
	for (int i = 0; i < POPULATION_SIZE; i++) {
		if (swarm[i].fitness < global_best_fitness) {
			for (int j = 0; j < PARAMS_SIZE; j++) {
				global_best_solution[j] = swarm[i].solution_array[j];
			}
			global_best_fitness = swarm[i].fitness;
		}
	}
}
