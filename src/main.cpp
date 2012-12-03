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
#define PI 3.14159265
#define POPULATION_SIZE 8
#define MAX_ITERATIONS 500
#define PARAMS_SIZE 2
#define UPPER_BOUND 1
#define LOWER_BOUND 0
#define C1 1
#define C2 1
#define MAX_SPEED 2
#define CONSTRICTION 0.729844

double inertia_bound[2] = { 0.4, 0.9 };
double w = 0.0;
double bounds_matrix[PARAMS_SIZE][2] = { { 0, 10 }, { 0, 10 } };

double global_best_solution[PARAMS_SIZE];
double global_best_fitness;

// métricas
double average = 0.0;
double variance = 0.0;
double standard_deviation = 0.0;

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
void calculate_metrics();

int main(int argc, char *argv[]) {
	clock_t time_start = clock();
	int iteration = 0;
	init();
	while (iteration < MAX_ITERATIONS) {
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
			temp += "X(" + number_to_String(i + 1) + ")=" + number_to_String(
					global_best_solution[i]) + " ";
		}
		cout << temp << endl;
		iteration++;
	}
	cout << "\nMelhor solução: f(x): " << global_best_fitness << endl;

	calculate_metrics();
	cout << "Média: " << average << endl;
	// cout << "Variância:" << variance << endl;
	cout << "Desvio padrão: " << standard_deviation << endl;

	string temp = "Valores: [";
	for (int i = 0; i < PARAMS_SIZE; i++) {
		temp += number_to_String(global_best_solution[i]);
		if (i != PARAMS_SIZE - 1) {
			temp += ", ";
		}
	}
	temp += "]";
	cout << temp << endl;

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
	for (int j = 0; j < PARAMS_SIZE; j++) {
		global_best_solution[j] = swarm[0].solution_array[PARAMS_SIZE];
	}
	global_best_fitness = swarm[0].fitness;
	get_global_best();
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
		/*
		 // padrão
		 swarm[index].speed_array[i] = swarm[index].speed_array[i]
		 + get_random_number() * (swarm[index].best_solution_array[i]
		 - swarm[index].solution_array[i]) + get_random_number()
		 * (global_best_solution[i] - swarm[index].solution_array[i]);
		 */

		/*
		 // coeficiente de inércia
		 swarm[index].speed_array[i] = w * swarm[index].speed_array[i]
		 + get_random_number() * (swarm[index].best_solution_array[i]
		 - swarm[index].solution_array[i]) + get_random_number()
		 * (global_best_solution[i] - swarm[index].solution_array[i]);
		 */

		// coeficiente de constricção
		swarm[index].speed_array[i] = CONSTRICTION
				* (swarm[index].speed_array[i] + get_random_number()
						* (swarm[index].best_solution_array[i]
								- swarm[index].solution_array[i])
						+ get_random_number() * (global_best_solution[i]
								- swarm[index].solution_array[i]));

		/*
		 // coeficiente de inércia + coef de fator individual e social
		 swarm[index].speed_array[i] = w * swarm[index].speed_array[i]
		 + C1 * get_random_number()
		 * (swarm[index].best_solution_array[i]
		 - swarm[index].solution_array[i])
		 + C2 * get_random_number()
		 * (global_best_solution[i]
		 - swarm[index].solution_array[i]);
		 */

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
	/*
	 // MIN f(x, y) = x^2 + y^2
	 swarm[index].fitness = pow(swarm[index].solution_array[0], 2) + pow(
	 swarm[index].solution_array[1], 2);
	 */

	/*
	 // MIN f(x,y) = x^2-x*y+y^2-3*y
	 swarm[index].fitness = pow(swarm[index].solution_array[0], 2)
	 - swarm[index].solution_array[0] * swarm[index].solution_array[1]
	 + pow(swarm[index].solution_array[1], 2) - 3
	 * swarm[index].solution_array[1];
	 */

	/*
	 // MIN f(x, y) = (x-2)^4 + (x - 2y)^2
	 swarm[index].fitness = pow(swarm[index].solution_array[0] - 2, 4)
	 + pow(
	 swarm[index].solution_array[0] - 2
	 * swarm[index].solution_array[1], 2);
	 */

	/*
	 // MIN f(x,y) = 100*(y-x^2)^2+(1 -x)^2
	 swarm[index].fitness = 100 * pow(
	 swarm[index].solution_array[1] - pow(
	 swarm[index].solution_array[0], 2), 2) + pow(
	 1 - swarm[index].solution_array[0], 2);
	 */

	/*
	 // MAX f(x) = x * sen(10*PI*x) + 1
	 swarm[index].fitness = swarm[index].solution_array[0] * sin(
	 10 * PI * swarm[index].solution_array[0]) + 1;
	 */

	// MIN f(x,y) = x*sen(4*x) + 1.1*y*sen(2*y)
	swarm[index].fitness = swarm[index].solution_array[0] * sin(
			4 * swarm[index].solution_array[0]) + 1.1
			* swarm[index].solution_array[1] * sin(
			2 * swarm[index].solution_array[1]);
}

void get_local_best_at(int index) {
	if (swarm[index].best_fitness > swarm[index].fitness) {
		for (int j = 0; j < PARAMS_SIZE; j++) {
			swarm[index].best_solution_array[j]
					= swarm[index].solution_array[j];
		}
		swarm[index].best_fitness = swarm[index].fitness;
	}
}

void get_global_best() {
	for (int i = 0; i < POPULATION_SIZE; i++) {
		//-------------------------------------------------- MAX > / MIN <
		if (swarm[i].fitness < global_best_fitness) {
			for (int j = 0; j < PARAMS_SIZE; j++) {
				global_best_solution[j] = swarm[i].solution_array[j];
			}
			global_best_fitness = swarm[i].fitness;
		}
	}
}

void calculate_metrics() {
	// Calcular a média
	double sum = 0;
	for (int i = 0; i < POPULATION_SIZE; i++) {
		sum += swarm[i].best_fitness;
	}
	average = (double) sum / (double) POPULATION_SIZE;
	// Calcuar a variância
	sum = 0;
	for (int i = 0; i < POPULATION_SIZE; i++) {
		sum += pow(swarm[i].best_fitness - average, 2);
	}
	variance = (double) sum / (double) POPULATION_SIZE;
	// Calculando o desvio padrão
	standard_deviation = pow(variance, 0.5);
}
