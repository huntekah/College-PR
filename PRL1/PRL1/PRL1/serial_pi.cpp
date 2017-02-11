#include <stdio.h>

#include <time.h>
#include <iostream>
#include "measures.h"
#include <omp.h>

long long num_steps = 1000000000;
double step;

double find_serial_pi() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;

	double x, pi, sum = 0.0;
	double step = 0.0;
	int i;
	step = 1. / (double)steps;
	for (i = 0; i<steps; i++)
	{
		x = (i + .5)*step;
		sum = sum + 4.0 / (1. + x*x);
	}

	pi = sum*step;
	return pi;
}

//data races on 'sum'
double find_parallel_pi_data_races() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;

	double x, pi, sum = 0.0;
	double step = 0.0;
	int i;
	step = 1. / (double)steps;

#pragma omp parallel for private(x) shared(sum)
	for (i = 0; i<steps; i++)
	{
		x = (i + .5)*step;
		sum = sum + 4.0 / (1. + x*x);
	}

	pi = sum*step;
	return pi;
}

double find_parallel_pi_atomic() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;

	double x, pi, sum = 0.0;
	double step = 0.0;
	int i;
	step = 1. / (double)steps;

#pragma omp parallel for private(x) shared(sum)
	for (i = 0; i<steps; i++)
	{
		x = (i + .5)*step;
		#pragma omp atomic
		sum += 4.0 / (1. + x*x);
	}

	pi = sum*step;
	return pi;
}

double find_parallel_pi_reduction() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;

	double x, pi, sum = 0.0;
	double step = 0.0;
	int i;
	step = 1. / (double)steps;

#pragma omp parallel for private(x) reduction(+:sum)
	for (i = 0; i<steps; i++)
	{
		x = (i + .5)*step;
		sum += 4.0 / (1. + x*x);
	}

	pi = sum*step;
	return pi;
}

double find_parallel_pi_arrays() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;
	//steps = 16;

	double x, pi, sum = 0.0;
	double step = 0.0;
	int i;
	step = 1. / (double)steps;

	double* sums;

#pragma omp parallel private(x) shared(sums)
	{
#pragma omp single
	{
	sums = new double[omp_get_num_threads()];
	for (int i = 0; i < omp_get_num_threads(); i++) {
		sums[i] = 0.0;
	}
	}
#pragma omp for
	for (i = 0; i < steps; i++)
	{
		x = (i + .5)*step;
		sums[omp_get_thread_num()] += 4.0 / (1. + x*x);

	}

	}
	for (int i = 0; i < omp_get_max_threads(); i++) {
		sum += sums[i];
	}
	pi = sum *step;
	delete [] sums;
	return pi;
}

double find_parallel_pi_local_sum() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;
	//steps = 16;

	int i;
	double x = 0.0, pi = 0.0, sum = 0.0;
	step = 1. / (double)steps;
#pragma omp parallel
	{
		// Lokalna suma
		double lsum = 0;
		// Kazdy watek dostaje część iteracji do zrobienia, kazdy pracuje na swojej lokalnej sumie
#pragma omp for private(x)
		for (i = 0; i < steps; i++)
		{
			x = (i + .5)*step;
			lsum += 4.0 / (1. + x*x);
		}

#pragma omp atomic
		// Scalenie lokalnych lsum do zmiennej globalnej sum
		sum += lsum;
	}
	pi = sum*step;

	return pi;
}

double find_parallel_pi_volatile() {
	printf("function: %s\n", __FUNCTION__);

	int steps = num_steps;
	//4 wątki
	omp_set_num_threads(4);
	double x = 0.0, pi = 0.0, sum = 0.0;
	// Wskazówka dla kompilatora, że jak jest odwołanie do zmiennej sums[id] ma być za kazdym razem pobierana i modyfikowana w pamięci
	// Efekt spowolnienia będzie do zaobserwowania:
	volatile double sums[4] = { 0,0,0,0 };
	int i;
	step = 1. / (double)steps;

	// Kazdy z watkow dopisuje do swojej prywatnej sumy - w osobnej komore tablicy

#pragma omp parallel
	{
		// Pobierz ID

		int id = omp_get_thread_num();

#pragma omp for private(x)
		for (i = 0; i < steps; i++)
		{
			// Pobierz do rejestru sums[i]
			x = (i + .5)*step;
			sums[id] += 4.0 / (1. + x*x);
		}


	}
	sum = sums[0] + sums[1] + sums[2] + sums[3];

	pi = sum*step;

	return pi;
}


int main(int argc, char* argv[])
{
	std::pair<double, double> result = measure_foo(&find_serial_pi);
	tell_results(result);
	
	result = measure_foo(&find_parallel_pi_data_races);
	tell_results(result);

	result = measure_foo(&find_parallel_pi_atomic);
	tell_results(result);

	result = measure_foo(&find_parallel_pi_reduction);
	tell_results(result);
	
	result = measure_foo(&find_parallel_pi_arrays);
	tell_results(result);

	result = measure_foo(&find_parallel_pi_local_sum);
	tell_results(result);

	result = measure_foo(&find_parallel_pi_volatile);
	tell_results(result);

	std::system("Pause");
	return 0;
}