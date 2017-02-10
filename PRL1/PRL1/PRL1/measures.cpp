#include "measures.h"

double foo()
{
	return 1.0;
}

std::pair<double, double> measure_foo(double(*foo)())
{
	clock_t start, stop;
	start = clock();
	double pi = foo();
	stop = clock();
	//
	std::pair<double, double> result(pi, (double)(stop - start));
	return result;
}

void tell_results(std::pair<double, double> results) {
	double pi = results.first;
	double time = results.second;
	printf("Wartosc liczby PI wynosi %15.12f\n", pi);
	printf("Czas przetwarzania wynosi %f sekund\n\n", (time / 1000.0));
}