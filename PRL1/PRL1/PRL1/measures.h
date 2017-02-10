#pragma once

#include <stdio.h>
#include <time.h>
#include <utility>
double foo();

std::pair<double,double> measure_foo(double(*foo)());

void tell_results(std::pair<double, double>);