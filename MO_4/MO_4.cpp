#include <iostream>
#include "RandomSearch.h"

int main()
{
	SimpleRandomSearch method;
	method.input();
	int MethodType = 0;
	std::cout << "Enter method(1 - Simple search; 2 - Method hyperquadrate; 3 - algorithm_1; 4 - algorithm_2; 5 - algorithm_3):\n";
	std::cin >> MethodType;

	if (MethodType == 1)
	{
		method.log_file.open("simple.txt");
		method.Search();
		std::cout << "Simple search:\n";
		std::cout << "f(x,y): " << method.func_0 << '\n';
		std::cout << "x: " << method.x[0] << '\n' << "y: " << method.x[1] << '\n';
	}
	else if (MethodType == 2)
	{
		method.hyperquadrate();
		std::cout << "Method hyperquadrate:\n";
		std::cout << "f(x,y): " << method.func_0 << '\n';
		std::cout << "x: " << method.x[0] << '\n' << "y: " << method.x[1] << '\n';
		std::cout << "Count_use_function: " << method.func.count_use_f << '\n';
	}
	else if (MethodType == 3)
	{
		method.log_file.open("alg1.txt");
		method.algorithm_1();
		method.log_file.close();

		std::cout << "algorithm_1:\n";
		std::cout << "f(x,y): " << method.func_0 << '\n';
		std::cout << "x: " << method.x[0] << '\n' << "y: " << method.x[1] << '\n';
		std::cout << "Count_use_function: " << method.func.count_use_f << '\n';
	}
	else if (MethodType == 4)
	{
		method.log_file.open("alg2.txt");
		method.algorithm_2();
		method.log_file.close();

		std::cout << "algorithm_2:\n";
		std::cout << "f(x,y): " << method.func_0 << '\n';
		std::cout << "x: " << method.x[0] << '\n' << "y: " << method.x[1] << '\n';
		std::cout << "Count_use_function: " << method.func.count_use_f << '\n';
	}
	else if (MethodType == 5)
	{
		method.log_file.open("alg3.txt");
		method.algorithm_3();
		method.log_file.close();

		std::cout << "algorithm_3:\n";
		std::cout << "f(x,y): " << method.func_0 << '\n';
		std::cout << "x: " << method.x[0] << '\n' << "y: " << method.x[1] << '\n';
		std::cout << "Count_use_function: " << method.func.count_use_f << '\n';
	}
	else
	{
		std::cout << "Not correct method:\n";
		return 0;
	}

}
