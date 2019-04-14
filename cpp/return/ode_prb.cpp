# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <vector>
# include <string>
# include <cmath>
# include <Eigen/Dense>
# include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

# include "ode.hpp"
# include "ur3bp.hpp"

int main()
{
	double pi = 3.141592653589793;
	bool cflag;
	bool crash;
	vector<double> inputs;
	vector<double> outputs;
	Vector3d r;
	Vector3d v;
	Vector3d rm;
	Vector3d vm;
	double p;
	double c;
	double ti;
	double tf;


	cout << "\n";
	cout << "UR3BP\n\n";
	clock_t begin_time = clock();

	cflag = true;
	crash = false;

	if (!file_to_vector("lunar_return_inputs.txt", inputs))
		crash = true;

	if (crash)
	{
		cout << "Invalid inputs.\n";
	}
	else
	{
		r << inputs[0], inputs[1], inputs[2];
		v << inputs[3], inputs[4], inputs[5];
		rm << inputs[6], inputs[7], inputs[8];
		vm << inputs[9], inputs[10], inputs[11];
		c = inputs[12];
		p = inputs[13];

		mrv(r, v, rm, vm, p, ti);
		//ur3bp_w(r, v, rm, vm, tf, tf);
		ur3bp_r(r, v, rm, vm, ti, tf, c, cflag);
		//ur3bp_w(r, v, rm, vm, tf, tf);

		//ofstream output_file;
		//output_file.open("lunar_return_outputs.txt");
		//for (size_t i = 0; i < outputs.size(); ++i)
		//{
		//	output_file << outputs[i];
		//}
		//output_file.close();
		cout << "\n";
		cout << "Execution Time: ";
		cout << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n\n";
		system("PAUSE");
		if (cflag)
		{
			return 0;
		}
		else
		{
			return 3;
		}

	}
	return 0;
}
