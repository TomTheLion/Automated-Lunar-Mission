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

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

double fcalc(double e, double dt, double f);

int main()
{
	double rl;
	double rr;
	double rp;
	double il;
	double ir;
	double dt;
	double re;
	double tf;
	double tfl;
	double tfr;
	double em;
	double im;
	double lpe;
	double fm;
	double fm0;
	double pi = 3.141592653589793;
	bool cflag;
	bool crash;
	vector<double> moon_op;
	vector<double> traj_co;
	Vector3d r;
	Vector3d v;
	Vector3d rm;
	Vector3d vm;
	Vector3d iaxis;
	Vector6d c;

	//cout << "\n";
	//cout << "UR3BP\n\n";
	clock_t begin_time = clock();

	cflag = true;
	crash = false;

	if (!file_to_vector("lunar_orbital_parameters.txt", moon_op))
		crash = true;
	if (!file_to_vector("trajectory_constraints.txt", traj_co))
		crash = true;

	if (crash)
	{
		cout << "Invalid inputs.\n";
	}
	else
	{
		rl = traj_co[0];
		rr = traj_co[1];
		rp = traj_co[2];
		il = traj_co[3];
		ir = traj_co[4];
		dt = traj_co[5];

		em = moon_op[0];
		im = moon_op[1];
		lpe = moon_op[2];
		fm0 = moon_op[3];

		re = (rl + rr) / 2.0;
		rl = pow(rl, 2.0);
		rr = pow(rr, 2.0);

		c << rl, rr, 1.0, 1.0, cos(il), cos(ir);
		fm = fcalc(em, dt, fm0);

		cr3bp_p(r, v, tf, rp);
		if (!cflag)
			return 2;
		ur3bp_i(r, v, rm, vm, iaxis, em, im, lpe, fm);

		tfl = tf;
		tfr = tf;

		ur3bp(r, v, rm, vm, iaxis, im, tfl, tfr, c, cflag);
		//ur3bp_w(r, v, rm, vm, tfl, tfr);

		write_orbital_parameters(r, v, tfl, 0.981889592516314, "orbital_parameters.txt");
		write_orbital_parameters(rm, vm, tfl, 1.0, "orbital_parameters_moon.txt");
		//cout << "\n";
		//cout << "Execution Time: ";
		system("PAUSE");
		//cout << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n\n";
		if (cflag)
		{
			return 0;
		}
		else
		{
			return 3;
		}

	}
	return 1;
}

double fcalc(double e, double dt, double f)
{
	double pi = 3.141592653589793;
	double E = 2 * atan2(sin(f / 2) * sqrt(1 - e), cos(f / 2) * sqrt(1 + e));
	double M = E - e * sin(E) + dt;

	while (abs(E - e * sin(E) - M) > 1e-8)
	{
		E = M + e * sin(E);
	}
	f = 2 * atan2(sin(E / 2) * sqrt(1 + e), cos(E / 2) * sqrt(1 - e));
	return f;
}