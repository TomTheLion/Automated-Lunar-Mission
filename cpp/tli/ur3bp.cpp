# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <vector>
# include <cmath>
# include <Eigen/Dense>
# include <Eigen/Geometry>

// earth / moon
//#define m1 0.0121532665607235
//#define m2 0.9878467334392770

// kerbin / mun
#define m1 0.018110407483686
#define m2 0.981889592516314

using namespace std;
using namespace Eigen;

# include "ode.hpp"
# include "ur3bp.hpp"

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

bool file_to_vector(string filename, vector<double> &vector)
{
	ifstream in(filename.c_str());
	string str;

	if (!in)
	{
		cerr << "Cannot open the File : " << filename << endl;
		return false;
	}

	while (getline(in, str))
	{
		if (str.size() > 0)
			vector.push_back(stod(str));
	}

	in.close();
	return true;
}

double safe_acos(double x)
{
	if (x > 1)
		return 0;
	else if (x < -1)
		return 3.141592653589793;
	else
		return acos(x);
}

void print_format(string label, double value)
{
	int m;
	int p;

	m = (int)log10(value);
	p = 12;

	if (m > 0)
		p -= m;

	cout << setprecision(p);
	cout << fixed;
	cout << "\n";
	cout << label;
	cout << value;
	return;
}

void y_init_cr3bp(Vector3d r, Vector3d v, double y[])
{
	int i;
	int j;

	y[0] = r(0);
	y[1] = r(1);
	y[2] = r(2);
	y[3] = v(0);
	y[4] = v(1);
	y[5] = v(2);

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (i == j)
				y[j + 6 * i + 6] = 1;
			else
				y[j + 6 * i + 6] = 0;
		}
	}
	return;
}

void y_init_ur3bp(Vector3d r, Vector3d v, Vector3d rm, Vector3d vm, double y[])
{
	int i;
	int j;

	y[0] = rm(0);
	y[1] = rm(1);
	y[2] = rm(2);
	y[3] = vm(0);
	y[4] = vm(1);
	y[5] = vm(2);

	y[6] = r(0);
	y[7] = r(1);
	y[8] = r(2);
	y[9] = v(0);
	y[10] = v(1);
	y[11] = v(2);

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (i == j)
				y[j + 6 * i + 12] = 1;
			else
				y[j + 6 * i + 12] = 0;
		}
	}
	return;
}

void write_orbital_parameters(Vector3d r, Vector3d v, double dt, double mu, string file_name)
{
	double r_norm;
	double v_norm;
	double h_norm;
	double e_norm;
	double n_norm;
	double lpe;
	double sma;
	double i;
	double lan;
	double pi = 3.141592653589793;
	Vector3d k;
	Vector3d h;
	Vector3d n;
	Vector3d e;
	
	k = Vector3d::UnitZ();
	r_norm = r.norm();
	v_norm = v.norm();
	h = r.cross(v);
	h_norm = h.norm();
	i = safe_acos(h(2) / h_norm);
	n = k.cross(h);
	n_norm = n.norm();
	e = 1 / mu * (v.cross(h) - mu * r / r_norm);
	e_norm = e.norm();
	sma = pow(h_norm, 2) / (1 - pow(e_norm, 2)) / mu;

	if (n_norm == 0)
		lan = 0;
	else
	{
		lan = safe_acos(n(0) / n_norm);
		if (n(1) < 0)
			lan = 2 * pi - lan;
	}

	if (e_norm == 0)
		lpe = 0;
	else
	{
		if (n_norm == 0)
			lpe = atan2(e(1), e(0));
		else
		{
			lpe = safe_acos(n.dot(e) / (n_norm * e_norm));
			if (e(2) < 0) { lpe = 2 * pi - lpe; }
		}
	}

	ofstream output_file;
	output_file.open(file_name);
	output_file << setprecision(16) << sma << '\n' << e_norm << '\n' << i * 180 / pi << '\n' << lan * 180 / pi << '\n' << lpe * 180 / pi << '\n' << dt << '\n';
	return;
}

void dcr3bp(double t, double y[], double yp[])
{
	int i;
	double x1 = y[0] + m1;
	double x2 = y[0] - m2;
	double d2 = pow(x1, 2) + pow(y[1], 2) + pow(y[2], 2);
	double d3 = d2 * sqrt(d2);
	double d5 = d2 * d3;
	double r2 = pow(x2, 2) + pow(y[1], 2) + pow(y[2], 2);
	double r3 = r2 * sqrt(r2);
	double r5 = r2 * r3;
	double g1, g2, g3, g4, g5, g6, g7, g8, g9;

	yp[0] = y[3];
	yp[1] = y[4];
	yp[2] = y[5];
	yp[3] = y[0] + 2 * y[4] - m2 * x1 / d3 - m1 * x2 / r3;
	yp[4] = y[1] - 2 * y[3] - y[1] * (m2 / d3 + m1 / r3);
	yp[5] = -y[2] * (m2 / d3 + m1 / r3);

	g1 = 1 - m2 / d3 - m1 / r3 + 3 * m2 * x1 * x1 / d5 + 3 * m1 * x2 * x2 / r5;
	g2 = 3 * m2 * x1 * y[1] / d5 + 3 * m1 * x2 * y[1] / r5;
	g3 = 3 * m2 * x1 * y[2] / d5 + 3 * m1 * x2 * y[2] / r5;
	g4 = g2;
	g5 = 1 - m2 / d3 - m1 / r3 + 3 * m2 * y[1] * y[1] / d5 + 3 * m1 * y[1] * y[1] / r5;
	g6 = 3 * m2 * y[1] * y[2] / d5 + 3 * m1 * y[1] * y[2] / r5;
	g7 = g3;
	g8 = g6;
	g9 = 1 - m2 / d3 - m1 / r3 + 3 * m2 * y[2] * y[2] / d5 + 3 * m1 * y[2] * y[2] / r5;

	for (i = 0; i < 18; i++)
	{
		yp[i + 6] = y[i + 24];
	}
	for (i = 0; i < 6; i++)
	{
		yp[i + 24] = g1 * y[i + 6] + g2 * y[i + 12] + g3 * y[i + 18] + 2 * y[i + 30];
		yp[i + 30] = g4 * y[i + 6] + g5 * y[i + 12] + g6 * y[i + 18] - 2 * y[i + 24];
		yp[i + 36] = g7 * y[i + 6] + g8 * y[i + 12] + g9 * y[i + 18];
	}
	return;
}

void dur3bp(double t, double y[], double yp[])
{
	const int n = 3;
	int l = 1;
	int k = n - 1;
	double mu[n]{ m1, 0, m2 };

	int m = 6;
	int km = k * m;
	int kp1m = (k + 1) * m;
	int kp2m = (k + 2) * m;
	int kp3m = (k + 3) * m;
	int kp4m = (k + 4) * m;
	int kp5m = (k + 5) * m;

	int i;
	int im;
	int ij;
	int j;
	int jm;
	int ji;
	int ki;
	int kj;
	int kl;
	int lj;

	double r2[n * n];
	double r3[n * n];
	double rx[n * n];
	double ry[n * n];
	double rz[n * n];
	double ypk[3]{ 0, 0, 0 };

	double rkl3;
	double rkl5;
	double rxkl;
	double rykl;
	double rzkl;

	double a41;
	double a42;
	double a43;
	double a52;
	double a53;
	double a63;

	double temp1;
	double temp2;
	double temp3;
	double temp4;
	double temp5;
	double temp6;

	for (i = 0; i < n; i++)
	{
		im = i * m;
		for (j = 0; j < n; j++)
		{
			jm = j * m;
			ij = j + i * n;
			if (j > i)
			{
				if (j != k)
				{
					rx[ij] = y[jm + 0] - y[im + 0];
					ry[ij] = y[jm + 1] - y[im + 1];
					rz[ij] = y[jm + 2] - y[im + 2];
				}
				else
				{
					rx[ij] = -y[im + 0];
					ry[ij] = -y[im + 1];
					rz[ij] = -y[im + 2];
				}
				r2[ij] = rx[ij] * rx[ij] + ry[ij] * ry[ij] + rz[ij] * rz[ij];
				r3[ij] = std::sqrt(r2[ij]) * r2[ij];
			}
			else if (j < i)
			{
				ji = i + j * n;
				rx[ij] = -rx[ji];
				ry[ij] = -ry[ji];
				rz[ij] = -rz[ji];
				r2[ij] = r2[ji];
				r3[ij] = r3[ji];
				if (i == k)
				{
					kj = j + k * n;
					temp1 = -mu[j] / r3[kj];
					ypk[0] += rx[kj] * temp1;
					ypk[1] += ry[kj] * temp1;
					ypk[2] += rz[kj] * temp1;
				}
			}
		}
	}

	for (i = 0; i < n - 1; i++)
	{
		im = i * m;
		ki = i + k * n;
		temp2 = -mu[k] / r3[ki];
		yp[im + 0] = y[im + 3];
		yp[im + 1] = y[im + 4];
		yp[im + 2] = y[im + 5];
		yp[im + 3] = temp2 * rx[ki] + ypk[0];
		yp[im + 4] = temp2 * ry[ki] + ypk[1];
		yp[im + 5] = temp2 * rz[ki] + ypk[2];
		for (j = 0; j < n - 1; j++)
		{
			if (j != i)
			{
				jm = j * m;
				ij = j + i * n;
				temp3 = mu[j] / r3[ij];
				yp[im + 3] += temp3 * rx[ij];
				yp[im + 4] += temp3 * ry[ij];
				yp[im + 5] += temp3 * rz[ij];
			}
		}
	}

	temp4 = mu[l] + mu[k];
	kl = l + k * n;
	rkl3 = r3[kl];
	rkl5 = r2[kl] * rkl3;
	rxkl = rx[kl];
	rykl = ry[kl];
	rzkl = rz[kl];

	a41 = temp4 * (3 * rxkl * rxkl / rkl5 - 1 / rkl3);
	a52 = temp4 * (3 * rykl * rykl / rkl5 - 1 / rkl3);
	a63 = temp4 * (3 * rzkl * rzkl / rkl5 - 1 / rkl3);
	a42 = temp4 * (3 * rxkl * rykl / rkl5);
	a43 = temp4 * (3 * rxkl * rzkl / rkl5);
	a53 = temp4 * (3 * rykl * rzkl / rkl5);

	for (j = 0; j < n - 1; j++)
	{
		if (j != l)
		{
			lj = j + l * n;
			temp5 = 3 * mu[j] / (r2[lj] * r3[lj]);
			temp6 = mu[j] / r3[lj];
			a41 += temp5 * rx[lj] * rx[lj] - temp6;
			a52 += temp5 * ry[lj] * ry[lj] - temp6;
			a63 += temp5 * rz[lj] * rz[lj] - temp6;
			a42 += temp5 * rx[lj] * ry[lj];
			a43 += temp5 * rx[lj] * rz[lj];
			a53 += temp5 * ry[lj] * rz[lj];
		}
	}

	for (i = 0; i < 18; i++)
	{
		yp[i + km] = y[i + kp3m];
	}
	for (i = 0; i < 6; i++)
	{
		yp[i + kp3m] = a41 * y[i + km] + a42 * y[i + kp1m] + a43 * y[i + kp2m];
		yp[i + kp4m] = a42 * y[i + km] + a52 * y[i + kp1m] + a53 * y[i + kp2m];
		yp[i + kp5m] = a43 * y[i + km] + a53 * y[i + kp1m] + a63 * y[i + kp2m];
	}
	return;
}

void cr3bp_p(Vector3d &r, Vector3d &v, double &tf, double rp)
{
	double energy;
	double a;
	double e;
	double E;
	double f;
	double F;
	double h;
	double M;
	double vp;
	double theta;
	double r_norm;
	double rx;
	double ry;
	double vx;
	double vy;
	double rs = pow(m1 / m2, 0.4);
	double pi = 3.141592653589793;

	energy = 0.5 - m1 / rs;
	vp = pow(2.0 * (energy + m1 / rp), 0.5);
	h = rp * vp;
	a = -0.5 * m1 / energy;
	e = 1.0 - rp / a;
	f = safe_acos((a / rs * (1.0 - e * e) - 1.0) / e);
	F = 2.0 * atanh(sqrt((e - 1.0) / (e + 1.0)) * tan(f / 2.0));
	M = e * sinh(F) - F;
	tf = pow(h, 3.0) / pow(m1, 2.0) / pow(e * e - 1.0, 1.5) * M;
	theta = pi / 2.0 - f + safe_acos(h / rs);

	vx = -cos(theta);
	vy = 1.0 - sin(theta);
	rx = 1.0 + rs * cos(f);
	ry = -rs * sin(f);
	r_norm = sqrt(rx * rx + ry * ry);
	energy = (vx * vx + vy * vy) / 2.0 - m2 / r_norm;
	h = rx * vy - ry * vx;
	e = sqrt(1.0 + 2.0 * energy * h * h / m2 / m2);
	a = -m2 / 2.0 / energy;
	f = safe_acos((a / r_norm * (1.0 - e * e) - 1.0) / e);
	E = 2.0 * atan(tan(f / 2.0) * sqrt((1.0 - e) / (1.0 + e)));
	M = E - e * sin(E);
	tf += M * sqrt(pow(a, 3.0) / m2);
	r << 1.0 - m1 + rp, 0.0, 0.0;
	v << 0.0, -vp, 0.0;
	return;
}

void ur3bp_i(Vector3d &r, Vector3d &v, Vector3d &rm, Vector3d &vm, Vector3d &iaxis, double em, double im, double lpe, double fm)
{
	double rm_norm;
	double vm_norm;
	double hm_norm;
	double phim;
	double vsign;
	double pi = 3.141592653589793;
	Vector3d hm;
	Vector3d zaxis;

	if (fm < pi)
		vsign = 1;
	else
		vsign = -1;
	
	hm_norm = sqrt(1.0 - em * em);;
	rm_norm = hm_norm * hm_norm / (1.0 + em * cos(fm));
	vm_norm = sqrt(2 / rm_norm - 1);
	phim = safe_acos(hm_norm / (rm_norm * vm_norm));

	rm << rm_norm, 0, 0;
	vm << vsign * vm_norm * sin(phim), vm_norm * cos(phim), 0;
	zaxis = Vector3d::UnitZ();
	iaxis = Vector3d::UnitX();

	r(0) -= (1 - m1);
	r(0) += rm(0);
	v(1) += vm(1);

	iaxis = AngleAxisd(-(lpe + fm), zaxis) * iaxis;
	rm = AngleAxisd(im, iaxis) * rm;
	vm = AngleAxisd(im, iaxis) * vm;
	r = AngleAxisd(im, iaxis) * r;
	v = AngleAxisd(im, iaxis) * v;
	return;
}

void ur3bp(Vector3d &r, Vector3d &v, Vector3d &rm, Vector3d &vm, Vector3d iaxis, double im, double &tfl, double &tfr, Vector6d c, bool &cflag)
{
	bool phase1;
	int iflag;
	int iter;
	int iwork1[5];
	int iwork2[5];
	int neqn = 48;
	int wi_yp = 100 + neqn * 4 - 1;
	double abserr;
	double relerr;
	double t1;
	double t2;
	double tout1;
	double tout2;
	double rm_norm;
	double rp_norm;
	double vp_norm;
	double raxy;
	double raxz;
	double vaxy;
	double inc_mul = 0.0;
	double beta = 1.0;
	double fn;
	double fnold;
	double pi = 3.141592653589793;
	double *work1;
	double *work2;
	double *y1;
	double *y2;
	Vector3d rl_fail;
	Vector3d rr_fail;
	Vector3d vl_fail;
	Vector3d vr_fail;
	Matrix6d df;
	Vector6d f;
	Vector6d x;
	Vector6d xold;

	phase1 = true;
	iter = 0;
	abserr = 1e-7;
	relerr = 0.0;

	rm_norm = rm.norm();
	rp_norm = r.norm() - rm_norm;
	vp_norm = v.norm() ;
	raxy = 0.0;
	raxz = -1e-4;
	vaxy = 0.0;

	y1 = new double[neqn];
	y2 = new double[neqn];
	work1 = new double[100 + 21 * neqn];
	work2 = new double[100 + 21 * neqn];

	while (true)
	{
		iter++;
		t1 = 0.0;
		t2 = 0.0;
		tout1 = tfr;
		tout2 = tfl;
		y_init_ur3bp(r, v, rm, vm, y1);
		y_init_ur3bp(r, v, rm, vm, y2);

		iflag = 1;
		ode(dur3bp, neqn, y1, t1, tout1, relerr, abserr, iflag, work1, iwork1);
		iflag = -1;
		ode(dur3bp, neqn, y2, tout2, t2, relerr, abserr, iflag, work2, iwork2);

		ddd(y2 + 6, work2 + wi_yp + 6, y1 + 6, work1 + wi_yp + 6, df, f, c, rp_norm, vp_norm, im, raxy, raxz, vaxy, v, iaxis, inc_mul);

		fn = f.norm();

		if (iter == 1 && iflag != 2)
		{
			cout << "CR3BP - Fatal error - ODE returned IFLAG = " << iflag << " on first iteration\n";
			cflag = false;
			break;
		}
		else if (iter > 1 && (iflag != 2 || fn > fnold))
		{
			beta *= 0.25;
			vp_norm = (1.0 - beta) * xold(0) + beta * x(0);
			vaxy = (1.0 - beta) * xold(1) + beta * x(1);
			raxy = (1.0 - beta) * xold(2) + beta * x(2);
			raxz = (1.0 - beta) * xold(3) + beta * x(3);
			tfl = (1.0 - beta) * xold(4) + beta * x(4);
			tfr = (1.0 - beta) * xold(5) + beta * x(5);

			r(0) = rp_norm * cos(raxy) * cos(raxz) + rm_norm;
			r(1) = rp_norm * sin(raxy) * cos(raxz);
			r(2) = rp_norm * sin(raxz);

			v(0) = sin(raxy) * cos(vaxy) + cos(raxy) * sin(raxz) * sin(vaxy);
			v(1) = -cos(raxy) * cos(vaxy) + sin(raxy) * sin(raxz) * sin(vaxy);
			v(2) = -cos(raxz) * sin(vaxy);

			v *= vp_norm;

			r = AngleAxisd(im, iaxis) * r;
			v = AngleAxisd(im, iaxis) * v;
		}
		else
		{
			beta = 1.0;
			fnold = f.norm();
			xold << vp_norm, vaxy, raxy, raxz, tfl, tfr;

			x = df.fullPivLu().solve(df * xold - f);

			vp_norm = x(0);
			vaxy = x(1);
			raxy = x(2);
			raxz = x(3);
			tfl = x(4);
			tfr = x(5);

			r(0) = rp_norm * cos(raxy) * cos(raxz) + rm_norm;
			r(1) = rp_norm * sin(raxy) * cos(raxz);
			r(2) = rp_norm * sin(raxz);

			v(0) = sin(raxy) * cos(vaxy) + cos(raxy) * sin(raxz) * sin(vaxy);
			v(1) = -cos(raxy) * cos(vaxy) + sin(raxy) * sin(raxz) * sin(vaxy);
			v(2) = -cos(raxz) * sin(vaxy);

			v *= vp_norm;

			r = AngleAxisd(im, iaxis) * r;
			v = AngleAxisd(im, iaxis) * v;
		}
		if (f.norm() < 1e-8 && !phase1 || iter > 150)
		{
			if (iter <= 150)
			{
				cflag = true;
				cout << "UR3BP Phase2 Converged in " << iter << " Iterations\n";
			}
			else
			{
				cflag = false;
				cout << "UR3BP Failed to Converge with the Specified Parameters";
				rl_fail << y2[6], y2[7], y2[8];
				vl_fail << y2[9], y2[10], y2[11];
				rr_fail << y1[6], y1[7], y1[8];
				vr_fail << y1[9], y1[10], y1[11];
				print_format("rl:   ", rl_fail.norm());
				print_format("rr:   ", rr_fail.norm());
				print_format("phil: ", safe_acos(rl_fail.cross(vl_fail).norm() / (rl_fail.norm() * vl_fail.norm())) * 180 / pi);
				print_format("phir: ", safe_acos(rr_fail.cross(vr_fail).norm() / (rr_fail.norm() * vr_fail.norm())) * 180 / pi);
				print_format("il:   ", safe_acos(rl_fail.cross(vl_fail)(2) / rl_fail.cross(vl_fail).norm()) * 180 / pi);
				print_format("ir:   ", safe_acos(rr_fail.cross(vr_fail)(2) / rr_fail.cross(vr_fail).norm()) * 180 / pi);
				cout << "\n";
			}

			ofstream output_file;
			output_file.open("trajectory_solution.txt");
			output_file << setprecision(16) << vp_norm << '\n' << vaxy << '\n' << rp_norm << '\n' << raxy << '\n' << raxz << '\n' << tfl << '\n' << tfr << '\n';
			output_file.close();

			rm(0) = y2[0];
			rm(1) = y2[1];
			rm(2) = y2[2];
			vm(0) = y2[3];
			vm(1) = y2[4];
			vm(2) = y2[5];
			r(0) = y2[6];
			r(1) = y2[7];
			r(2) = y2[8];
			v(0) = y2[9];
			v(1) = y2[10];
			v(2) = y2[11];
			break;
		}
		if (f.norm() < 1e-8 && phase1)
		{
			cout << "UR3BP Phase1 Converged in " << iter << " Iterations\n";
			phase1 = false;
			iter = 0;
			inc_mul = 1.0;
		}
	}
	delete[] work1;
	delete[] y1;
	delete[] work2;
	delete[] y2;
	return;
}

void ur3bp_w(Vector3d r, Vector3d v, Vector3d rm, Vector3d vm, double tf1, double tf2)
{
	int i;
	int iflag;
	int iwork[5];
	int neqn = 48;
	int wi_yp = 100 + neqn * 4 - 1;
	int step_num = 2048;
	double abserr;
	double relerr;
	double t;
	double tout;
	double pi = 3.141592653589793;
	double *work;
	double *y;

	ofstream output_file1a;
	ofstream output_file1b;
	ofstream output_file2a;
	ofstream output_file2b;
	output_file1a.open("output1a.txt");
	output_file1b.open("output1b.txt");
	output_file2a.open("output2a.txt");
	output_file2b.open("output2b.txt");

	abserr = 1e-7;
	relerr = 0.0;
	iflag = 1;

	t = 0.0;
	work = new double[100 + 21 * neqn];
	y = new double[neqn];
	y_init_ur3bp(r, v, rm, vm, y);

	output_file1a << setprecision(16) << t << ',' << y[6] << ',' << y[7] << ',' << y[8] << ',' << y[9] << ',' << y[10] << ',' << y[11] << '\n';
	output_file2a << setprecision(16) << t << ',' << y[0] << ',' << y[1] << ',' << y[2] << ',' << y[3] << ',' << y[4] << ',' << y[5] << '\n';

	for (i = 1; i <= step_num; i++)
	{
		tout = (double)(i)* tf1 / (double)(step_num);
		ode(dur3bp, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);
		output_file1a << setprecision(16) << t << ',' << y[6] << ',' << y[7] << ',' << y[8] << ',' << y[9] << ',' << y[10] << ',' << y[11] << '\n';
		output_file2a << setprecision(16) << t << ',' << y[0] << ',' << y[1] << ',' << y[2] << ',' << y[3] << ',' << y[4] << ',' << y[5] << '\n';
	}

	output_file1b << setprecision(16) << t << ',' << y[6] << ',' << y[7] << ',' << y[8] << ',' << y[9] << ',' << y[10] << ',' << y[11] << '\n';
	output_file2b << setprecision(16) << t << ',' << y[0] << ',' << y[1] << ',' << y[2] << ',' << y[3] << ',' << y[4] << ',' << y[5] << '\n';

	for (i = 1; i <= step_num; i++)
	{
		tout = tf1 + (double)(i)* tf2 / (double)(step_num);
		ode(dur3bp, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);
		output_file1b << setprecision(16) << t << ',' << y[6] << ',' << y[7] << ',' << y[8] << ',' << y[9] << ',' << y[10] << ',' << y[11] << '\n';
		output_file2b << setprecision(16) << t << ',' << y[0] << ',' << y[1] << ',' << y[2] << ',' << y[3] << ',' << y[4] << ',' << y[5] << '\n';
	}

	output_file1a.close();
	output_file1b.close();
	output_file2a.close();
	output_file2b.close();

	delete[] work;
	delete[] y;
	return;
}

void ddd(double yl[], double ypl[], double yr[], double ypr[], Matrix6d &df, Vector6d &f, Vector6d &c, double rp_norm, double vp_norm, double im, double raxy, double raxz, double vaxy, Vector3d v, Vector3d iaxis, double inc)
{
	double* stml = yl + 6;
	double* stmr = yr + 6;

	Vector3d rl;
	Vector3d vl;
	rl << yl[0], yl[1], yl[2];
	vl << yl[3], yl[4], yl[5];

	Vector3d rr;
	Vector3d vr;
	rr << yr[0], yr[1], yr[2];
	vr << yr[3], yr[4], yr[5];

	Vector3d dvp_dvp_norm;
	dvp_dvp_norm = v / vp_norm;

	Vector3d dvp_dvaxy;
	Vector3d dvp_draxy;
	Vector3d dvp_draxz;
	Vector3d drp_draxy;
	Vector3d drp_draxz;

	dvp_dvaxy(0) = -sin(raxy) * sin(vaxy) + cos(raxy) * sin(raxz) * cos(vaxy);
	dvp_dvaxy(1) = cos(raxy) * sin(vaxy) + sin(raxy) * sin(raxz) * cos(vaxy);
	dvp_dvaxy(2) = -cos(raxz) * cos(vaxy);
	dvp_dvaxy *= vp_norm;

	dvp_draxy(0) = cos(raxy) * cos(vaxy) - sin(raxy) * sin(raxz) * sin(vaxy);
	dvp_draxy(1) = sin(raxy) * cos(vaxy) + cos(raxy) * sin(raxz) * sin(vaxy);
	dvp_draxy(2) = 0;
	dvp_draxy *= vp_norm;

	dvp_draxz(0) = cos(raxy) * cos(raxz) * sin(vaxy);
	dvp_draxz(1) = sin(raxy) * cos(raxz) * sin(vaxy);
	dvp_draxz(2) = sin(raxz) * sin(vaxy);
	dvp_draxz *= vp_norm;

	drp_draxy(0) = -rp_norm * sin(raxy) * cos(raxz);
	drp_draxy(1) = rp_norm * cos(raxy) * cos(raxz);
	drp_draxy(2) = 0;
	
	drp_draxz(0) = -rp_norm * cos(raxy) * sin(raxz);
	drp_draxz(1) = -rp_norm * sin(raxy) * sin(raxz);
	drp_draxz(2) = rp_norm * cos(raxz);

	dvp_dvaxy = AngleAxisd(im, iaxis) * dvp_dvaxy;
	dvp_draxy = AngleAxisd(im, iaxis) * dvp_draxy;
	dvp_draxz = AngleAxisd(im, iaxis) * dvp_draxz;
	drp_draxy = AngleAxisd(im, iaxis) * drp_draxy;
	drp_draxz = AngleAxisd(im, iaxis) * drp_draxz;

	Vector3d drlx_dr, drly_dr, drlz_dr;
	Vector3d dvlx_dr, dvly_dr, dvlz_dr;
	Vector3d drlx_dv, drly_dv, drlz_dv;
	Vector3d dvlx_dv, dvly_dv, dvlz_dv;
	drlx_dr << stml[0 * 6 + 0], stml[0 * 6 + 1], stml[0 * 6 + 2];
	drly_dr << stml[1 * 6 + 0], stml[1 * 6 + 1], stml[1 * 6 + 2];
	drlz_dr << stml[2 * 6 + 0], stml[2 * 6 + 1], stml[2 * 6 + 2];
	dvlx_dr << stml[3 * 6 + 0], stml[3 * 6 + 1], stml[3 * 6 + 2];
	dvly_dr << stml[4 * 6 + 0], stml[4 * 6 + 1], stml[4 * 6 + 2];
	dvlz_dr << stml[5 * 6 + 0], stml[5 * 6 + 1], stml[5 * 6 + 2];
	drlx_dv << stml[0 * 6 + 3], stml[0 * 6 + 4], stml[0 * 6 + 5];
	drly_dv << stml[1 * 6 + 3], stml[1 * 6 + 4], stml[1 * 6 + 5];
	drlz_dv << stml[2 * 6 + 3], stml[2 * 6 + 4], stml[2 * 6 + 5];
	dvlx_dv << stml[3 * 6 + 3], stml[3 * 6 + 4], stml[3 * 6 + 5];
	dvly_dv << stml[4 * 6 + 3], stml[4 * 6 + 4], stml[4 * 6 + 5];
	dvlz_dv << stml[5 * 6 + 3], stml[5 * 6 + 4], stml[5 * 6 + 5];

	Vector3d drrx_dr, drry_dr, drrz_dr;
	Vector3d dvrx_dr, dvry_dr, dvrz_dr;
	Vector3d drrx_dv, drry_dv, drrz_dv;
	Vector3d dvrx_dv, dvry_dv, dvrz_dv;
	drrx_dr << stmr[0 * 6 + 0], stmr[0 * 6 + 1], stmr[0 * 6 + 2];
	drry_dr << stmr[1 * 6 + 0], stmr[1 * 6 + 1], stmr[1 * 6 + 2];
	drrz_dr << stmr[2 * 6 + 0], stmr[2 * 6 + 1], stmr[2 * 6 + 2];
	dvrx_dr << stmr[3 * 6 + 0], stmr[3 * 6 + 1], stmr[3 * 6 + 2];
	dvry_dr << stmr[4 * 6 + 0], stmr[4 * 6 + 1], stmr[4 * 6 + 2];
	dvrz_dr << stmr[5 * 6 + 0], stmr[5 * 6 + 1], stmr[5 * 6 + 2];
	drrx_dv << stmr[0 * 6 + 3], stmr[0 * 6 + 4], stmr[0 * 6 + 5];
	drry_dv << stmr[1 * 6 + 3], stmr[1 * 6 + 4], stmr[1 * 6 + 5];
	drrz_dv << stmr[2 * 6 + 3], stmr[2 * 6 + 4], stmr[2 * 6 + 5];
	dvrx_dv << stmr[3 * 6 + 3], stmr[3 * 6 + 4], stmr[3 * 6 + 5];
	dvry_dv << stmr[4 * 6 + 3], stmr[4 * 6 + 4], stmr[4 * 6 + 5];
	dvrz_dv << stmr[5 * 6 + 3], stmr[5 * 6 + 4], stmr[5 * 6 + 5];

	Vector3d drl_dvp_norm;
	Vector3d drl_dvaxy;
	Vector3d drl_draxy;
	Vector3d drl_draxz;
	Vector3d drl_dt;

	drl_dvp_norm << drlx_dv.dot(dvp_dvp_norm), drly_dv.dot(dvp_dvp_norm), drlz_dv.dot(dvp_dvp_norm);
	drl_dvaxy << drlx_dv.dot(dvp_dvaxy), drly_dv.dot(dvp_dvaxy), drlz_dv.dot(dvp_dvaxy);
	drl_draxy << drlx_dr.dot(drp_draxy) + drlx_dv.dot(dvp_draxy), drly_dr.dot(drp_draxy) + drly_dv.dot(dvp_draxy), drlz_dr.dot(drp_draxy) + drlz_dv.dot(dvp_draxy);
	drl_draxz << drlx_dr.dot(drp_draxz) + drlx_dv.dot(dvp_draxz), drly_dr.dot(drp_draxz) + drly_dv.dot(dvp_draxz), drlz_dr.dot(drp_draxz) + drlz_dv.dot(dvp_draxz);
	drl_dt << -ypl[0], -ypl[1], -ypl[2];

	Vector3d dvl_dvp_norm;
	Vector3d dvl_dvaxy;
	Vector3d dvl_draxy;
	Vector3d dvl_draxz;
	Vector3d dvl_dt;

	dvl_dvp_norm << dvlx_dv.dot(dvp_dvp_norm), dvly_dv.dot(dvp_dvp_norm), dvlz_dv.dot(dvp_dvp_norm);
	dvl_dvaxy << dvlx_dv.dot(dvp_dvaxy), dvly_dv.dot(dvp_dvaxy), dvlz_dv.dot(dvp_dvaxy);
	dvl_draxy << dvlx_dr.dot(drp_draxy) + dvlx_dv.dot(dvp_draxy), dvly_dr.dot(drp_draxy) + dvly_dv.dot(dvp_draxy), dvlz_dr.dot(drp_draxy) + dvlz_dv.dot(dvp_draxy);
	dvl_draxz << dvlx_dr.dot(drp_draxz) + dvlx_dv.dot(dvp_draxz), dvly_dr.dot(drp_draxz) + dvly_dv.dot(dvp_draxz), dvlz_dr.dot(drp_draxz) + dvlz_dv.dot(dvp_draxz);
	dvl_dt << -ypl[3], -ypl[4], -ypl[5];

	Vector3d drr_dvp_norm;
	Vector3d drr_dvaxy;
	Vector3d drr_draxy;
	Vector3d drr_draxz;
	Vector3d drr_dt;

	drr_dvp_norm << drrx_dv.dot(dvp_dvp_norm), drry_dv.dot(dvp_dvp_norm), drrz_dv.dot(dvp_dvp_norm);
	drr_dvaxy << drrx_dv.dot(dvp_dvaxy), drry_dv.dot(dvp_dvaxy), drrz_dv.dot(dvp_dvaxy);
	drr_draxy << drrx_dr.dot(drp_draxy) + drrx_dv.dot(dvp_draxy), drry_dr.dot(drp_draxy) + drry_dv.dot(dvp_draxy), drrz_dr.dot(drp_draxy) + drrz_dv.dot(dvp_draxy);
	drr_draxz << drrx_dr.dot(drp_draxz) + drrx_dv.dot(dvp_draxz), drry_dr.dot(drp_draxz) + drry_dv.dot(dvp_draxz), drrz_dr.dot(drp_draxz) + drrz_dv.dot(dvp_draxz);
	drr_dt << ypr[0], ypr[1], ypr[2];

	Vector3d dvr_dvp_norm;
	Vector3d dvr_dvaxy;
	Vector3d dvr_draxy;
	Vector3d dvr_draxz;
	Vector3d dvr_dt;

	dvr_dvp_norm << dvrx_dv.dot(dvp_dvp_norm), dvry_dv.dot(dvp_dvp_norm), dvrz_dv.dot(dvp_dvp_norm);
	dvr_dvaxy << dvrx_dv.dot(dvp_dvaxy), dvry_dv.dot(dvp_dvaxy), dvrz_dv.dot(dvp_dvaxy);
	dvr_draxy << dvrx_dr.dot(drp_draxy) + dvrx_dv.dot(dvp_draxy), dvry_dr.dot(drp_draxy) + dvry_dv.dot(dvp_draxy), dvrz_dr.dot(drp_draxy) + dvrz_dv.dot(dvp_draxy);
	dvr_draxz << dvrx_dr.dot(drp_draxz) + dvrx_dv.dot(dvp_draxz), dvry_dr.dot(drp_draxz) + dvry_dv.dot(dvp_draxz), dvrz_dr.dot(drp_draxz) + dvrz_dv.dot(dvp_draxz);
	dvr_dt << ypr[3], ypr[4], ypr[5];

	df(0, 0) = ddhs(rl, drl_dvp_norm);
	df(0, 1) = ddhs(rl, drl_dvaxy);
	df(0, 2) = ddhs(rl, drl_draxy);
	df(0, 3) = ddhs(rl, drl_draxz);
	df(0, 4) = ddhs(rl, drl_dt);
	df(0, 5) = 0.0;

	df(1, 0) = ddhs(rr, drr_dvp_norm);
	df(1, 1) = ddhs(rr, drr_dvaxy);
	df(1, 2) = ddhs(rr, drr_draxy);
	df(1, 3) = ddhs(rr, drr_draxz);
	df(1, 4) = 0.0;
	df(1, 5) = ddhs(rr, drr_dt);

	df(2, 0) = ddphi(rl, vl, drl_dvp_norm, dvl_dvp_norm);
	df(2, 1) = ddphi(rl, vl, drl_dvaxy, dvl_dvaxy);
	df(2, 2) = ddphi(rl, vl, drl_draxy, dvl_draxy);
	df(2, 3) = ddphi(rl, vl, drl_draxz, dvl_draxz);
	df(2, 4) = ddphi(rl, vl, drl_dt, dvl_dt);
	df(2, 5) = 0.0;

	df(3, 0) = ddphi(rr, vr, drr_dvp_norm, dvr_dvp_norm);
	df(3, 1) = ddphi(rr, vr, drr_dvaxy, dvr_dvaxy);
	df(3, 2) = ddphi(rr, vr, drr_draxy, dvr_draxy);
	df(3, 3) = ddphi(rr, vr, drr_draxz, dvr_draxz);
	df(3, 4) = 0.0;
	df(3, 5) = ddphi(rr, vr, drr_dt, dvr_dt);

	df(4, 0) = ddi(rl, vl, drl_dvp_norm, dvl_dvp_norm);
	df(4, 1) = ddi(rl, vl, drl_dvaxy, dvl_dvaxy);
	df(4, 2) = ddi(rl, vl, drl_draxy, dvl_draxy);
	df(4, 3) = ddi(rl, vl, drl_draxz, dvl_draxz);
	df(4, 4) = ddi(rl, vl, drl_dt, dvl_dt);
	df(4, 5) = 0.0;

	df(5, 0) = ddi(rr, vr, drr_dvp_norm, dvr_dvp_norm);
	df(5, 1) = ddi(rr, vr, drr_dvaxy, dvr_dvaxy);
	df(5, 2) = ddi(rr, vr, drr_draxy, dvr_draxy);
	df(5, 3) = ddi(rr, vr, drr_draxz, dvr_draxz);
	df(5, 4) = 0.0;
	df(5, 5) = ddi(rr, vr, drr_dt, dvr_dt);

	f(0) = rl.squaredNorm() - c(0);
	f(1) = rr.squaredNorm() - c(1);
	f(2) = rl.cross(vl).squaredNorm() / (rl.squaredNorm() * vl.squaredNorm()) - c(2);
	f(3) = rr.cross(vr).squaredNorm() / (rr.squaredNorm() * vr.squaredNorm()) - c(3);
	f(4) = inc * (rl.cross(vl)(2) / rl.cross(vl).norm() - c(4));
	f(5) = inc * (rr.cross(vr)(2) / rr.cross(vr).norm() - c(5));
	return;
}

double ddhs(Vector3d r, Vector3d dr)
{
	return 2 * r.dot(dr);
}

double ddphi(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv)
{
	double r2;
	double v2;
	double dr2;
	double dv2;
	double r2v2;
	double dr2v2;
	double hm2;
	double dhm2;
	Vector3d h;
	Vector3d dh;

	r2 = r.squaredNorm();
	v2 = v.squaredNorm();
	dr2 = 2.0 * r.dot(dr);
	dv2 = 2.0 * v.dot(dv);
	r2v2 = r2 * v2;
	dr2v2 = dr2 * v2 + dv2 * r2;

	h = r.cross(v);
	dh = dr.cross(v) + r.cross(dv);
	hm2 = h.squaredNorm();
	dhm2 = 2.0 * h.dot(dh);

	return (dhm2 * r2v2 - dr2v2 * hm2) / pow(r2v2, 2.0);
}

double ddi(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv)
{
	double hz;
	double hm;
	double dhz;
	double dhm;
	Vector3d h;
	Vector3d dh;

	h = r.cross(v);
	dh = dr.cross(v) + r.cross(dv);

	hz = h(2);
	hm = h.norm();
	dhz = dh(2);
	dhm = h.dot(dh) / hm;

	return (dhz * hm - dhm * hz) / pow(hm, 2.0);
}