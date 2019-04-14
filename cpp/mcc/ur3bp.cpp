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

typedef Matrix<double, 4, 4> Matrix4d;
typedef Matrix<double, 4, 1> Vector4d;

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

void ur3bp_g(Vector3d &r, Vector3d &v, Vector3d &rm, Vector3d &vm, double tf, double c, bool &cflag)
{
	int iflag;
	int iter;
	int iwork[5];
	int neqn = 48;
	int wi_yp = 100 + neqn * 4 - 1;
	double abserr;
	double relerr;
	double t;
	double tout;
	double beta = 1.0;
	double fn;
	double fnold;
	double pi = 3.141592653589793;
	double *work;
	double *y;
	double *stm;
	double df;
	double f;
	double x;
	double xold;

	double rp0, rp1, drp_dv_norm;
	double dphi_dt, dtf, adtf, rdf;
	double *yf, *ypf;
	Vector3d rf, vf, drf_dt, dvf_dt;

	iter = 0;
	abserr = 1e-7;
	relerr = 0.0;

	y = new double[neqn];
	stm = y + 12;
	work = new double[100 + 21 * neqn];

	while (true)
	{
		iter++;
		t = 0.0;
		y_init_ur3bp(r, v, rm, vm, y);
		while (true)
		{
			tout = tf;

			iflag = 1;
			ode(dur3bp, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

			yf = y + 6;
			ypf = work + wi_yp + 6;

			rf << yf[0], yf[1], yf[2];
			vf << yf[3], yf[4], yf[5];
			drf_dt << ypf[0], ypf[1], ypf[2];
			dvf_dt << ypf[3], ypf[4], ypf[5];

			dphi_dt = ddphi(rf, vf, drf_dt, dvf_dt);

			f = rf.cross(vf).squaredNorm() / (rf.squaredNorm() * vf.squaredNorm()) - 1;
			dtf = f / dphi_dt;
			adtf = abs(dtf);
			rdf = rf.dot(vf);
			cout << "\n" << "1st" << "\n" << rdf * dtf;
			if (adtf > 0.01)
			{
				dtf = 0.01 * dtf / adtf;
			}
			if (dtf * rdf < 0)
			{
				dtf = 0.01 * rdf / abs(rdf);
			}

			

			tf -= dtf;
			if (abs(f) < 1e-8)
			{
				rp0 = rf.norm();
				break;
			}		
		}

		t = 0.0;
		y_init_ur3bp(r, v / v.norm() * (v.norm() + 1e-8), rm, vm, y);
		while (true)
		{
			tout = tf;

			iflag = 1;
			ode(dur3bp, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

			yf = y + 6;
			ypf = work + wi_yp + 6;

			rf << yf[0], yf[1], yf[2];
			vf << yf[3], yf[4], yf[5];
			drf_dt << ypf[0], ypf[1], ypf[2];
			dvf_dt << ypf[3], ypf[4], ypf[5];

			dphi_dt = ddphi(rf, vf, drf_dt, dvf_dt);

			f = rf.cross(vf).squaredNorm() / (rf.squaredNorm() * vf.squaredNorm()) - 1;
			dtf = f / dphi_dt;
			adtf = abs(dtf);
			rdf = rf.dot(vf);
			cout << "\n" << rdf * dtf;

			if (adtf > 0.01)
			{
				dtf = 0.01 * dtf / adtf;
			}
			if (dtf * rdf < 0)
			{
				dtf = 0.01 * rdf / abs(rdf);
			}

			tf -= dtf;
			if (abs(f) < 1e-8)
			{
				rp1 = rf.norm();
				break;
			}
		}
	
		df = (rp1 - rp0) / 1e-8;

		f = rp0 - c;
		fn = abs(f);

		if (iter == 1 && iflag != 2)
		{
			cout << "CR3BP_G - Fatal error - ODE returned IFLAG = " << iflag << " on first iteration\n";
			cflag = false;
			break;
		}
		else if (iter > 1 && (iflag != 2 || fn > fnold))
		{
			beta *= 0.25;
			v = v / v.norm() * ((1 - beta) * xold + beta * x);
		}
		else
		{
			beta = 1.0;
			xold = v.norm();
			fnold = abs(f);

			x = xold - f / df;
			v = v / v.norm() * x;
		}

		if (abs(f) < 1e-8)
		{
			cflag = true;
			cout << "UR3BP_G Converged in " << iter << " Iterations\n";
			ofstream output_file;
			output_file.open("mid_course_correction_maneuver.txt");
			output_file << setprecision(16) << v(0) << '\n' << v(1) << '\n' << v(2) << '\n';
			output_file.close();
			break;
		}
		if (iter > 100)
		{
			cout << "UR3BP_G Failed to Converge in 100 Iterations\n";
			cout << "Remaining error: " << setprecision(16) << f << "\n";
			ofstream output_file;
			output_file.open("mid_course_correction_maneuver.txt");
			output_file << setprecision(16) << v(0) << '\n' << v(1) << '\n' << v(2) << '\n';
			output_file.close();
			cflag = false;
			break;
		}
	}
	delete[] work;
	delete[] y;
	return;
}

double rp(Vector3d r, Vector3d v)
{
	double en, a, hs, e, rp;

	en = v.squaredNorm() / 2 - m2 / r.norm();
	a = -m2 / 2 / en;
	hs = r.cross(v).squaredNorm();
	e = pow(1 - hs / a / m2, 0.5);
	rp = a * (1 - e);

	return rp;
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