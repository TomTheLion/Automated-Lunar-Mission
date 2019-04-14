typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

bool file_to_vector(string filename, vector<double> &vector);
void print_format(string label, double value);
double safe_acos(double x);
void y_init_cr3bp(Vector3d r, Vector3d v, double y[]);
void y_init_ur3bp(Vector3d r, Vector3d v, Vector3d rm, Vector3d vm, double y[]);
void write_orbital_parameters(Vector3d r, Vector3d v, double dt, double mu, string file_name);
void dcr3bp(double t, double y[], double yp[]);
void dur3bp(double t, double y[], double yp[]);
void cr3bp_p(Vector3d &r, Vector3d &v, double &tf, double rp);
void ur3bp_i(Vector3d &r, Vector3d &v, Vector3d &rm, Vector3d &vm, Vector3d &iaxis, double em, double im, double lpe, double fm);void ur3bp(Vector3d &r, Vector3d &v, Vector3d &rm, Vector3d &vm, Vector3d iaxis, double im, double &tfl, double &tfr, Vector6d c, bool &cflag);
void ur3bp_w(Vector3d r, Vector3d v, Vector3d rm, Vector3d vm, double tf1, double tf2);
void ddd(double yl[], double ypl[], double yr[], double ypr[], Matrix6d &df, Vector6d &f, Vector6d &c, double rp_norm, double vp_norm, double im, double raxy, double raxz, double vaxy, Vector3d v, Vector3d iaxis, double inc);
double ddhs(Vector3d r, Vector3d dr);
double ddphi(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv);
double ddi(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv);