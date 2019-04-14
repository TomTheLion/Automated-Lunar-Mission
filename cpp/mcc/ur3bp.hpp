typedef Matrix<double, 4, 4> Matrix4d;
typedef Matrix<double, 4, 1> Vector4d;

bool file_to_vector(string filename, vector<double> &vector);
void y_init_ur3bp(Vector3d r, Vector3d v, Vector3d rm, Vector3d vm, double y[]);
void dur3bp(double t, double y[], double yp[]);
void ur3bp_g(Vector3d &r, Vector3d &v, Vector3d &rm, Vector3d &vm, double tf, double c, bool &cflag);
double rp(Vector3d r, Vector3d v);
double dt(Vector3d r, Vector3d v);
void ur3bp_w(Vector3d r, Vector3d v, Vector3d rm, Vector3d vm, double tf1, double tf2);
double ddphi(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv);