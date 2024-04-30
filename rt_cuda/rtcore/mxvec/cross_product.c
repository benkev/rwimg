//
// Vector cross product
// c = a x b
//
inline void cross_product(double a[3], double b[3],double c[3]) {
  c[0] =  a[1]*b[2] - a[2]*b[1];
  c[1] = -a[0]*b[2] + a[2]*b[0];
  c[2] =  a[0]*b[1] - a[1]*b[0];
}
