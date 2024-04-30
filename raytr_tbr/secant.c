#include <stdio.h>
#include <math.h>
  
double f(double x)
{
    return cos(x) - x*x*x;
}
  
double SecantMethod(double xn_1, double xn, double e, int m)
{
    int n;
    double d;
    for (n = 1; n <= m; n++)
    {
        d = (xn - xn_1) / (f(xn) - f(xn_1)) * f(xn);
        if (fabs(d) < e)
            return xn;
	printf("xn = %g, xn_1 = %g\n", xn, xn_1);
        xn_1 = xn;
        xn = xn - d;
    }
    return xn;
}
  
int main(void)
{
    printf("%0.15f\n", SecantMethod(0, 1, 5E-11, 100));

    printf("%i\n", 12345);
    
    return 0;

}
