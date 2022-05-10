# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

# define GLFW_INCLUDE_GLU
# include <GLFW/glfw3.h>
# include <unistd.h>

# define PI 3.14159265358979323846

static const double R = 6371220.;
static const double g = 9.81;
static const double Gamma = 1e-7;
static const double Omega = 2 * PI / 86400;

static const double gaussTriangleXsi[3]    = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double gaussTriangleEta[3]    = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double gaussTriangleWeight[3] = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };
static const double gaussEdgeXsi[2]        = { 0.577350269189626,-0.577350269189626 };
static const double gaussEdgeWeight[2]     = { 1.000000000000000, 1.000000000000000 }; 

void    tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName);
void    tsunamiAnimate(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName);
double  tsunamiInitialConditionOkada(double x, double y);
void    tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub);
void    tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem);

