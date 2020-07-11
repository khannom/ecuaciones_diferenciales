#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
using namespace std;

double f1(double t, double x)
{
	return t+x;
}

double sol1(double t)
{
	return -t+2*exp(t)-1;
}

double f2(double t, double x)
{
	return (1/pow(t,4))-(2/t)*x;
}

double sol2(double t)
{
	return -1/pow(t,3)+2/pow(t,2);
}

double f3(double t, double x)
{
	return 0.25-0.1*x;
}

double sol3(double t)
{
	return -5*exp(-t/10)/2+(5.0/2);
}

void compareTable(const vector<double>& t,const vector<double>& x,const vector<double>& tt,const vector<double>& xx)
{
	system("clear");
	cout<<"valor t:\tvalor runge-kutta(4to orden)\t\tvalor real"<<endl<<endl;
	for (int i = 0; i < t.size(); ++i)
	{
		cout<<t[i]<<"\t\t\t"<<x[i]<<"\t\t\t\t"<<xx[i]<<endl;
	}
}

void comparePlot(const vector<double>& t,const vector<double>& x,const vector<double>& tt,const vector<double>& xx)
{
	ofstream rkdat("rk.dat");
	for (int i = 0; i < t.size(); ++i)
	{
		rkdat<<t[i]<<"\t"<<x[i]<<"\n";
	}
	rkdat.close();
	ofstream sol("sol.dat");
	for (int i = 0; i < tt.size(); ++i)
	{
		sol<<tt[i]<<"\t"<<xx[i]<<"\n";
	}
	sol.close();
	ofstream gnu("exe.gnu");
	string txt = "unset key\n";
	gnu<<"unset key\n";
	gnu<<"set xrange [0:11]\n";
	gnu<<"plot \"rk.dat\" lt 3 lc rgb \"blue\", \"sol.dat\" lt 7 w l";
	gnu.close();
	system("gnuplot -p exe.gnu");
}





void rungeKutta4(double a, double b, double x0, double h, double (*function)(double, double), double (*realSolution)(double))
{
	int n = (b-a)/h+1;
	vector<double> x;
	x.reserve(n);
	x.push_back(x0);
	vector<double> t;
	t.reserve(n);
	t.push_back(a);
	double k1,k2,k3,k4;
	for (double i = a+h, idx = 1; i < b; i += h, ++idx)
	{
		t.push_back(t[idx-1]+h);
		k1 = function(t[idx-1],x[idx-1]);
		k2 = function(t[idx-1]+0.5*h,x[idx-1]+0.5*k1*h);
		k3 = function(t[idx-1]+0.5*h,x[idx-1]+0.5*k2*h);
		k4 = function(t[idx-1]+h,x[idx-1]+k3*h);
		x.push_back(x[idx-1]+h*(k1+2*k2+2*k3+k4)/6);
	}
	vector<double> tt;
	vector<double> xx;
	double hh = 0.1;
	int nn = (b-a)/hh+1;
	tt.reserve(nn);
	tt.reserve(nn);
	for (double i = t.front(), idx = 0; i <= t.back(); i += hh, ++idx)
	{
		tt.push_back(i);
		xx.push_back(realSolution(tt[idx]));
	}
	compareTable(t,x,tt,xx);
	comparePlot(t,x,tt,xx);
}



int main()
{
	//Ejemplo 1
	//rungeKutta4(0,10,1,0.1,f1,sol1);
	//Ejemplo 2
	//rungeKutta4(1,11,1,0.1,f2,sol2);
	//Ejemplo 3
	rungeKutta4(0,10,0,0.1,f3,sol3);
	return 0;
}
