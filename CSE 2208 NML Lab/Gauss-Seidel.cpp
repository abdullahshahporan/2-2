#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

int main() {
	float a1, b1, c1, d1;
	float a2, b2, c2, d2;
	float a3, b3, c3, d3;
	cin >> a1 >> b1 >> c1 >> d1;
	cin >> a2 >> b2 >> c2 >> d2;
	cin >> a3 >> b3 >> c3 >> d3;
	float xp = 0, yp = 0, zp = 0, xn, yn, zn, ex = 0, ey = 0, ez = 0;
	int steps = 20;
	cout << "x\t\ty\t\tz\t\tEx\t\tEy\t\tEz" << endl;
	do {
		cout <<fixed<<setprecision(4)<< xp << "\t\t"<< yp << "\t\t" << zp << "\t\t" << abs(ex) << "\t\t" << abs(ey) << "\t\t" << abs(ez) << endl;
		xn = ((1 / a1) * (d1 - b1 * yp - c1 * zp));
		ex = xn - xp;
		xp = xn;
		yn = ((1 / b2) * (d2 - a2 * xp - c2 * zp));
		ey = yn - yp;
		yp = yn;
		zn = ((1 / c3) * (d3 - a3 * xp - b3 * yp));
		ez = zn - zp;
		zp = zn;
	} while (steps--);

	return 0;
}
