#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>

using namespace std;

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592;

// check if any particles are touching
bool anytouch(int N, double * x, double * y, double * R_eff, double Lx, double Ly)
{
    double dx, dy, Dnm, dnm2;    
	for(int nn = 0; nn < N-1; nn++)
	{
		for(int mm = nn + 1; mm < N; mm++)
		{
			//vector pointing from particle i to j
            dx = x[mm] - x[nn];
            dx = dx - round(dx/Lx)*Lx;
            dy = y[mm] - y[nn];
            dy = dy - round(dy/Ly)*Ly;
			//sum of 2 radii
			Dnm = R_eff[nn] + R_eff[mm];
            dnm2 = dx*dx + dy*dy;
			//they touch if the distance between them is less than Dnm square
			if(dnm2 < Dnm * Dnm) return true;
		}
	}
	return false;
}

void space_equally_circle(int n, double * xval, double * yval, double * lengths)
{
    // This function finds the equally spaced point on a unit circle of radius 1
    double circumference = 2.0 * pi;
    double dtheta = circumference / n;
    for (int i=0;i<n;i++)
    {
        xval[i] = cos(dtheta * i);
        yval[i] = sin(dtheta * i);
    }
    double diff = sqrt((xval[1] - xval[0]) * (xval[1] - xval[0]) + (yval[0] - yval[1]) * (yval[0] - yval[1]));
    for (int i=0;i<n;i++)
    {
        lengths[i] = diff;
    }
}

void printBasicNmers(int n, double * xval, double * yval, double * lengths)
{
    printf("Basic nmers:\n");
    printf("i, xval, yval, lengths\n");
    for (int i=0;i<n;i++)
    {
        printf("%d, %1.6f, %1.6f, %1.6f\n",i,xval[i],yval[i],lengths[i]);
    }
}

void printNmers(int n, int Nc, double * x_shape, double * y_shape, double * r_shape, double * th_shape)
{
    printf("x_shape, y_shape, r_shape, th_shape\n");
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            printf("%1.6f, %1.6f, %1.6f, %1.6f\n",x_shape[j*n+i],y_shape[j*n+i],r_shape[j*n+i],th_shape[j*n+i]);
        }
        printf("\n");
    }

}

void printParticles(int Nc, double * m, double * I)
{
    printf("m, I\n");
    for (int j=0;j<Nc;j++)
    {
        printf("%1.6f, %1.6f\n",m[j],I[j]);
    }
}

void nmersBuild(int n, int Nc, double G, double * xval, double * yval, double * x_shape, double * y_shape, double * r_shape, double * th_shape)
{
    // fill the shape arrays
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            x_shape[j*n+i] = xval[i];
            y_shape[j*n+i] = yval[i];
        }
    }
    // scale with particle size
    for (int j=Nc/2;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            x_shape[j*n+i] = x_shape[j*n+i] * G;
            y_shape[j*n+i] = y_shape[j*n+i] * G;
        }
    }
    // fill r and th array
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            r_shape[j*n+i] = sqrt(x_shape[j*n+i] * x_shape[j*n+i] + y_shape[j*n+i] * y_shape[j*n+i]);
            th_shape[j*n+i] = atan2(y_shape[j*n+i],x_shape[j*n+i]);
        }
    }
    // recalculate x,y array
    for (int j=0;j<Nc;j++)
    {
        for (int i=0;i<n;i++)
        {
            x_shape[j*n+i] = r_shape[j*n+i] * cos(th_shape[j*n+i]);
            y_shape[j*n+i] = r_shape[j*n+i] * sin(th_shape[j*n+i]);
        }
    }

}

int inpolygon(int nvert, double * vertx, double * verty, double testx, double testy)
{
    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) 
    {
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
        (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
        c = !c;
    }
    return c;
}

void inpart(int n ,long int MCpoints, int * in, double * x, double * y, double * xp, double * yp, double R)
{
    double R2 = R * R;
    bool tempbool;
    // inpolygon
    for (long int i=0;i<MCpoints;i++)
    {
        in[i] = inpolygon(n,xp,yp,x[i],y[i]);
    }
    for (int j=0;j<n;j++)
    {
        for (long int i=0;i<MCpoints;i++)
        {
            tempbool = ((x[i] - xp[j]) * (x[i] - xp[j]) + (y[i] - yp[j]) * (y[i] - yp[j])) < R2;
            in[i] = in[i] | tempbool;
        }
    }
}

void nmersProperty(int n, int Nc, int long MCpoints, double * xlist, double * ylist, double * x, double * y, int * in, double * m, double * I, double * x_shape, double * y_shape, double Rl, double Rs)
{
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);

    // long int MCpoints = 1e5;
    // double * xlist = new double[n];
    // double * ylist = new double[n];
    // double * Rlist = new double[n];
    // double * x = new double[MCpoints];
    // double * y = new double[MCpoints];
    // int * in = new int[MCpoints];
    double Lx_mc, Ly_mc, xmin, ymin, Atemp, Itemp, temp, count, sum, xcm, ycm, xsum, ysum, temp_max, temp_min;

    // first deal with small particles
    for (int i=0;i<n;i++)
    {
        xlist[i] = x_shape[i];
        ylist[i] = y_shape[i];
    }
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((xlist[i] + Rs) > temp_max) temp_max = (xlist[i] + Rs);
        if ((xlist[i] - Rs) < temp_min) temp_min = (xlist[i] - Rs);
    }
    Lx_mc = temp_max - temp_min;
    xmin = temp_min;
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((ylist[i] + Rs) > temp_max) temp_max = (ylist[i] + Rs);
        if ((ylist[i] - Rs) < temp_min) temp_min = (ylist[i] - Rs);
    }
    Ly_mc = temp_max - temp_min;
    ymin = temp_min;
    for (long int i=0;i<MCpoints;i++)
    {
        x[i] = Lx_mc * distribution(generator) + xmin;
        y[i] = Ly_mc * distribution(generator) + ymin;
    }
    inpart(n, MCpoints, in, x, y, xlist, ylist, Rs);
    // Area
    sum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) sum += 1.0;
    Atemp = (sum / double(MCpoints)) * Lx_mc * Ly_mc;
    // CM
    xsum = 0;
    ysum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) xsum += x[i];
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) ysum += y[i];
    xcm = xsum / sum;
    ycm = ysum / sum;
    // I
    temp = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) temp += ((x[i] - xcm)*(x[i] - xcm) + (y[i] - ycm)*(y[i] - ycm));
    Itemp = Atemp * temp / sum;
    for (int i=0;i<Nc/2;i++)
    {
        m[i] = Atemp;
        I[i] = Itemp;
    }
    // done with small particles

    // then deal with large particles
    for (int i=0;i<n;i++)
    {
        xlist[i] = x_shape[(Nc/2)*n+i];
        ylist[i] = y_shape[(Nc/2)*n+i];
    }
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((xlist[i] + Rl) > temp_max) temp_max = (xlist[i] + Rl);
        if ((xlist[i] - Rl) < temp_min) temp_min = (xlist[i] - Rl);
    }
    Lx_mc = temp_max - temp_min;
    xmin = temp_min;
    temp_max = -INFINITY;
    temp_min = INFINITY;
    for (int i=0;i<n;i++)
    {
        if ((ylist[i] + Rl) > temp_max) temp_max = (ylist[i] + Rl);
        if ((ylist[i] - Rl) < temp_min) temp_min = (ylist[i] - Rl);
    }
    Ly_mc = temp_max - temp_min;
    ymin = temp_min;
    for (long int i=0;i<MCpoints;i++)
    {
        x[i] = Lx_mc * distribution(generator) + xmin;
        y[i] = Ly_mc * distribution(generator) + ymin;
    }
    inpart(n, MCpoints, in, x, y, xlist, ylist, Rl);
    // Area
    sum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) sum += 1.0;
    Atemp = (sum / double(MCpoints)) * Lx_mc * Ly_mc;
    // CM
    xsum = 0;
    ysum = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) xsum += x[i];
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) ysum += y[i];
    xcm = xsum / sum;
    ycm = ysum / sum;
    // I
    temp = 0;
    for (int long i=0;i<MCpoints;i++) if (in[i] == 1) temp += ((x[i] - xcm)*(x[i] - xcm) + (y[i] - ycm)*(y[i] - ycm));
    Itemp = Atemp * temp / sum;
    for (int i=Nc/2;i<Nc;i++)
    {
        m[i] = Atemp;
        I[i] = Itemp;
    }

    // delete[] xlist, ylist, x, y, in;
}

double average(int N, double * numList)
{
    double average;
    double temp = 0.0;
    for (int i=1;i<N;i++)
    {
        temp += numList[i];
    }
    average = temp / N;
    return average;
}

double arrminf(int N, double * arr)
{
    double temp_min = INFINITY;
    for (int i=0;i<N;i++) if (arr[i] < temp_min) temp_min = arr[i];
    return temp_min;
}

double arrmaxf(int N, double * arr)
{
    double temp_max = -INFINITY;
    for (int i=0;i<N;i++) if (arr[i] > temp_max) temp_max = arr[i];
    return temp_max;
}

double arrsumf(int N, double * arr)
{
    double temp_sum = 0;
    for (int i=0;i<N;i++) temp_sum += arr[i];
    return temp_sum;
}

// Set CM velocity to 0
void CMVelocityZeroing(int N, double * m, double * vx)
{
    double M = 0;
    double P = 0;
    for (int i = 0; i < N; i++)
    {
        M += m[i];
        P += (m[i] * vx[i]);
    }
    double CM_v = P / M;
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] - CM_v;
    }
}

// Velocity Verlet position integration
void VV_pos_integration(int N, double * x, double * vx, double * ax_old, double dt, double half_dt2)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = x[i] + vx[i] * dt + ax_old[i] * half_dt2;
    }
}

// Velocity Verlet velocity integration
void VV_vel_integration(int N, double * vx, double * ax, double * ax_old, double dt)
{
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] + (ax[i] + ax_old[i]) * dt / 2.0;
    }
}

void getAcceleration(int N, double * Fx, double * ax, double * m)
{
    for (int i = 0; i < N; i++) ax[i] = Fx[i] / m[i];
}

void copyAcceleration(int N, double * ax, double * ax_old)
{
    for (int i = 0; i < N; i++) ax_old[i] = ax[i];
}

void correctionFIRE(int Nc, double * vx, double * Fx, double a)
{
    double normF, normv;
    normF = 0;
    normv = 0;
    for (int i=0;i<Nc;i++)
    {
        normF += Fx[i] * Fx[i];
        normv += vx[i] * vx[i];
    }
    normF = sqrt(normF);
    normv = sqrt(normv);
    for (int i=0;i<Nc;i++)
    {
        vx[i] = (1.0 - a) * vx[i] + a * (Fx[i]/normF) * normv;
    }
}

// Set noncontacting particles to zero velocity
void NonContactZeroing(int N, double * vx, int * Cn)
{
    for (int i = 0; i < N; i++)
    {
        if (Cn[i] == 0)
        {
            vx[i] = 0.0;
        }
    }
}

// Calculate Ek
// double getEk(int N, double * m, double * vx)
// {
//     double K = 0;
//     for (int i = 0; i < N; i++)
//     {
//         K += m[i] * vx[i] * vx[i];
//     }
//     K = K / 2.0;
//     return K;
// }

// Calculate max Ek
double getEk(int N, double * m, double * I, double * vx, double * vy, double * w)
{
    double * K = new double[N];
    for (int i = 0; i < N; i++) K[i] = 0.0;
    for (int i = 0; i < N; i++)
    {
        K[i] += (m[i] * vx[i] * vx[i] + m[i] * vy[i] * vy[i] + I[i] * w[i] * w[i])/2;
    }
    double Kmax = 0.0;
    for (int i = 0; i < N; i++) if (K[i] > Kmax) Kmax = K[i];

    return Kmax;
}

double bumpy_2D_shrjam_getU(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape,
double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double U, dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn;
    U = 0;
    int nn, mm, nnn, mmm;
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm=sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    U += K * ((DDnm - ddnm) * (DDnm - ddnm)/2);
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    U = U / (K * Nc);
    return U;
}

int bumpy_2D_shrjam_getC(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape,
double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn;
    int C = 0;
    int nn, mm, nnn, mmm;
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm=sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    C += 1;
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    return C;
}

double bumpy_2D_stress(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape,
double * Fx, double * Fy, double * stress, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double P, dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn, F;
    int nn, mm, nnn, mmm;
    // Zeroing Fx, Fy and T
    for (int i=0; i < Nc; i++)
    {        
        Fx[i] = 0.0;
        Fy[i] = 0.0;
    }
    for (int i=0;i<4;i++) stress[i] = 0.0;
    // Overlap detection
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm = sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    F = -K*(DDnm/ddnm-1);
                                    Fx[nn] += F * dxx;
                                    Fx[mm] -= F * dxx;
                                    Fy[nn] += F * dyy;
                                    Fy[mm] -= F * dyy;
                                    stress[0] -= F * dxx * dx;
                                    stress[1] -= 0.5 * F * (dxx*dy + dyy*dx);
                                    stress[2] -= 0.5 * F * (dyy*dx + dxx*dy);
                                    stress[3] -= F * dyy * dy;
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    for (int i=0;i<4;i++) stress[i] /= (Lx * Ly);
    P = (stress[0] + stress[3]) / 2.0;
    return P;
}

double bumpy_2D_Force(int Nc, int N, int n, double * x, double * y, double * th, double * r_shape, double * th_shape,
double * Fx, double * Fy, double * T, int * Cn, double K, double * R_eff, double * Dn, double Lx, double Ly, double gam)
{
    double U, dx, dy, dxx, dyy, im, Dnm, dnm, DDnm, ddnm, rymmm, rynnn, rxmmm, rxnnn, F;
    U = 0;
    int nn, mm, nnn, mmm;
    // Zeroing Fx, Fy and T
    for (int i = 0; i < Nc; i++)
    {        
        Fx[i] = 0.0;
        Fy[i] = 0.0;
        T[i] = 0.0;
        Cn[i] = 0;
    }
    // Overlap detection
    for (nn=0;nn<Nc;nn++)
    {
        for (mm=nn+1;mm<Nc;mm++)
        {
            dy = y[mm] - y[nn];
            im = round(dy/Ly);
            dy = dy - im * Ly;
            Dnm = R_eff[nn] + R_eff[mm];
            if (abs(dy) < Dnm)
            {
                dx = x[mm] - x[nn];
                dx = dx - round(dx/Lx-im*gam)*Lx - im*gam*Lx;
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    for (nnn=0;nnn<n;nnn++)
                    {
                        for (mmm=0;mmm<n;mmm++)
                        {
                            rymmm = r_shape[mm*n+mmm] * sin(th[mm]+th_shape[mm*n+mmm]);
                            rynnn = r_shape[nn*n+nnn] * sin(th[nn]+th_shape[nn*n+nnn]);
                            dyy = (y[mm] + rymmm) - (y[nn]+rynnn);
                            im = round(dyy/Ly);
                            dyy = dyy-im*Ly;
                            DDnm = (Dn[nn]+Dn[mm])/2;
                            if (abs(dyy) < Dnm)
                            {
                                rxmmm = r_shape[mm*n+mmm] * cos(th[mm]+th_shape[mm*n+mmm]);
                                rxnnn = r_shape[nn*n+nnn] * cos(th[nn]+th_shape[nn*n+nnn]);
                                dxx = (x[mm] + rxmmm) - (x[nn]+rxnnn);
                                dxx = dxx-round(dxx/Lx-im*gam)*Lx-im*gam*Lx;
                                ddnm = sqrt(dxx * dxx + dyy * dyy);
                                if (ddnm < DDnm)
                                {
                                    Cn[nn] += 1;
                                    Cn[mm] += 1;
                                    F = -K*(DDnm/ddnm-1);
                                    Fx[nn] += F * dxx;
                                    Fx[mm] -= F * dxx;
                                    Fy[nn] += F * dyy;
                                    Fy[mm] -= F * dyy;
                                    T[nn] += F * (rxnnn * dyy - rynnn * dxx);
                                    T[mm] -= F * (rxmmm * dyy - rymmm * dxx);
                                    U += K * ((DDnm - ddnm) * (DDnm - ddnm)/2);
                                }
                            }

                        }
                    }
                }

            }
        }
    }
    U = U / (K * Nc);
    return U;
}

double bumpy_2D_comp_VV(int Nc, int N, int n, double gam, double * Dn, double * r_shape, double * th_shape, double * m, double * I, double * R_eff, double * x, double * y, 
double * th, double * Fx, double * Fy, double * T, int * Cn, double * vx, double * vy, double * w, double * ax, double * ay, double * alph, double * ax_old, double * ay_old, double * alph_old, double dt, long int Nt, double K, double Lx, double Ly, double Utol, double dphi)
{
    // double * ax = new double[Nc];
    // double * ay = new double[Nc];
    // double * alph = new double[Nc];

    double phitot, rsc, U, im, dt2, half_dt2, Ek_trans, Ek_rot, Ek_tot;
    int i, j;
    long int nt;
    // zeroing initial acceleration
    for (i=0;i<Nc;i++)
    {
        ax[i] = 0.0;
        ay[i] = 0.0;
        alph[i] = 0.0;
    }
    dt2 = dt * dt;
    half_dt2 = dt2 / 2.0;
    phitot = arrsumf(Nc, m) / (Lx * Ly);
    rsc = sqrt(1.0 + dphi / phitot);
    // Grow
    for (i=0;i<Nc;i++)
    {
        Dn[i] = Dn[i] * rsc;
        R_eff[i] = R_eff[i] * rsc;
        m[i] = m[i] * rsc * rsc;
        I[i] = I[i] * rsc * rsc * rsc * rsc;
    }
    for (i=0;i<N;i++) r_shape[i] = r_shape[i] * rsc;
    phitot = arrsumf(Nc, m) / (Lx * Ly);

    U = bumpy_2D_shrjam_getU(Nc, N, n, x, y, th, r_shape, th_shape, K, R_eff, Dn, Lx, Ly, gam);
    // printf("U=%e.\n",U);
    // FIRE COEFF. //
    int nfiremin, cut;
    double finc, fdec, astart, a, fa, dtmax, P;
    // bool fire = true;
    bool fire = false;
    nfiremin = 50;
    finc = 1.1;
    fdec = 0.5;
    astart = 0.1;
    a = 0.1;
    fa = 0.99;
    dtmax = 10.0 * dt;
    cut = 1;

    double B = 0.1;
    double Kmax;

    nt = 0;
    while (U > 0 && nt < Nt)
    {
        // MOD
        for (int i=0;i<Nc;i++)
        {
            im = floor(y[i] / Ly);
            x[i] = fmod(x[i]-im*gam*Lx,Lx);
            if (x[i] < 0) x[i] += Lx;
            y[i] = fmod(y[i],Ly);
            if (y[i] < 0) y[i] += Ly;
            th[i] = fmod(th[i],2*pi);
            if (th[i] < 0) th[i] += 2*pi;
        }
        // zeroing
        CMVelocityZeroing(Nc, m, vx);
        CMVelocityZeroing(Nc, m, vy);

        // Ek
        // Ek_trans = getEk(Nc, m, vx);
        // Ek_trans = Ek_trans + getEk(Nc, m, vy);
        // Ek_rot = getEk(Nc, I, w);
        // Ek_tot = Ek_trans + Ek_rot;
        // Ek_tot = Ek_tot / Nc;
        Kmax = getEk(Nc, m, I, vx, vy, w);

        // first step in Verlet integration
        VV_pos_integration(Nc, x, vx, ax_old, dt, half_dt2);
        VV_pos_integration(Nc, y, vy, ay_old, dt, half_dt2);
        VV_pos_integration(Nc, th, w, alph_old, dt, half_dt2);

        U = bumpy_2D_Force(Nc,N,n,x,y,th,r_shape,th_shape,Fx,Fy,T,Cn,K,R_eff,Dn,Lx,Ly,gam);

        // Damping
        if (~fire)
        {
            for (int i = 0; i < Nc; i++) Fx[i] = Fx[i] - B * 2 * sqrt(m[i]/pi) * vx[i];
            for (int i = 0; i < Nc; i++) Fy[i] = Fy[i] - B * 2 * sqrt(m[i]/pi) * vy[i];
            for (int i = 0; i < Nc; i++) T[i] = T[i] - B * 2 * sqrt(m[i]/pi) * w[i]; 
        }

        getAcceleration(Nc, Fx, ax, m);
        getAcceleration(Nc, Fy, ay, m);
        getAcceleration(Nc, T, alph, I);

        // FIRE
        ///*
        P = 0.0;
        if (fire)
        {
            for (int i=0;i<Nc;i++)
            {
                P += (vx[i] * Fx[i] + vy[i] * Fy[i] + w[i] * T[i]);
            }
            correctionFIRE(Nc, vx, Fx, a);
            correctionFIRE(Nc, vy, Fy, a);
            correctionFIRE(Nc, w, T, a);
            if (P < 0)
            {
                for (int i=0;i<Nc;i++)
                {
                    vx[i] = 0.0;
                    vy[i] = 0.0;
                    w[i] = 0.0;
                }
                cut = nt;
                dt = dt * fdec;
                a = astart;
            }
            else if (P >= 0 && nt - cut > nfiremin)
            {
                if (dt * finc < dtmax) dt = dt * finc;
                else dt = dtmax;
                a = a * fa;
            }
        }

        VV_vel_integration(Nc, vx, ax, ax_old, dt);
        VV_vel_integration(Nc, vy, ay, ay_old, dt);
        VV_vel_integration(Nc, w, alph, alph_old, dt);

        copyAcceleration(Nc, ax, ax_old);
        copyAcceleration(Nc, ay, ay_old);
        copyAcceleration(Nc, alph, alph_old);

        if (~fire)
        {
            NonContactZeroing(Nc, vx, Cn);
            NonContactZeroing(Nc, vy, Cn);
            NonContactZeroing(Nc, w, Cn); 
        }
        
        if (nt > 0 && Kmax < 1e-24) 
        {
            // printf("Break by Ek=%e at %ld.\n",Kmax,nt);
            // print("\n");
            break;
        }

        nt += 1;
    }

    // printf("U=%e.\n",U);
    // delete[] ax, ay, alph;
    return U;
}

int main(int argc, char **argv)
{
    // TO COMPILE (MAC or LINUX?)
    // clang++ -std=c++14 -stdlib=libc++ Bumpyjam_v0.cpp 
    // TO COMPILE (WINDOWS w/ VS_community 2017+)
    // cl -EHsc -std:c++14 Bumpyjam_v0.cpp 
    // TO COMPILE (WINDOWS w/ clang - not perfect, only works for -m32)
    // clang++ -std=c++14 -m32 -stdlib=libc++ Bumpyjam_v0.cpp

    // TO RUN (in the pseudo-sense...)
    // ./compiled_program Nc n /mu phi_target seed
    // you can just do the following:
    // ./compiled_program 6 16 1 1 1

    //Save output file
    string space = "_";
	string fileext = ".out";

    string debugfile = "";
    string debugfileext = ".xy";

    // string filename = "E:/Dropbox/Yale/C++/Bumpy/output_"; // Windows
    string filename = "/Users/philipwang/Dropbox/Yale/C++/Bumpy/output_"; // Mac
    filename.append(argv[1]);
	filename.append(space);
	filename.append(argv[2]);
	filename.append(space);
	filename.append(argv[3]);
	filename.append(space);
	filename.append(argv[4]);
    filename.append(space);
	filename.append(argv[5]);
    debugfile.append(filename);
    filename.append(fileext);
    debugfile.append(debugfileext);
    char *filechar = new char[filename.length() + 1];
	strcpy(filechar, filename.c_str());
	FILE *out = fopen(filechar, "w+");
    printf("File char is: %s\n", filechar);

    char *debugfilechar = new char[debugfile.length() + 1];
    strcpy(debugfilechar, debugfile.c_str());
    FILE *debug = fopen(debugfilechar, "w+");

    int Nc;// = 6;
    int n;// = 16; // n-mers
    double mu;// = 1.0;
    double phi_target;// = 1.00;
    int seed;// = 1;

    Nc = atoi(argv[1]);
    n = atoi(argv[2]);
    mu = atof(argv[3]);
    phi_target = atof(argv[4]);
    seed = atoi(argv[5]);

    double temp, U, dphi, dG, Utol, tol, gam, phitot;
    int count, count_max, C;
    int N = n * Nc;

    mt19937 generator(seed);

    // create particles
    int Nsmall = Nc / 2;
    int Nbig = Nc / 2;

    // create basic nmers
    double * xval = new double[n];
    double * yval = new double[n];
    double * lengths = new double[n];
    space_equally_circle(n, xval, yval, lengths);
    // printBasicNmers(n, xval, yval, lengths);

    double l = average(n, lengths);
    double D = (l/2) * sqrt(1 + mu * mu) / mu;
    for (int i=0;i<n;i++)
    {
        xval[i] = xval[i] / D;
        yval[i] = yval[i] / D;
        lengths[i] = lengths[i] / D;
    }
    double rad = 1.0 / D;
    D = 1.0;
    double E0 = 1.0;
    double K = 2.0 * E0 * D * D;
    double G = 1.4;

    // set spatial quantities
    double Lx = round(5 * G * sqrt(Nc) * (D + rad));
    double Ly = round(5 * G * sqrt(Nc) * (D + rad));

    int long Nt = 1e6;
    double N_per_coll = 100.0;
    double dt = 2.0 * pi * sqrt(1.0 / K) / N_per_coll;

    // set intermediate MD arrays
    double * Dn0 = new double[Nc];
    double * R_eff0 = new double[Nc];
    double * m0 = new double[Nc];
    double * I0 = new double[Nc]; // moment of inertia of disks along center = 1/2 * m * r * r    
    double * x0 = new double[Nc];
    double * y0 = new double[Nc];
    double * th0 = new double[Nc];
    double * vx0 = new double[Nc];
    double * vy0 = new double[Nc];
    double * w0 = new double[Nc];
    double * ax_old0 = new double[Nc];
    double * ay_old0 = new double[Nc];
    double * alph_old0 = new double[Nc];
    // build nmers
    double * x_shape0 = new double [N];
    double * y_shape0 = new double [N];
    double * r_shape0 = new double [N];
    double * th_shape0 = new double [N];

    // set MD arrays
    double * Dn = new double[Nc];
    double * R_eff = new double[Nc];
    double * m = new double[Nc];
    double * I = new double[Nc]; // moment of inertia of disks along center = 1/2 * m * r * r    
    double * x = new double[Nc];
    double * y = new double[Nc];
    double * th = new double[Nc];
    double * vx = new double[Nc];
    double * vy = new double[Nc];
    double * w = new double[Nc];
    double * ax_old = new double[Nc];
    double * ay_old = new double[Nc];
    double * alph_old = new double[Nc];

    // local MD arrays
    double * ax = new double[Nc];
    double * ay = new double[Nc];
    double * alph = new double[Nc];

    // build nmers
    double * x_shape = new double [N];
    double * y_shape = new double [N];
    double * r_shape = new double [N];
    double * th_shape = new double [N];
    // set nmer local variables
    long int MCpoints = 1e5;
    double * xlist = new double[n];
    double * ylist = new double[n];
    // double * Rlist = new double[n];
    double * x_mc = new double[MCpoints];
    double * y_mc = new double[MCpoints];
    int * in_mc = new int[MCpoints];

    // stress tensor
    double * stress = new double [4];
    double P, Ptol;

    for(int i = 0; i < Nc / 2; i++) R_eff[i] = D / 1.0;
	for(int i = Nc / 2; i < Nc; i++) R_eff[i] = G / 1.0;
    for(int i = 0; i < Nc / 2; i++) Dn[i] = D; // small
	for(int i = Nc / 2; i < Nc; i++) Dn[i] = G; // big

    uniform_real_distribution<double> distribution(0.0, Lx); // Generate number between 0, Lx    
    // seed random initial particle positions
	do
	{
		for(int i = 0; i < Nc; i++) x[i] = distribution(generator);
        for(int i = 0; i < Nc; i++) y[i] = distribution(generator);
	}while(anytouch(Nc, x, y, R_eff, Lx, Ly));//re-seed until no particles touch

    printf("Particles ready!\n");

    nmersBuild(n, Nc, G, xval, yval, x_shape, y_shape, r_shape, th_shape);
    nmersProperty(n, Nc, MCpoints, xlist, ylist, x_mc, y_mc, in_mc, m, I, x_shape, y_shape, G/2.0, D/2.0);
    // printNmers(n, Nc, x_shape, y_shape, r_shape, th_shape);
    // printf("%e\n",-INFINITY);
    // printParticles(Nc, m, I);
    temp = 0;
    for(int i = 0; i < Nc; i++) temp += m[i];
    phitot = temp / (Lx * Ly);
    // printf("Lx=%e,Lx*Ly=%e,phi_ini=%e.\n",Lx,Lx*Ly,phitot);
    
    for (int i=0;i<Nc;i++)
    {
        temp = -INFINITY;
        for (int j=0;j<n;j++) if (r_shape[i*n+j] > temp) temp = r_shape[i*n+j]; // get max
        R_eff[i] = temp + Dn[i] / 2;
    }

    // update dt
    double m_mean = average(Nc, m);
    dt = 2.0 * pi * sqrt(m_mean / K) / N_per_coll;

    // First print out a file to plot the unjam states

    double * Fx = new double[Nc];
    double * Fy = new double[Nc];
    double * T = new double[Nc];
    int * Cn = new int[Nc];

    count_max = 1e6;
    gam = 0;
    Utol = 1e-14;
    Ptol = 1e-7;
    tol = 1e-7;
    dphi = 1e-3;
    count = 0;
    C = 0;
    U = 0;
    P = 0;

    double Kmax = 0;

    // Print initial states
    fprintf(debug, "%d\n", Nc);
    fprintf(debug, "%d\n", n);
    fprintf(debug, "%d\n", N);
    fprintf(debug, "%3.4f\n", mu);
    fprintf(debug, "%03.16f\n", Lx);
    fprintf(debug, "%03.16f\n", Ly);
    fprintf(debug, "%03.20f\n", gam);
    for (int j = 0; j < Nc; j++)
	{
		fprintf(debug, "%03.20f", Dn[j]); // this is Dn!!!
		fprintf(debug, ", %03.20f", x[j]);
		fprintf(debug, ", %03.20f", y[j]);
        fprintf(debug, ", %03.20f", th[j]);
        for (int k=0;k<n;k++) fprintf(debug, ", %03.20f", r_shape[j*n+k]);
        for (int k=0;k<n;k++) fprintf(debug, ", %03.20f", th_shape[j*n+k]);
        fprintf(debug, ", %d", Cn[j]);
		fprintf(debug, "\n");
	}

    fclose(debug);

    // Compression jamming
    // while (P < Ptol || P > (1.0 + tol) * Ptol)
    // need to fix bug here
    while (U < Utol || U > 2.0 * Utol)
    {
        if (U < Utol)
        // if (P < Ptol)
        {
            dphi = fabs(dphi);
            // copy all inputs
            copy(Dn, Dn + Nc, Dn0);
            copy(R_eff, R_eff + Nc, R_eff0);
            copy(m, m + Nc, m0);
            copy(I, I + Nc, I0);
            copy(x, x + Nc, x0);
            copy(y, y + Nc, y0);
            copy(th, th + Nc, th0);
            copy(vx, vx + Nc, vx0);
            copy(vy, vy + Nc, vy0);
            copy(w, w + Nc, w0);
            copy(ax_old, ax_old + Nc, ax_old0);
            copy(ay_old, ay_old + Nc, ay_old0);
            copy(alph_old, alph_old + Nc, alph_old0);
            copy(x_shape, x_shape + N, x_shape0);
            copy(y_shape, y_shape + N, y_shape0);
            copy(r_shape, r_shape + N, r_shape0);
            copy(th_shape, th_shape + N, th_shape0);
        }
        // else if (P > (1.0 + tol) * Ptol)
        // else if (U > 2 * Utol)
        else if (U > 2 * Utol && C >= (2 * Nc - 1))
        {
            dphi = -fabs(dphi) / 2;
            // copy all inputs
            copy(Dn0, Dn0 + Nc, Dn);
            copy(R_eff0, R_eff0 + Nc, R_eff);
            copy(m0, m0 + Nc, m);
            copy(I0, I0 + Nc, I);
            copy(x0, x0 + Nc, x);
            copy(y0, y0 + Nc, y);
            copy(th0, th0 + Nc, th);
            copy(vx0, vx0 + Nc, vx);
            copy(vy0, vy0 + Nc, vy);
            copy(w0, w0 + Nc, w);
            copy(ax_old0, ax_old0 + Nc, ax_old);
            copy(ay_old0, ay_old0 + Nc, ay_old);
            copy(alph_old0, alph_old0 + Nc, alph_old);
            copy(x_shape0, x_shape0 + N, x_shape);
            copy(y_shape0, y_shape0 + N, y_shape);
            copy(r_shape0, r_shape0 + N, r_shape);
            copy(th_shape0, th_shape0 + N, th_shape);
        }
        U = bumpy_2D_comp_VV(Nc, N, n, gam, Dn, r_shape, th_shape, m, I, R_eff, x, y, th, Fx, Fy, T, Cn, vx, vy, w, ax, ay, alph, ax_old, ay_old, alph_old, dt, Nt, K, Lx, Ly, Utol, dphi);
        P = bumpy_2D_stress(Nc, N, n, x, y, th, r_shape, th_shape, Fx, Fy, stress, K, R_eff, Dn, Lx, Ly, gam);
        Kmax = getEk(Nc, m, I, vx, vy, w);
        temp = 0;
        count += 1;
        // get C
        C = bumpy_2D_shrjam_getC(Nc, N, n, x, y, th, r_shape, th_shape, K, R_eff, Dn, Lx, Ly, gam);
        for(int i = 0; i < Nc; i++) temp += m[i];
        phitot = temp / (Lx * Ly);
        if (count % 10 == 0) printf("Step %d, phi=%1.7f, dphi=%e, C=%d, U/K/N=%e, P=%e\n", count, phitot, dphi, C, U, P);
        if (count > count_max) break;
    }

    printf("U/K/N=%e, Kmax=%e\n", U, Kmax);

    // print final compression jammed state
    fprintf(out, "%d\n", Nc);
    fprintf(out, "%d\n", n);
    fprintf(out, "%d\n", N);
    fprintf(out, "%3.4f\n", mu);
    fprintf(out, "%03.16f\n", Lx);
    fprintf(out, "%03.16f\n", Ly);
    fprintf(out, "%03.20f\n", gam);
    for (int j = 0; j < Nc; j++)
	{
		// fprintf(out, "%3d, ", j);
		fprintf(out, "%03.20f", Dn[j]); // this is Dn!!!
		fprintf(out, ", %03.20f", x[j]);
		fprintf(out, ", %03.20f", y[j]);
        fprintf(out, ", %03.20f", th[j]);
        for (int k=0;k<n;k++) fprintf(out, ", %03.20f", r_shape[j*n+k]);
        for (int k=0;k<n;k++) fprintf(out, ", %03.20f", th_shape[j*n+k]);
        fprintf(out, ", %d", Cn[j]);
		fprintf(out, "\n");
	}
    printf("\n");
    
    fclose(out);

    delete[] xval, yval, lengths;
    delete[] Dn0, R_eff0, m0, I0, x0, y0, th0, vx0, vy0, w0, ax_old0, ay_old0, alph_old0;
    delete[] x_shape0, y_shape0, r_shape0, th_shape0;
    delete[] ax, ay, alph;
    delete[] Dn, R_eff, m, I, x, y, th, vx, vy, w, ax_old, ay_old, alph_old;
    delete[] x_shape, y_shape, r_shape, th_shape;
    
    delete[] xlist, ylist, x_mc, y_mc, in_mc;

    return 0;
}