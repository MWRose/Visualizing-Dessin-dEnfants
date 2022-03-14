#include <iostream>
#include <complex>
#include <vector>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
using namespace std;
using cmplx = complex<long double>;

const double h = 0.000001;
const int n = 5;
const int black[2]{-1, 2};
const int white[2]{-2, 1};

cmplx bMap(cmplx t)
{
    return (pow(t, n));
}

cmplx dt(cmplx t)
{
    return (cmplx(5, 0) * pow(t, 4));
}

double hcalc(cmplx t)
{
    return (abs(dt(t)) * h);
}

void euler(vector<cmplx> &points, double initialQ, cmplx initialE, function<double(cmplx)> hcalc, double start, double stop)
{
    int test = 0;
    double q = initialQ;
    cmplx e = initialE;
    points.push_back(initialE);
    while (q < stop)
    {
        double new_h = hcalc(e);
        double nextq = q + new_h;
        cmplx nexte = e + cmplx(1, 0) / dt(e) * cmplx(new_h, 0);
        points.push_back(nexte);
        q = nextq;
        e = nexte;
        test++;
    }

    q = initialQ;
    e = initialE;

    while (q > start)
    {
        double new_h = hcalc(e);
        double nextq = q - new_h;
        cmplx nexte = e - cmplx(1, 0) / dt(e) * cmplx(new_h, 0);
        points.push_back(nexte);
        q = nextq;
        e = nexte;
        test++;
    }

    cout << test << "\n";
}

int rk2(vector<cmplx> &points, double initialQ, cmplx initialE, function<double(cmplx)> hcalc, double start, double stop)
{
    int test = 0;
    double q = initialQ;
    cmplx e = initialE;
    points.push_back(initialE);
    while (q < stop)
    {
        cmplx m = cmplx(1, 0) / dt(e);
        double new_h = hcalc(e);

        double nextq = q + new_h;
        cmplx eHat = e + m * cmplx(new_h, 0);
        cmplx n = cmplx(1, 0) / dt(eHat);
        cmplx nexte = e + (m + n) / cmplx(2, 0) * cmplx(new_h, 0);

        points.push_back(nexte);
        q = nextq;
        e = nexte;
        test++;
    }
    while (q > start)
    {
        cmplx m = cmplx(1, 0) / dt(e);
        double new_h = hcalc(e);

        double nextq = q - new_h;
        cmplx eHat = e - m * cmplx(new_h, 0);
        cmplx n = cmplx(1, 0) / dt(eHat);
        cmplx nexte = e - (m + n) / cmplx(2, 0) * cmplx(new_h, 0);

        points.push_back(nexte);
        q = nextq;
        e = nexte;
        test++;
    }
    cout << test << "\n";
}

int computePoints(vector<cmplx> &points, vector<cmplx> &midpoints, function<double(cmplx)> hcalc)
{
    for (int i = 0; i < midpoints.size(); i++)
    {
        int test = 0;
        cmplx mp = midpoints[i];
        rk2(points, 0.5, mp, hcalc, 0, 1);
    }
}

pcl::PointXYZ stereoProject(cmplx c)
{
    double x = 2 * c.real() / (pow(abs(c), 2) + 1);
    double y = 2 * c.imag() / (pow(abs(c), 2) + 1);
    double z = (pow(abs(c), 2) - 1) / (pow(abs(c), 2) + 1);
    pcl::PointXYZ basic_point;
    basic_point.x = x;
    basic_point.y = y;
    basic_point.z = z;
    return (basic_point);
}

int createCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, vector<cmplx> &points)
{
    for (auto item : points)
    {
        cloud->points.push_back(stereoProject(item));
    }
    cloud->width = cloud->size();
    cloud->height = 1;
}

int main(int argc, char *argv[])
{
    vector<cmplx> midpoints = {
        cmplx(0.269014918521185704585410204277269179132628921019573277086367947, 0.827942785987195408535269672744287830767876052483616065811332165),
        cmplx(-0.70429020016924777415354521301714222862219163319357480120729773, 0.51169678248036692509229257895244838117096130879945932574129400),
        cmplx(-0.70429020016924777415354521301714222862219163319357480120729773, -0.51169678248036692509229257895244838117096130879945932574129400),
        cmplx(0.269014918521185704585410204277269179132628921019573277086367947, -0.827942785987195408535269672744287830767876052483616065811332165),
        cmplx(0.8705505632961241391362700174797460989791254243480030482418595685, 0)
        };
    vector<complex<long double>> points;
    pcl::PointCloud<pcl::PointXYZ>::Ptr basic_cloud_ptr(new pcl::PointCloud<pcl::PointXYZ>);
    computePoints(points, midpoints, hcalc);
    createCloud(basic_cloud_ptr, points);
    pcl::io::savePCDFileASCII("test_pcd.pcd", *basic_cloud_ptr);
    cout << basic_cloud_ptr->points.size() << " points created\n";
    return (0);
}