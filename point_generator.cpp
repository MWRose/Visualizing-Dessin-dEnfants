#include <iostream>
#include <complex>
#include <vector>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
using namespace std;
using cmplx = complex<long double>;

const double h = 0.000001;
const int n = 3;
const int black[2]{-1, 2};
const int white[2]{-2, 1};

// Midpoints for simple map
// const cmplx midpoints[3]{
//     cmplx (-1 * sqrt(3), 0),
//     cmplx (0, 0),
//     cmplx (sqrt(3), 0)
// };



// Simple
// cmplx bMap(cmplx t){
//     return(pow((cmplx(1.0,0) + t), 2) * (cmplx(2, 0) - t) / cmplx(4, 0));
// }
cmplx bMap(cmplx t){
    return(pow(t, n));
}


cmplx dt(cmplx t, double h){
    return(cmplx (3, 0) * pow(t, 2));
}

int computePoints(vector<cmplx> & points, vector<cmplx> & midpoints){
    for (int i = 0; i < midpoints.size(); i++){
        int test = 0;
        cmplx mp = midpoints[i];
        double q = 0.5;
        cmplx e = mp;
        points.push_back(mp);
        while (q < 1.0){
            double nextq = q + abs(dt(e, h)) * h;
            cmplx nexte = e + abs(dt(e, h)) / dt(e, h) * cmplx(h, 0);
            points.push_back(nexte);
            q = nextq;
            e = nexte;
            test++;
        }

        q = 0.5;
        e = mp;

        while (q > 0){
            double nextq = q - abs(dt(e, -1 * h)) * h;
            cmplx nexte = e - abs(dt(e, -1 * h)) / dt(e, -1 * h) * cmplx(h, 0);
            points.push_back(nexte);
            q = nextq;
            e = nexte;
            test++;
        }
        cout << test << "\n";
    }
    return(0);
}

pcl::PointXYZ stereoProject(cmplx c){
    double x = 2 * c.real() / (pow(abs(c), 2) + 1);
    double y = 2 * c.imag() / (pow(abs(c), 2) + 1);
    double z = (pow(abs(c), 2) - 1) / (pow(abs(c), 2) + 1);
    pcl::PointXYZ basic_point;
    basic_point.x = x;
    basic_point.y = y;
    basic_point.z = z;
    return(basic_point);
}

int createCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, vector<cmplx> & points){
    for (auto item: points){
        cloud->points.push_back(stereoProject(item));
    }
    cloud->width = cloud->size ();
    cloud->height = 1;
}


int main(int argc, char* argv[])
{  
    vector<cmplx> midpoints = {
    cmplx (-0.396850262992049868687926409, 0.6873648184993013131917),
    cmplx (-0.39685026299204986868792640981, -0.687364818499301313191739598443),
    cmplx (0.7937005259840997373758528196, 0)
    };
    vector<complex<long double>> points;
    pcl::PointCloud<pcl::PointXYZ>::Ptr basic_cloud_ptr (new pcl::PointCloud<pcl::PointXYZ>);
    computePoints(points, midpoints);
    createCloud(basic_cloud_ptr, points);
    pcl::io::savePCDFileASCII ("test_pcd.pcd", *basic_cloud_ptr);
    cout << basic_cloud_ptr->points.size() << " points created\n";
    return(0) ;
}