#include <iostream>
#include <ginac/ginac.h>
#include <stdlib.h>
#include <complex>
#include <vector>

using namespace std;
using namespace GiNaC;
using cmplx = complex<long double>;

const numeric h = 0.0001;
const numeric q0 = 0.5;
const numeric ERROR = 0.001;

numeric rk2(vector<numeric> &points, numeric initialQ, numeric initialE, numeric h, numeric start, numeric stop, ex dt, symbol t, symbol s)
{

    int test = 0;
    numeric q = initialQ;
    ex e = initialE;
    // points.push_back(initialE);

    while (q < stop)
    {

        ex m = dt.subs(lst{t == e, s == q});
        numeric newH = h;
        numeric nextq = q + newH;
        
        ex eHat = e + m * newH;
        ex n = dt.subs(lst{t == eHat, s == q});
        ex nexteEx = e + (m + n) / 2 * newH;
        q = nextq;
        e = nexteEx;
        test++;
    }
    return(ex_to<numeric>(e));
    // cout << e<< " TEST\n";
    q = initialQ;
    e = initialE;

    while (q > start)
    {
        ex m = dt.subs(lst{t == e, s == q});
        numeric newH = h;

        numeric nextq = q - newH;
        ex eHat = e - m * newH;
        ex n = dt.subs(lst{t == eHat, s == q});
        ex nexteEx = e + (m + n) / 2 * newH;
        numeric nexte = ex_to<numeric>(evalf(e + (m + n) / 2 * newH));
        points.push_back(nexte);
        q = nextq;
        e = nexte;
        test++;
    }

    // cout << points[points.size() - 1]<< "\n";
}


// int rk2(vector<numeric> &points, numeric initialQ, numeric initialE, numeric h, numeric start, numeric stop, ex dt, symbol t, symbol s)
// {
//     int test = 0;
//     numeric q = initialQ;
//     numeric e = initialE;
//     points.push_back(initialE);
//     ex copydt = dt;

//     while (q < stop)
//     {   
//         dt = copydt;
//         numeric nextq = q + h;
//         ex nexte = e + (dt.subs(lst{s == q, t == e}) * h);
//         e = ex_to<numeric>(nexte);
//         points.push_back(e);
//         q = nextq;
//         test++;
//     }

//     q = initialQ;
//     e = initialE;

//     while (q > start)
//     {
//         numeric nextq = q - h;
//         ex nexte = e - dt.subs(lst{t == e, s == q}) * h;
//         e = ex_to<numeric>(evalf(nexte));
//         points.push_back(e);
//         q = nextq;
//         test++;
//     }
//     cout << test << "\n";
//     cout << points[points.size() - 1]<< "\n";
// }

int main(int argc, char **argv)
{
    symbol t, s;
    symtab table;
    table["t"] = t;
    table["s"] = s;

    parser reader(table);
    ex bmap = reader(argv[1]);
    int degree = bmap.degree(t);

    // Define helper functions
    ex ft = bmap.numer();
    ex gt = bmap.denom();
    ex bt = ft - (q0 * gt);

    // Get random polynomial of known roots
    ex at = 1;
    std::srand(std::time(nullptr));
    vector<numeric> roots{};
    for (int i = 0; i < degree; i++)
    {   

        numeric randC = (double) rand() / RAND_MAX  + ((double) rand() / RAND_MAX) * I;
        roots.push_back(randC);
        ex e = t - (randC);
        at = at * e;
    }

    // Derivitive with respect to t
    ex dt = (at - bt) / ((1 - s) * at.diff(t) + s * bt.diff(t));

    // Compute the points
    vector<numeric> points{};
    for (auto t0 : roots){
        cout << t0 << "\n";
        numeric root = rk2(points, 0, t0, h, 0, 1, dt, t, s);
        cout << "ROOT: " << root << "\n";
        cout << "CHECK: " << abs(evalf(bmap.subs(t == root))) << "\n";
    }
    
    cout << bt<< endl;
    return 0;
}
