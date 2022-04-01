#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>
#include <math.h>
#include <bitset>
#include <map>
#include <deque>
#include <unordered_map>
#include <stack>
#include <unistd.h>
#include <queue>
#include <time.h>

#include <fstream>
#include <stdio.h>
#include <string.h>

using namespace std;

using ll = long long;
using ld = long double;

vector<vector<ll>> input1, input2;

struct point{
    int x;
    int y;
    point(int x, int y) : x(x), y(y){};
};

vector<vector<ll>> trivial_mat_mult(int x1, int y1, int x2, int y2, int len) {
    vector<vector<ll>> ans(len);
    for(int i = 0; i < len; i++) {
        ans[i].resize(len, 0);
        for(int j = 0; j < len; j++) {
            for(int pos = 0; pos < len; pos++) {
                ans[i][j] += input1[y1 + i][x1 + pos] * input2[y2 + pos][x2 + j];
            }
        }
    }
    return ans;
}

void add_elements_inside(point p1, point p2, int len, int factor, int mat_number) {
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            if (mat_number == 1) {
                input1[p1.y + j][p1.x + i] += factor * input1[p2.y + j][p2.x + i];
            } else{
                input2[p1.y + j][p1.x + i] += factor * input2[p2.y + j][p2.x + i];
            }
        }
    }
}

vector<vector<ll>> strassen_mat_mult(point p1, point p2, int len, int split_point) {
    if(len <= split_point) {
        return trivial_mat_mult(p1.x, p1.y, p2.x, p2.y, len);
    }
    int x1 = p1.x;
    int x2 = p2.x;
    int y1 = p1.y;
    int y2 = p2.y;

    int mid_x1 = x1 + len / 2;
    int mid_y1 = y1 + len / 2;
    int mid_x2 = x2 + len / 2;
    int mid_y2 = y2 + len / 2;

    point A(x1, y1);
    point B(mid_x1, y1);
    point C(x1, mid_y1);
    point D(mid_x1, mid_y1);

    point E(x2, y2);
    point F(mid_x2, y2);
    point G(x2, mid_y2);
    point H(mid_x2, mid_y2);

    //we indicated additional operations for "transforming back" to out initial matrix with (*)
    // those give us theoretical bound of 25 rather than 16

    add_elements_inside(F, H, len / 2, -1, 2);
    vector<vector<ll>> prod1 = strassen_mat_mult(A, F, len / 2, split_point);
    add_elements_inside(F, H, len / 2, +1, 2); // (*)

    add_elements_inside(A, B, len / 2, +1, 1);
    vector<vector<ll>> prod2 = strassen_mat_mult(A, H, len / 2, split_point);
    add_elements_inside(A, B, len / 2, -1, 1); // (*)

    add_elements_inside(C, D, len / 2, +1, 1);
    vector<vector<ll>> prod3 = strassen_mat_mult(C, E, len / 2, split_point);
    add_elements_inside(C, D, len / 2, -1, 1); // (*)

    add_elements_inside(G, E, len / 2, -1, 2);
    vector<vector<ll>> prod4 = strassen_mat_mult(D, G, len / 2, split_point);
    add_elements_inside(G, E, len / 2, +1, 2); // (*)

    add_elements_inside(A, D, len / 2, +1, 1);
    add_elements_inside(E, H, len / 2, +1, 2);
    vector<vector<ll>> prod5 = strassen_mat_mult(A, E, len / 2, split_point);
    add_elements_inside(A, D, len / 2, -1, 1); // (*)
    add_elements_inside(E, H, len / 2, -1, 2); // (*)

    add_elements_inside(B, D, len / 2, -1, 1);
    add_elements_inside(G, H, len / 2, +1, 2);
    vector<vector<ll>> prod6 = strassen_mat_mult(B, G, len / 2, split_point);
    add_elements_inside(B, D, len / 2, +1, 1); // (*)
    add_elements_inside(G, H, len / 2, -1, 2); // (*)

    add_elements_inside(C, A, len / 2, -1, 1);
    add_elements_inside(E, F, len / 2, +1, 2);
    vector<vector<ll>> prod7 = strassen_mat_mult(C, E, len / 2, split_point);
    add_elements_inside(C, A, len / 2, +1, 1); // (*)
    add_elements_inside(E, F, len / 2, -1, 2); // (*)

    vector<vector<ll>> ans(len);
    for(int i = 0; i < len; i++){
        ans[i].resize(len, 0);
        for (int j = 0; j < len; ++j) {
            if (i >= len / 2 && j < len / 2) {
                ans[i][j] += prod3[i - len / 2][j] + prod4[i - len / 2][j];
                continue;
            }
            if (i >= len / 2 && j >= len / 2) {
                ans[i][j] += prod1[i - len / 2][j - len / 2] - prod3[i - len / 2][j - len / 2] +
                             prod5[i - len / 2][j - len / 2] + prod7[i - len / 2][j - len / 2];
                continue;
            }
            if (i < len / 2 && j < len / 2) {
                ans[i][j] += -prod2[i][j] + prod4[i][j] + prod5[i][j] + prod6[i][j];
                continue;
            }
            if (i < len / 2 && j >= len / 2) {
                ans[i][j] += prod1[i][j - len / 2] + prod2[i][j - len / 2];
                continue;
            }
        }
    }
    return ans;
}

ll count_triangles(ld p, int n) {
    input1.clear(), input2.clear();
    input1.resize(n), input2.resize(n);
    for(int i = 0; i < n; i++) {
        input1[i].resize(n, 0), input2[i].resize(n, 0);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; j++) {
            auto indicator = (long double) rand() / RAND_MAX;
            if (indicator < p) {
                input1[i][j] = 1;
            } else {
                input1[i][j] = 0;
            }
            input1[j][i] = input1[i][j];

            input2[i][j] = input1[i][j];
            input2[j][i] = input1[i][j];
        }
    }
    vector<vector<ll>> squared = strassen_mat_mult(point(0, 0), point(0,0), n, 64);
    for(int i = 0; i < n; i++) {
        for (int j = 0; j < n; ++j) {
            input2[i][j] = squared[i][j];
        }
    }
    vector<vector<ll>> cubed = strassen_mat_mult(point(0, 0), point(0,0), n, 64);
    ll ans = 0;
    for (int i = 0; i < n; ++i) {
        ans += cubed[i][i];
    }
    return ans / 6;
}

int find_optimal_split(int dimension) {
    int n = 1;
    while (n < dimension) {
        n *= 2;
    }
    input1.resize(n), input2.resize(n);
    for (int i = 0; i < n; ++i) {
        input1[i].resize(n, 0), input2[i].resize(n, 0);
        for (int j = 0; j < n; j++) {
            input1[i][j] = rand() % 5;
            input2[i][j] = rand() % 5;
        }
    }

    clock_t tStart_trivial = clock();
    vector<vector<ll>> res_trivial = trivial_mat_mult(0, 0, 0, 0, n);
    double time_min = (double) (clock() - tStart_trivial) / CLOCKS_PER_SEC;
    int min_index = n;

    for (int i = 1; i < n; i *= 2) {
        clock_t tStart = clock();
        vector<vector<ll>> res_strassen = strassen_mat_mult(point(0, 0), point(0, 0), n, i);
        double time_strassen = (double) (clock() - tStart) / CLOCKS_PER_SEC;

        if (time_strassen <= time_min) {
            min_index = i;
            time_min = time_strassen;
        }
    }
    return min_index;
}

int find_first_cross(int dimension) {
    int n = 1;
    while (n < dimension) {
        n *= 2;
    }
    input1.resize(n), input2.resize(n);
    for (int i = 0; i < n; ++i) {
        input1[i].resize(n, 0), input2[i].resize(n, 0);
        for (int j = 0; j < n; j++) {
            input1[i][j] = rand() % 5;
            input2[i][j] = rand() % 5;
        }
    }

    clock_t tStart_trivial = clock();
    vector<vector<ll>> res_trivial = trivial_mat_mult(0, 0, 0, 0, n);
    double time_trivial = (double) (clock() - tStart_trivial) / CLOCKS_PER_SEC;

    int l = 2, r = n;
    while (r - l > 1) {
        int split = (r + l) / 2;
        clock_t tStart = clock();
        vector<vector<ll>> res_strassen = strassen_mat_mult(point(0, 0), point(0, 0), n, split);
        double time_strassen = (double) (clock() - tStart) / CLOCKS_PER_SEC;
        if (time_strassen <= time_trivial) {
            r = split;
        } else {
            l = split;
        }
    }
    return r;
}

int main(int argc, char **argv) {
    srand(239);

    if(!strcmp(argv[1], "0")) {
        int dimension = stoi(argv[2]);
        int n = 1;
        while (n < dimension) {
            n *= 2;
        }

        ifstream input_file(argv[3]);
        if (!input_file.is_open()) {
            cout << "Incorrect file name" << endl;
            return 0;
        }

        input1.resize(n), input2.resize(n);
        for (int i = 0; i < n; ++i) {
            input1[i].resize(n, 0);
            if(i < dimension) {
                for (int j = 0; j < dimension; j++) {
                    input_file >> input1[i][j];
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            input2[i].resize(n, 0);
            if (i < dimension) {
                for (int j = 0; j < dimension; j++) {
                    input_file >> input2[i][j];
                }
            }
        }
        input_file.close();

        int split = min(n / 8, 64);
        if (split < 3) {
            split = 3;
        }
        vector<vector<ll>> res_strassen = strassen_mat_mult(point(0, 0), point(0, 0), n, split);
        for (int i = 0; i < dimension; ++i) {
            cout << res_strassen[i][i] << endl;
        }
        return 0;
    }

    //for triangle
    if(!strcmp(argv[1], "tr")) {
        ld p = 0.01;
        if(argc >= 3) {
            p = stod(argv[2]);
        }
        int tests = 5;
        if(argc >= 4) {
            tests = stoi(argv[3]);
        }
        int n = 1024;
        cout << "theoretically: " << (n * (n - 1) * (n - 2) / 6) * p * p * p << endl;
        cout << "experimentally: " << endl;
        ld total_res = 0;
        for(int i = 0; i < tests; i++) {
            ll res = count_triangles(p, n);
            cout << res;
            if(i < tests - 1) {
                cout << ", ";
            }
            total_res += res;
        }
        cout << endl << "average: " << total_res / tests << endl;
    }

    if(!strcmp(argv[1], "test")) {
        int n = 256;
        int tests = 5;
        while(tests--) {
            input1.resize(n), input2.resize(n);
            for (int i = 0; i < n; ++i) {
                input1[i].resize(n, 0), input2[i].resize(n, 0);
                for (int j = 0; j < n; j++) {
                    input1[i][j] = rand() % 1000;
                    input2[i][j] = rand() % 1000;
                }
            }
            vector<vector<ll>> res_strassen = strassen_mat_mult(point(0, 0), point(0, 0), n, 32);
            vector<vector<ll>> res_trivial = trivial_mat_mult(0, 0, 0, 0, n);
            for (int i = 0; i < n; ++i) {
                for(int j = 0; j < n; j++) {
                    assert(res_trivial[i][j] == res_strassen[i][j]);
                }
            }
        }
        cout << "test successful" << endl;
        return 0;
    }

    //for finding optimal split
    if(!strcmp(argv[1], "spl")) {
        if (argc >= 3) {
            cout << find_optimal_split(stoi(argv[2])) << endl;
            return 0;
        } else {
            vector<int> n_vec = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
            cout << "optimal values for different n: " << endl;
            for(int n: n_vec) {
                cout << "n = " << n << ": " << find_optimal_split(n) << endl;
            }
        }
    }
}
