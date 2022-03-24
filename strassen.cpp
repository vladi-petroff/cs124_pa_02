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

    add_elements_inside(F, H, len / 2, -1, 2);
    vector<vector<ll>> prod1 = strassen_mat_mult(A, F, len / 2, split_point);
    add_elements_inside(F, H, len / 2, +1, 2);

    add_elements_inside(A, B, len / 2, +1, 1);
    vector<vector<ll>> prod2 = strassen_mat_mult(A, H, len / 2, split_point);
    add_elements_inside(A, B, len / 2, -1, 1);

    add_elements_inside(C, D, len / 2, +1, 1);
    vector<vector<ll>> prod3 = strassen_mat_mult(C, E, len / 2, split_point);
    add_elements_inside(C, D, len / 2, -1, 1);

    add_elements_inside(G, E, len / 2, -1, 2);
    vector<vector<ll>> prod4 = strassen_mat_mult(D, G, len / 2, split_point);
    add_elements_inside(G, E, len / 2, +1, 2);

    add_elements_inside(A, D, len / 2, +1, 1);
    add_elements_inside(E, H, len / 2, +1, 2);
    vector<vector<ll>> prod5 = strassen_mat_mult(A, E, len / 2, split_point);
    add_elements_inside(A, D, len / 2, -1, 1);
    add_elements_inside(E, H, len / 2, -1, 2);

    add_elements_inside(B, D, len / 2, -1, 1);
    add_elements_inside(G, H, len / 2, +1, 2);
    vector<vector<ll>> prod6 = strassen_mat_mult(B, G, len / 2, split_point);
    add_elements_inside(B, D, len / 2, +1, 1);
    add_elements_inside(G, H, len / 2, -1, 2);

    add_elements_inside(C, A, len / 2, -1, 1);
    add_elements_inside(E, F, len / 2, +1, 2);
    vector<vector<ll>> prod7 = strassen_mat_mult(C, E, len / 2, split_point);
    add_elements_inside(C, A, len / 2, +1, 1);
    add_elements_inside(E, F, len / 2, -1, 2);

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

ll create_graph(ld p, int n) {
    input1.resize(n), input2.resize(n);
    for (int i = 0; i < n; ++i) {
        input1[i].resize(n, 0);
        input2[i].resize(n, 0);
        for (int j = 0; j < n; j++) {
            if (j != i) {
                auto indicator = (long double) rand() / RAND_MAX;
                if (indicator < p) {
                    input1[i][j] = 1;
                } else {
                    input1[i][j] = 0;
                }
                input2[i][j] = input1[i][j];
            }
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

int main(int argc, char **argv) {

    //for triangle
    if(!strcmp(argv[1], "triangle")) {
        srand(239);
        ld p = 0.01;
        int tests = 10;
        int n = 512;
        cout << "theoretically: " << (n * (n - 1) * (n - 2) / 6) * p * p * p << endl;
        cout << "experimentally: " << endl;
        ld total_res = 0;
        while(tests--){
            ll res = create_graph(p, n);
            cout << res << endl;
            total_res += res;
        }
        cout << "average: " << total_res / 10 << endl;
    }

    //for finding optimal split
    if(!strcmp(argv[1], "1")) {
        int n = stoi(argv[2]);
        input1.resize(n), input2.resize(n);
        srand(239);
        for (int i = 0; i < n; ++i) {
            input1[i].resize(n, 0);
            input2[i].resize(n, 0);
            for (int j = 0; j < n; j++) {
                input1[i][j] = rand() % 5;
            }
            for (int j = 0; j < n; j++) {
                input2[i][j] = rand() % 5;
            }
        }

        clock_t tStart_trivial = clock();
        vector<vector<ll>> res_trivial = trivial_mat_mult(0, 0, 0, 0, n);
        double time_trivial = (double) (clock() - tStart_trivial) / CLOCKS_PER_SEC;
        cout << time_trivial << endl;

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
        cout << r << endl;
        return 0;
    }

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
}

/*
//flag 0 for initial task, flag 1 for runtime and values
if(argc != 5) {
    cout << "incorrect number of input variables (expected 3, got: " << argc - 1 << ")" << endl;
    return 0;
}

if(strcmp(argv[1], "0") && strcmp(argv[1], "1") ) {
    cout << "incorrect flag" << endl;
    return 0;
}

int n = stoi(argv[2]), trials = stoi(argv[3]), dim = stoi(argv[4]);
if(dim > 4 || dim == 2) {
    cout << "incorrect dimension (possible values: 0, 2, 3, 4)" << endl;
    return 0;
}
ld avrg = 0;
if(!strcmp(argv[1], "0")) {

    return 0;
}

if(!strcmp(argv[1], "1")) {
    clock_t tStart = clock();
    vector<ld> values;
    map<string, double> times = {{"generate", 0}, {"sort", 0}, {"mst", 0}};
    for(int test = 0; test < trials; test++) {
        clock_t tStart1 = clock();

        times["generate"] += (double) (clock() - tStart1) / CLOCKS_PER_SEC;
        //printf("Time for generating graph: %.2fs\n", (double) (clock() - tStart1) / CLOCKS_PER_SEC);

        clock_t tStart2 = clock();
        sort(g.begin(), g.end());
        times["sort"] += (double) (clock() - tStart2) / CLOCKS_PER_SEC;
        //printf("Time for sorting graph: %.2fs\n", (double) (clock() - tStart2) / CLOCKS_PER_SEC);

        clock_t tStart3 = clock();
        ld ans = kruskal_mst(g, n);
        times["mst"] += (double) (clock() - tStart3) / CLOCKS_PER_SEC;
        //printf("Time for find MST (Kruskal): %.2fs\n", (double) (clock() - tStart3) / CLOCKS_PER_SEC);
        values.push_back(ans);
    }

    cout << "values:" << endl;
    for (int i = 0; i < trials; ++i) {
        avrg += values[i];
        cout << values[i] << " ";
    }
    cout << endl << "Average weight = " << avrg / trials << endl;
    printf("Time statistics: \n");
    printf("Average time for generating graph: %.2fs\n", times["generate"] / trials);
    printf("Average time for sorting graph: %.2fs\n", times["sort"] / trials);
    printf("Average time for find MST (Kruskal): %.2fs\n", times["mst"] / trials);
    printf("Average time for all procedure: %.2fs\n", ((double)(clock() - tStart)/CLOCKS_PER_SEC) / trials);
    printf("Overall program time taken (for all trials): %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}
*/
