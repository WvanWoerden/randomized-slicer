//g++ short_path.cpp -O3 -o short_path -fopenmp
#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <math.h>
#include <fstream>
#include <omp.h>

using namespace std;

typedef pair<double, int> pii;
typedef pair<double, pair<int,int>> qpii;

double dijkstra( vector<vector<pii>> G, int s, int t ) {
    vector<double> d(G.size(),-1);
    vector<double> dp(G.size(), -1);
    vector<int> steps(G.size(), -1);
	d[s] = 0;
    dp[s] = 0;
    steps[s] = 0;
	priority_queue<qpii,vector<qpii>,greater<qpii>> Q;
	int b = s, c, step;
    double m;
	do {
		for( auto X : G[b] ) {
			pair<int,int> tmp = { steps[b]+1, X.second };
            if( d[X.second] < -0.5 and (dp[X.second] < -0.5 or d[b] + X.first < dp[X.second]) ) {
                Q.emplace( d[b] + X.first, tmp );
                dp[X.second] = d[b] + X.first;
            }
        }
        do {
			m = Q.top().first;
            step = Q.top().second.first;
			c = Q.top().second.second;
			Q.pop();
		} while( d[c] > -0.5 );
		d[c] = m;
        steps[c] = step;
		b = c;
	} while( b != t );
	return d[t];
}

double P( double a, double u){
    if( std::fabs(u) > 1 )
        return -999;
    else {
        double tmp = a*sqrt(1-u*u);
        double p= log2(min(tmp, double(1)));
        return p;
    }
}

double U(double a, double x, double y) {
    return (a*a+x*x-y*y)/(2*x*a);
}

double B(double a){
    return sqrt(a*a*a*a/(4*a*a-4));
}

vector<vector<pii>> make_graph_cvp(double a, int grid, double c=1){
    double b = B(a);
    int n = grid+1;
    vector<double> S(n);
    for( int i = 0; i <= grid; i++ ) {
        S[i] = sqrt( c*c+(b*b-c*c)*i/double(grid) );
    }

    vector<vector<pii>> G(n);
    for( int i = 0; i < n; i++ ) {
        G[i].reserve(grid);
        for( int j = 0; j < n; j++ ) {
            if( i==j)
                continue;
            double x = S[i];
            double y = S[j];
            double u = U(a,x,y);
            double p = max(-P(a,u), double(0));
            G[i].push_back( {p, j} );
        }
    }
    return G;
}

vector<vector<pii>> make_graph_bdd(double a, double d, int grid){
    double b = B(a);
    double a2 = a*a;
    double b2 = b*b;
    double d2 = d*d;

    int n = grid+1;
    vector<double> S(n);
    for( int i = 0; i <= grid; i++ ) {
        S[i] = 1+(b-1)*i/double(grid);
    }

    vector<vector<pii>> G(n+1);
    for( int i = 0; i < n; i++ ) {
        G[i].reserve(grid);
        for( int j = 0; j < n; j++ ) {
            if( i==j)
                continue;
            double x = S[i];
            double y = S[j];
            double u = U(a,x,y);
            double p = max(-P(a,u), double(0));
            G[i].push_back( {p, j} );
        }
        
        if( S[i] * S[i] <= a2-d2 ) {
            G[i].push_back( {0.0, n} );
        }
        else if( S[i] < a + d ) {
            double s2 = S[i]*S[i];
            double p = (-s2*s2+2*s2*(d2+a2) - (a2-d2)*(a2-d2))/(4*s2*d2);
            p = -log2( min(p , double(1.)))/2.;
            G[i].push_back( {p, n });
        }
    }
    return G;
}

double Bnns(double a) {
    double gamma_max = sqrt((1-4*a*a)/(5-8*a));
    return sqrt(2-2*gamma_max);
}

double Pnns( double a, double x, double y) {
    double g_x = (1-(x*x)/2.);
    double g_y = (1-(y*y)/2.);
    return log2(min(4/3. * (1. - (a*a+g_y*g_y - 2*a*g_x*g_y)/(1-g_x*g_x)), 1.));
}

vector<vector<pii>> make_graph_nns(double a, int grid){
    double b = Bnns(a);
    int n = grid+1;
    vector<double> S(n);
    for( int i = 0; i <= grid; i++ ) {
        S[i] = 1+(b-1)*i/double(grid);
    }

    vector<vector<pii>> G(n);
    for( int i = 0; i < n; i++ ) {
        G[i].reserve(grid);
        for( int j = 0; j < n; j++ ) {
            if( i==j)
                continue;
            double x = S[i];
            double y = S[j];
            double p = max(-Pnns(a,x,y), double(0));
            G[i].push_back( {p, j} );
        }
    }
    return G;
}

int main() {
    
    int datapoints = 1000;
    int grid = 5000;

    omp_set_num_threads(6);

    // Which experiments to run
    bool cvp = true;
    bool bdd = true;
    bool nns = true;

    // ------- delta-CVP --------
    if(cvp) {
        double a1 = 1.;
        double a2 = 1.41421;
        const vector<double> deltas = {1, sqrt(4/3.), 1.2,1.3,1.5,2};
#pragma omp parallel for
for( int it = 0; it < deltas.size(); it++ ) {
            double delta = deltas[it];
            ofstream myfile;
            string filepath = "data/"+std::to_string(delta)+"-cvp-"+std::to_string(datapoints)+"-"+std::to_string(grid)+".data";
            myfile.open (filepath);
           
            cerr << "Running "+to_string(delta)+"-cvp experiments.." << endl;
            for( int i = 0; i < datapoints; i++ ) {
                double a = (a1 * i + a2 * (datapoints-i))/(datapoints);
                auto G = make_graph_cvp( a, grid, delta );
                double plog = dijkstra(G, grid, 0);
                myfile << a << " " << pow(2.,-plog) << endl;
            }
            myfile.close();
        }
    }
    
    // ------- delta-BDD ---------
    if(bdd){
        double a1 = 1.;
        double a2 = 1.41421;
        const vector<double> deltas = {0,0.2,0.4,0.6, 0.8, 1.};
#pragma omp parallel for
for( int it = 0; it < deltas.size(); it++ ) {
            double delta = deltas[it];
            ofstream myfile;
            string filepath = "data/"+std::to_string(delta)+"-bdd-"+std::to_string(datapoints)+"-"+std::to_string(grid)+".data";
            myfile.open (filepath);
           
            cerr << "Running "+to_string(delta)+"-bdd experiments.." << endl;
            for( int i = 0; i < datapoints; i++ ) {
                double a = (a1 * i + a2 * (datapoints-i))/(datapoints);
                auto G = make_graph_bdd( a, delta, grid);
                double plog = dijkstra(G, grid, grid+1);
                myfile << a << " " << pow(2.,-plog) << endl;
            }
            myfile.close();
        }
    }
    
    // ------- NNS ---------
    if(nns) {
        double a1 = 0.4;
        double a2 = 0.5;
        ofstream myfile;
        string filepath = "data/nns-"+std::to_string(datapoints)+"-"+std::to_string(grid)+".data";
        myfile.open (filepath);
       
        cerr << "Running nns experiments.." << endl;
        for( int i = 0; i < datapoints; i++ ) {
            double a = (a1 * i + a2 * (datapoints-i))/(datapoints);
            auto G = make_graph_nns( a, grid);
            double plog = dijkstra(G, grid, 0);
            myfile << a << " " << pow(2.,-plog) << endl;
        }
        myfile.close();
    }

    return 0;
}
