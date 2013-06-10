#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <sys/time.h>
using namespace std;

// flags
#define PROFILE 0
#define PRINT_FRAMES 1
#define PRINT_SCORE_UPDATES 0

// constant values
const int MAX_N = 512;
const double MAX_RUNNING_TIME = 9.5;
//const double MAX_RUNNING_TIME = 1.4;
const double INF = 1e100;
const double G = 9.8;
const double E = 0.0;
const double DECAY_PER_FRAME = 0.003;
const double TIME_PER_FRAME = 0.0005;
const int RESTART_FRAME = 7500;
const double BOUNCE_MARGIN = 0.00025;

#if PROFILE
class Profiler {
public:
    int last;
    long long acc[100];
    Profiler() {
        memset(acc, 0, sizeof(acc));
    }
    void start() {
        struct timeval  tv;
        gettimeofday(&tv, NULL);
        last = tv.tv_usec;
    }
    void end(int id) {
        struct timeval  tv;
        gettimeofday(&tv, NULL);
        int elapsed = (tv.tv_usec - last + 1000000) % 1000000;
        acc[id] += elapsed;
    }
    void report() {
        for (int i = 0; i < 100; i++) if (acc[i] > 0) {
            cerr << i << ": " << acc[i] / 1000.0 << "(ms), ";
        }
        cerr << endl;
        memset(acc, 0, sizeof(acc));
    }
};
Profiler prof;
#define PROF_START()    prof.start();
#define PROF_END(id)    prof.end(id);
#define PROF_REPORT()   prof.report();
#else
#define PROF_START()
#define PROF_END(id)
#define PROF_REPORT()
#endif
    
unsigned int randxor() {
    static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
    unsigned int t;
    t=(x^(x<<11));x=y;y=z;z=w; return( w=(w^(w>>19))^(t^(t>>8)) );
}

class Vec {
public:
    double x, y;
    Vec() {};
    Vec(double x, double y): x(x), y(y) {}
    Vec& operator+=(const Vec& o) { x += o.x; y += o.y; return *this; }
    Vec& operator-=(const Vec& o) { x -= o.x; y -= o.y; return *this; }
    Vec& operator*=(const double m) { x *= m; y *= m; return *this; }
    Vec& operator/=(const double m) { x /= m; y /= m; return *this; }
    void clear() { x=0, y=0; }
    double length() { return sqrt(x*x+y*y); }
    bool isSmallerThan(double r) { return x*x + y*y < r*r; }
    void normalize() { double l=length(); if(l!=0) {x/=l;y/=l;} }
};
//ostream& operator<<(ostream& o,const Vec& v) {o << "(" << v.x << "," << v.y << ")"; return o;}
Vec operator+(const Vec& a, const Vec& b) {return Vec(a.x+b.x,a.y+b.y);}
Vec operator-(const Vec& a, const Vec& b) {return Vec(a.x-b.x,a.y-b.y);}
Vec operator*(const Vec& a, const double m) {return Vec(a.x*m, a.y*m);}
Vec operator/(const Vec& a, const double m) {return Vec(a.x/m, a.y/m);}
double dot(const Vec& a, const Vec& b) {return a.x*b.x+a.y*b.y;}
double cross(const Vec& a, const Vec& b) {return a.x*b.y-a.y*b.x;}
int ccw(const Vec& a, const Vec& b) {double cp=cross(a, b); return cp ? (cp>0?1:-1) : 0;}

struct Circle
{
    // states;
    Vec pos, v;
    bool is_hover;

    // properties
    int index;
    double r, m;
    Vec o_pos;
    
    // precomputed values
    double inv_r2;
    double inv_m;
    double gravity;
};

bool greaterMass(const Circle& lhs, const Circle& rhs) {
    return lhs.m > rhs.m;
}

    
// state variables
int N;
Circle c[MAX_N];

// computing buffers and precomputed values
struct XS {
    double x;
    int index;
    bool is_left;
    bool operator<(const XS& rhs) const {
        return x < rhs.x;
    }
};
XS xs[MAX_N * 2];
struct Pair {
    int a, b;
};
int pairs[MAX_N * 8];
double m_mul_div_sum[MAX_N][MAX_N];

class CirclesSeparation {
private:
    double best_cost;
    vector<double> best_x, best_y;
    Circle* hover_circle;
    int frames;
    
public:
    CirclesSeparation() {
        hover_circle = NULL;
        frames = 0;
    }
    
    void setup(vector<double> x, vector<double> y, vector<double> r, vector<double> m) {
        srand(10);
        N = x.size();
        best_x = vector<double>(N);
        best_y = vector<double>(N);
        best_cost = INF;
        for (int i = 0; i < N; i++) {
            c[i].index = i;
            c[i].pos = c[i].o_pos = Vec(x[i], y[i]);
            c[i].r = r[i];
            c[i].m = m[i];
            c[i].is_hover = false;
            c[i].inv_r2 = 1.0 / r[i] / r[i];
            c[i].inv_m = 1.0 / m[i];
        }
        
        // determine initial position
        sort(c, c + N, greaterMass);

        // precompute fixed values
        for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
            c[i].gravity = G * TIME_PER_FRAME * min(3.0, 1 + 0.01 * c[i].m * c[i].inv_r2);
            m_mul_div_sum[i][j] = m_mul_div_sum[j][i] = c[i].m * c[j].m / (c[i].m + c[j].m);
        }
        restart();
    }
    
    void restart() {
        for (int i = 0; i < N; i++) {
            c[i].v.clear();
            c[i].pos = c[i].o_pos;
        }
        
        // determine initial position
        vector<pair<double, int> > scores(N);
        double totalArea = 0;
        for (int i = 0; i < N; i++) {
            double mass = c[i].m;
            double area = c[i].r * c[i].r;
            double jama = area / mass / mass;
            scores[i] = make_pair(jama, i);
            totalArea += 3.14*area;
        }
        
        sort(scores.begin(), scores.end());
        
        // adjust initial filling rate
        double thresh = 0.8 + 0.4 * (double)rand() / RAND_MAX;
        double filledArea = 0.0;
        vector<pair<double, int> > outers;
        for (int i = 0; i < N; i++) {
            int index = scores[i].second;
            if (filledArea > thresh) {
                outers.push_back(pair<double, int>(-c[i].m * c[i].inv_r2, index));
            } else {
                filledArea += c[index].r * c[index].r * 3.14;
            }
        }
        for (int i = 0; i < outers.size(); i++) {
            int index = outers[i].second;
            Vec d = c[index].pos - Vec(0.5, 0.5);
            d.normalize();
            c[index].pos += d * (3 + ((double)i / outers.size()) * 2);
        }
    }
    
    void update() {
        frames++;
        
        ///////////////////////////////////////////////////////
        // apply force
        PROF_START();
        for (int i = 0; i < N; i++) {
            Vec d = c[i].o_pos - c[i].pos;
            d.normalize();
            //double mul = G * TIME_PER_FRAME * min(3.0, 1 + 0.01 * c[i].m * c[i].inv_r2);
            c[i].v += d * c[i].gravity;
        }
        PROF_END(0);
        
        // shake
        static int shake_dir = 0;
        static const int DX[] = {-1, 1, 0, 0};
        static const int DY[] = {0 ,0, -1, 1};
        if (frames % 2000 == 0) {
            for (int i = 0; i < N; i++) {
                c[i].v += Vec(DX[shake_dir], DY[shake_dir]);
            }
            shake_dir = (shake_dir + 1) % 4;
        }
        

        ///////////////////////////////////////////////////////
        // broad phase
        PROF_START();
        int num_xs = 0;
        for (int i = 0; i < N; i++) {
            xs[num_xs].x = c[i].pos.x - c[i].r - BOUNCE_MARGIN;
            xs[num_xs].index = i;
            xs[num_xs].is_left = true;
            num_xs++;
            xs[num_xs].x = c[i].pos.x + c[i].r + BOUNCE_MARGIN;
            xs[num_xs].index = i;
            xs[num_xs].is_left = false;
            num_xs++;
        }
        sort(xs, xs + num_xs);
        PROF_END(1);
        
        PROF_START()
        int num_pairs = 0;
        for (int i = 0; i < num_xs; i++) {
            if (xs[i].is_left) {
                int a = xs[i].index;
                for (int j = i+1; j < num_xs; j++) {
                    int b = xs[j].index;
                    if (b == a) break;
                    if (xs[j].is_left && abs(c[a].pos.y - c[b].pos.y) < c[a].r + c[b].r + BOUNCE_MARGIN*2) {
                        pairs[num_pairs++] = (a << 16) | b;
                    }
                }
            }
        }
        PROF_END(2);

        ///////////////////////////////////////////////////////
        // detect collision and solve constraints
        PROF_START();
        const double CD = 0.2 / TIME_PER_FRAME;
        for (int k = 0; k < num_pairs; k++) {
            int i = pairs[k] >> 16;
            int j = pairs[k] & ((1<<16)-1);
            if (c[i].is_hover || c[j].is_hover) continue;
            Vec d = c[j].pos - c[i].pos;
            if (d.isSmallerThan(c[i].r + c[j].r + BOUNCE_MARGIN*2)) {
                Vec dv = c[j].v - c[i].v;
                double len = d.length();
                Vec norm = d / len;
                double CDD = CD * (c[i].r + c[j].r  + BOUNCE_MARGIN*2 - len);
                double C = m_mul_div_sum[i][j] * ((1 + E) * dot(dv, norm) - CDD);
                c[i].v += norm * C * c[i].inv_m;
                c[j].v -= norm * C * c[j].inv_m;
            }
        }
        PROF_END(3);
        
        // add friction
        for (int i = 0; i < N; i++) {
            c[i].v *= 1 - DECAY_PER_FRAME;
        }


        ///////////////////////////////////////////////////////
        // integrate
        PROF_START();
        for (int i = 0; i < N; i++) {
            c[i].pos += c[i].v * TIME_PER_FRAME;
        }
        PROF_END(4);
        
        ///////////////////////////////////////////////////////
        // others
        
        // kick circles which have overlap
        PROF_START();
    /*
        if (frames % 300 == 0) {
            bool need_shake = false;
            for (int i = N-1; i >= 0; i--) {
                for (int j = 0; j < i; j++) {
                    if ((c[i].pos - c[j].pos).isSmallerThan(c[i].r + c[j].r)) {
     
                        int moved = j;
                        if (c[j].m * c[j].inv_r2 > c[i].m * c[i].inv_r2) moved = i;
                        Vec d = c[moved].pos - Vec(0.5, 0.5);
                        d.normalize();
                        c[moved].pos += d * 1;
                    }
                }
            }
        OUTER:
            if (need_shake) 
        }
        */
        
        if (frames % 100 == 0)
            updateBest();

        
        if (frames % RESTART_FRAME == 0) {
            restart();
        }
        PROF_END(5);
        
        if (frames % 10000 == 0) {
            PROF_REPORT();
        }
    }
    
    double currentCost() {
        double cost = 0;
        for (int i = 0; i < N; i++) {
            Vec d = c[i].pos - c[i].o_pos;
            cost += d.length() * c[i].m;
        }
        return cost;
    }
                
    int currentFrame() {
        return frames;
    }
        
    void catch_circle(Vec pos) {
        for (int i = 0; i < N; i++) {
            if ((pos - c[i].pos).isSmallerThan(c[i].r)) {
                c[i].is_hover = true;
                c[i].v.clear();
                hover_circle = &c[i];
                break;
            }
        }
    }
    
    void move_circle(Vec pos) {
        if (hover_circle) {
            hover_circle->pos = pos;
        }
    }
    
    void release_circle() {
        if (hover_circle) {
            hover_circle->is_hover = false;
            hover_circle = NULL;
        }
    }
    
    vector<double> minimumWork(vector<double> x, vector<double> y, vector<double> r, vector<double> m)
    {
        clock_t end__ = clock() + MAX_RUNNING_TIME * CLOCKS_PER_SEC;
        
        setup(x, y, r, m);
        while (clock() < end__) for (int i=0; i<100; i++) update();
        
#if PRINT_FRAMES
        cerr << "frames: " << frames << endl;
#endif
        
        vector<double> res;
        for (int i = 0; i < N; i++) {
            res.push_back(best_x[i]);
            res.push_back(best_y[i]);
        }
        return res;
    }
    
private:
    
    bool isOverlap(int i, int j) {
        if (i == j) return false;
        return (c[i].pos - c[j].pos).isSmallerThan(c[i].r + c[j].r);
    }
    
    
    bool hasOverlap() {
        for (int i = 0; i < N; i++)
            for (int j = i+1; j < N; j++)
                if (isOverlap(i, j))
                    return true;
        return false;
    }
    
    void updateBest() {
        if (hasOverlap()) return;
        
        double cost = currentCost();
        if (best_cost > cost) {
#if PRINT_SCORE_UPDATES
            cerr << "cost: " << cost << endl;
#endif
            best_cost = cost;
            for (int i = 0; i < N; i++) {
                best_x[c[i].index] = c[i].pos.x;
                best_y[c[i].index] = c[i].pos.y;
            }
        }
    }
};
