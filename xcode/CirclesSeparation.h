#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <sys/time.h>
using namespace std;

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
bool operator==(const Vec& a, const Vec& b) {return a.x==b.x && a.y==b.y;}
bool operator!=(const Vec& a, const Vec& b) {return !(a==b); }
bool operator< (const Vec& a, const Vec& b) {return a.x<b.x || (!(b.x<a.x) && a.y<b.y);}
bool operator<=(const Vec& a, const Vec& b) {return !(b<a);}
bool operator> (const Vec& a, const Vec& b) {return b<a;}
bool operator>=(const Vec& a, const Vec& b) {return !(a<b);}
ostream& operator<<(ostream& o,const Vec& v) {o <<"(" << v.x << "," << v.y << ")"; return o;}
Vec operator+(const Vec& a, const Vec& b) {return Vec(a.x+b.x,a.y+b.y);}
Vec operator-(const Vec& a, const Vec& b) {return Vec(a.x-b.x,a.y-b.y);}
Vec operator*(const Vec& a, const double m) {return Vec(a.x*m, a.y*m);}
Vec operator/(const Vec& a, const double m) {return Vec(a.x/m, a.y/m);}
double dot(const Vec& a, const Vec& b) {return a.x*b.x+a.y*b.y;}
double cross(const Vec& a, const Vec& b) {return a.x*b.y-a.y*b.x;}
int ccw(const Vec& a, const Vec& b) {double cp=cross(a, b); return cp ? (cp>0?1:-1) : 0;}

const double INF = 1e100;
const double G = 9.8;
const double E = 0.0;
const double TIME_PER_FRAME = 0.001;
const int RESTART_FRAME = 12000;

#define MAX_CIRCLE_NUM 512
pair<double, pair<int, bool> > xs[MAX_CIRCLE_NUM * 2];
pair<int, int> pairs[MAX_CIRCLE_NUM * 8];
double m_mul_div_sum[MAX_CIRCLE_NUM][MAX_CIRCLE_NUM];

#define PROFILE 0
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

unsigned int randxor() {
    static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
    unsigned int t;
    t=(x^(x<<11));x=y;y=z;z=w; return( w=(w^(w>>19))^(t^(t>>8)) );
}

struct Circle {
    int index;
    double r, m;
    Vec o_pos, pos, v;
    bool is_hover;
    double inv_r2;
    double inv_m;
    double gravity;
};

bool greaterMass(const Circle& lhs, const Circle& rhs) {
    return lhs.m > rhs.m;
}

class CirclesSeparation {
    
private:
    int N;
    double best_cost;
    vector<double> best_x, best_y;
    vector<Circle> c;
    Circle* hover_circle;
    int frames;
    Profiler prof;
    
public:
    CirclesSeparation() {
        hover_circle = NULL;
        frames = 0;
    }
    
    void setup(vector<double> x, vector<double> y, vector<double> r, vector<double> m) {
        srand(10);
        N = x.size();
        c = vector<Circle>(N);
        best_x = vector<double>(N);
        best_y = vector<double>(N);
        best_cost = INF;
        for (int i = 0; i < x.size(); i++) {
            c[i].index = i;
            c[i].pos = c[i].o_pos = Vec(x[i], y[i]);
            c[i].r = r[i];
            c[i].m = m[i];
            c[i].is_hover = false;
            c[i].inv_r2 = 1.0 / r[i] / r[i];
            c[i].inv_m = 1.0 / m[i];
        }
        
        // determine initial position
        sort(c.begin(), c.end(), greaterMass);

        // precompute
        for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
            m_mul_div_sum[i][j] = m_mul_div_sum[j][i] = c[i].m * c[j].m / (c[i].m + c[j].m);
            c[i].gravity = G * TIME_PER_FRAME * min(3.0, 1 + 0.01 * c[i].m * c[i].inv_r2);
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
        
        ////////////////////////////////
        // apply force
#if PROFILE
        prof.start();
#endif
        for (int i = 0; i < N; i++) {
            Vec d = c[i].o_pos - c[i].pos;
            d.normalize();
            //double mul = G * TIME_PER_FRAME * min(3.0, 1 + 0.01 * c[i].m * c[i].inv_r2);
            c[i].v += d * c[i].gravity;
        }
#if PROFILE
        prof.end(0);
#endif
        
        ////////////////////////////////
        // broad phase
#if PRIFLE
        prof.start();
#endif
        int num_xs = 0;
        for (int i = 0; i < N; i++) {
            xs[num_xs].first = c[i].pos.x - c[i].r - 0.001;
            xs[num_xs].second.first = i;
            xs[num_xs].second.second = true;
            num_xs++;
            xs[num_xs].first = c[i].pos.x + c[i].r + 0.001;
            xs[num_xs].second.first = i;
            xs[num_xs].second.second = false;
            num_xs++;
        }
        sort(xs, xs + num_xs);
        
        int num_pairs = 0;
        for (int i = 0; i < num_xs; i++) {
            if (xs[i].second.second) {
                int a = xs[i].second.first;
                for (int j = i+1; j < num_xs; j++) {
                    if (xs[j].second.first == a) break;
                    if (xs[j].second.second) {
                        int b = xs[j].second.first;
                        if (abs(c[a].pos.y - c[b].pos.y) < c[a].r + c[b].r + 0.001) {
                            pairs[num_pairs++] = make_pair(a, b);
                        }
                    }
                }
            }
        }
#if PROFILE
        prof.end(1);
#endif
        
        ////////////////////////////////
        // detect collision and solve constraints
#if PROFILE
        prof.start();
#endif
        const double CD = 0.6 / TIME_PER_FRAME;
        for (int k = 0; k < num_pairs; k++) {
            int i = pairs[k].first;
            int j = pairs[k].second;
            if (c[i].is_hover || c[j].is_hover) continue;
            Vec d = c[j].pos - c[i].pos;
            if (d.isSmallerThan(c[i].r + c[j].r + 0.001)) {
                Vec dv = c[j].v - c[i].v;
                double len = d.length();
                Vec norm = d / len;
                double CDD = CD * (c[i].r + c[j].r  + 0.001 - len);
                double C = m_mul_div_sum[i][j] * ((1 + E) * dot(dv, norm) - CDD);
                c[i].v += norm * C * c[i].inv_m;
                c[j].v -= norm * C * c[j].inv_m;
            }
        }
#if PROFILE
        prof.end(2);
#endif
        
        // add friction
        for (int i = 0; i < N; i++) {
            c[i].v *= 0.98;
        }

        ////////////////////////////////
        // integrate
#if PROFILE
        prof.start();
#endif
        for (int i = 0; i < N; i++) {
            c[i].pos += c[i].v * TIME_PER_FRAME;
        }
#if PROFILE
        prof.end(4);
#endif
        
        // kick circles which have overlap
        if (frames % 300 == 0) {
            for (int i = N-1; i >= 0; i--) {
                for (int j = 0; j < i; j++) {
                    if ((c[i].pos - c[j].pos).isSmallerThan(c[i].r + c[j].r)) {
                        int moved = j;
                        if (c[j].m / c[j].r / c[j].r > c[i].m / c[i].r / c[i].r) moved = i;
                        Vec d = c[moved].pos - Vec(0.5, 0.5);
                        d.normalize();
                        c[moved].pos += d * 1;
                    }
                }
            }
        }
        
        if (frames % 100 == 0)
            updateBest();

        
        if (frames % RESTART_FRAME == 0) {
            restart();
        }
        
#if PROFILE
        if (frames % 10000 == 0) {
            prof.report();
        }
#endif
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
    
    const vector<Circle>& current_circles() const {
        return c;
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
        const double MAX_SEC = 9.5;
        //const double MAX_SEC = 1.4;
        clock_t end__ = clock() + MAX_SEC * CLOCKS_PER_SEC;
        
        setup(x, y, r, m);
        while (clock() < end__) for (int i=0; i<100; i++) update();
        
        //cerr << "frames: " << frames << endl;
        
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
            //cerr << "cost: " << cost << endl;
            best_cost = cost;
            for (int i = 0; i < N; i++) {
                best_x[c[i].index] = c[i].pos.x;
                best_y[c[i].index] = c[i].pos.y;
            }
        }
    }
};
