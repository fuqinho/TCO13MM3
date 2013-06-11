#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include <xmmintrin.h>
using namespace std;

// flags
#define PROFILE 0
#define PRINT_FRAMES 1
#define PRINT_SCORE_UPDATES 1

// constant values
const int MAX_N = 512;
const float MAX_RUNNING_TIME = 9.5;
//const float MAX_RUNNING_TIME = 1.4;
const float INF = 1e100;
const float G = 9.8;
const float E = 0.0;
const float DECAY_PER_FRAME = 0.002;
const float TIME_PER_FRAME = 0.0005;
const int RESTART_FRAME = 7500;
const float BOUNCE_MARGIN = 0.0002;

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
    float x, y;
    Vec() {};
    Vec(float x, float y): x(x), y(y) {}
    Vec& operator+=(const Vec& o) { x += o.x; y += o.y; return *this; }
    Vec& operator-=(const Vec& o) { x -= o.x; y -= o.y; return *this; }
    Vec& operator*=(const float m) { x *= m; y *= m; return *this; }
    Vec& operator/=(const float m) { x /= m; y /= m; return *this; }
    void clear() { x=0, y=0; }
    float length() { return sqrt(x*x+y*y); }
    bool isSmallerThan(float r) { return x*x + y*y < r*r; }
    void normalize() { float l=length(); if(l!=0) {x/=l;y/=l;} }
};
//ostream& operator<<(ostream& o,const Vec& v) {o << "(" << v.x << "," << v.y << ")"; return o;}
Vec operator+(const Vec& a, const Vec& b) {return Vec(a.x+b.x,a.y+b.y);}
Vec operator-(const Vec& a, const Vec& b) {return Vec(a.x-b.x,a.y-b.y);}
Vec operator*(const Vec& a, const float m) {return Vec(a.x*m, a.y*m);}
Vec operator/(const Vec& a, const float m) {return Vec(a.x/m, a.y/m);}
float dot(const Vec& a, const Vec& b) {return a.x*b.x+a.y*b.y;}
float cross(const Vec& a, const Vec& b) {return a.x*b.y-a.y*b.x;}
int ccw(const Vec& a, const Vec& b) {float cp=cross(a, b); return cp ? (cp>0?1:-1) : 0;}

struct Circle
{
    // states;
    bool is_hover;

    // properties
    int index;
    float m;
    
    // precomputed values
    float inv_r2;
    float inv_m;
    float gravity;
};

bool greaterMass(const Circle& lhs, const Circle& rhs) {
    return lhs.m > rhs.m;
}

struct InputStats {
    float total_area;
    float max_r;
    InputStats(): total_area(0.0), max_r(0.0) {}
};
    
// state variables
int N;
Circle c[MAX_N];
float best_cost;
double best_x[MAX_N];
double best_y[MAX_N];
double input_r[MAX_N];

float ball_x[MAX_N]__attribute__((aligned(16)));
float ball_y[MAX_N]__attribute__((aligned(16)));
float ball_r[MAX_N]__attribute__((aligned(16)));
float ball_m[MAX_N]__attribute__((aligned(16)));
float ball_vx[MAX_N]__attribute__((aligned(16)));
float ball_vy[MAX_N]__attribute__((aligned(16)));
float ball_ox[MAX_N]__attribute__((aligned(16)));
float ball_oy[MAX_N]__attribute__((aligned(16)));
int ball_index[MAX_N];

float dot(float* v1, float* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
}
void normalize(float* v) {
    float len = sqrt(v[0] * v[0] + v[1] * v[1]);
    if (len != 0.0f) {
        v[0] /= len;
        v[1] /= len;
    }
}

// computing buffers and precomputed values
struct XS {
    float x;
    int index;
    bool is_left;
    bool operator<(const XS& rhs) const {
        return x < rhs.x;
    }
};
XS xs[MAX_N * 2];
int pairs[MAX_N * 8];
float m_mul_div_sum[MAX_N][MAX_N];

class CirclesSeparation {
private:
    Circle* hover_circle_;
    int frames_;
    int total_frames_;
    int iteration_;
    InputStats input_stats_;
public:
    CirclesSeparation(): hover_circle_(NULL), frames_(0), total_frames_(0), iteration_(0)
    {
    }
    
    void setup(vector<double> x, vector<double> y, vector<double> r, vector<double> m) {
        srand(10);
        best_cost = INF;
        
        // read and analyze input data
        N = x.size();
        for (int i = 0; i < N; i++) {
            // save acculate radius
            input_r[i] = r[i];
            
            // properties
            c[i].index = i;
            c[i].m = m[i];
            c[i].is_hover = false;
            
            // precompute values
            c[i].inv_r2 = 1.0 / r[i] / r[i];
            c[i].inv_m = 1.0 / m[i];
            c[i].gravity = G * TIME_PER_FRAME * min(4.0, 1 + 0.008 * c[i].m * c[i].m * c[i].inv_r2);
            
            // data for analysis
            input_stats_.max_r = max(input_stats_.max_r, (float)r[i]);
            input_stats_.total_area += M_PI * r[i] * r[i];
        }
        
        // determine initial position
        sort(c, c + N, greaterMass);

        for (int i = 0; i < N; i++) {
            ball_x[i] = ball_ox[i] = x[c[i].index];
            ball_y[i] = ball_oy[i] = y[c[i].index];
            ball_m[i] = m[c[i].index];
            ball_r[i] = r[c[i].index];
            ball_vx[i] = 0.0f;
            ball_vy[i] = 0.0f;
            ball_index[i] = c[i].index;
        }
        
        // precompute fixed values
        for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
            m_mul_div_sum[i][j] = m_mul_div_sum[j][i] = c[i].m * c[j].m / (c[i].m + c[j].m);
        }
        
        restart();
    }
    
    void restart() {
        // clear states
        frames_ = 0;
        for (int i = 0; i < N; i++) {
            ball_vx[i] = ball_vy[i] = 0.0f;
            ball_x[i] = ball_ox[i];
            ball_y[i] = ball_oy[i];
        }
        
        // determine parameter
        
        // determine initial position
        float thresh = 0.8 + 0.4 * (float)rand() / RAND_MAX;
        
        vector<pair<float, int> > priorities(N);
        for (int i = 0; i < N; i++) {
            float priority =  ball_m[i] * c[i].inv_r2; // TODO efficiency
            priorities[i].first = priority;
            priorities[i].second = i;
        }
        sort(priorities.rbegin(), priorities.rend());
        
        float filledArea = 0.0;
        vector<int> outer;
        for (int i = 0; i < N; i++) {
            int index = priorities[i].second;
            if (filledArea > thresh) {
                outer.push_back(index);
            } else {
                filledArea += ball_r[index] * ball_r[index] * 3.14;
            }
        }
        for (int i = 0; i < outer.size(); i++) {
            int index = outer[i];
            float d[2];
            d[0] = ball_x[index] - 0.5f;
            d[1] = ball_y[index] - 0.5f;
            normalize(d);
            ball_x[index] += d[0] * (2 + ((float)i / outer.size()) * 2);
            ball_y[index] += d[1] * (2 + ((float)i / outer.size()) * 2);
        }
    }
    
    void update() {
        
        frames_++;
        total_frames_++;
        
        ///////////////////////////////////////////////////////
        // apply force
        PROF_START();
        for (int i = 0; i < N; i++) {
            float d[2];
            d[0] = ball_ox[i] - ball_x[i];
            d[1] = ball_oy[i] - ball_y[i];
            normalize(d);
            ball_vx[i] += d[0] * c[i].gravity;
            ball_vy[i] += d[1] * c[i].gravity;
        }
        
        // shake
        static int shake_dir = 0;
        static const float DX[] = {-1, 1, 0, 0};
        static const float DY[] = {0 ,0, -1, 1};
        if (frames_ % 2500 == 0) {
            if (hasOverlap()) {
                expandBalls();
            } else {
                for (int i = 0; i < N; i++) {
                    ball_vx[i] += DX[shake_dir];
                    ball_vy[i] += DY[shake_dir];
                }
                shake_dir = (shake_dir + 1) % 4;
            }
        }
        PROF_END(0);

        ///////////////////////////////////////////////////////
        // broad phase
        PROF_START();
        int num_xs = 0;
        for (int i = 0; i < N; i++) {
            xs[num_xs].x = ball_x[i] - ball_r[i] - BOUNCE_MARGIN;
            xs[num_xs].index = i;
            xs[num_xs].is_left = true;
            num_xs++;
            xs[num_xs].x = ball_x[i] + ball_r[i] + BOUNCE_MARGIN;
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
                    if (xs[j].is_left && abs(ball_y[a] - ball_y[b]) < ball_r[a] + ball_r[b] + BOUNCE_MARGIN*2) {
                        pairs[num_pairs++] = (a << 16) | b;
                    }
                }
            }
        }
        PROF_END(2);

        ///////////////////////////////////////////////////////
        // detect collision and solve constraints
        PROF_START();
        
        float CD = 0.6 / TIME_PER_FRAME;
        if (frames_ < 300) CD = 0.03 / TIME_PER_FRAME;
        for (int k = 0; k < num_pairs; k++) {
            int i = pairs[k] >> 16;
            int j = pairs[k] & ((1<<16)-1);
            if (c[i].is_hover || c[j].is_hover) continue;
            float norm[2];
            norm[0] = ball_x[j] - ball_x[i];
            norm[1] = ball_y[j] - ball_y[i];
            float norm_dist2 = norm[0] * norm[0] + norm[1] * norm[1];
            float minimum_distance = ball_r[i] + ball_r[j] + BOUNCE_MARGIN * 2.0f;

            if (norm_dist2 < minimum_distance * minimum_distance) {
                float dv[2];
                dv[0] = ball_vx[j] - ball_vx[i];
                dv[1] = ball_vy[j] - ball_vy[i];
                float len = sqrt(norm_dist2);
                norm[0] /= len;
                norm[1] /= len;
                float CDD = CD * (minimum_distance - len);
                float C = m_mul_div_sum[i][j] * ((1 + E) * dot(dv, norm) - CDD);
                
                ball_vx[i] += norm[0] * C * c[i].inv_m;
                ball_vy[i] += norm[1] * C * c[i].inv_m;
                ball_vx[j] -= norm[0] * C * c[j].inv_m;
                ball_vy[j] -= norm[1] * C * c[j].inv_m;
            }
        }
        
        PROF_END(3);
        
        // add friction
        for (int i = 0; i < N; i++) {
            ball_vx[i] *= 1 - DECAY_PER_FRAME;
            ball_vy[i] *= 1 - DECAY_PER_FRAME;
        }

        
        ///////////////////////////////////////////////////////
        // integrate
        PROF_START();
        for (int i = 0; i < N; i++) {
            ball_x[i] += ball_vx[i] * TIME_PER_FRAME;
            ball_y[i] += ball_vy[i] * TIME_PER_FRAME;
        }
        PROF_END(4);
        
        ///////////////////////////////////////////////////////
        // others
        
        
        PROF_START();
        if (frames_ % 100 == 0) {
            updateBest();
        }
        if (frames_ >= RESTART_FRAME) {
            iteration_++;
            restart();
        }
        PROF_END(5);

        if (total_frames_ % 10000 == 0) {
            PROF_REPORT();
        }
    }
    
    float currentCost() {
        float cost = 0;
        for (int i = 0; i < N; i++) {
            float d[2];
            d[0] = ball_x[i] - ball_ox[i];
            d[1] = ball_y[i] - ball_oy[i];
            cost += sqrt(d[0]*d[0]+d[1]*d[1]) * ball_m[i];
        }
        return cost;
    }
                
    int currentFrame() {
        return frames_;
    }
        
    void catch_circle(Vec pos) {
        /*
        for (int i = 0; i < N; i++) {
            if ((pos - c[i].pos).isSmallerThan(c[i].r)) {
                c[i].is_hover = true;
                c[i].v.clear();
                hover_circle_ = &c[i];
                break;
            }
        }
         */
    }
    
    void move_circle(Vec pos) {
        if (hover_circle_) {
            //hover_circle_->pos = pos;
        }
    }
    
    void release_circle() {
        if (hover_circle_) {
            hover_circle_->is_hover = false;
            hover_circle_ = NULL;
        }
    }
    
    vector<double> minimumWork(vector<double> x, vector<double> y, vector<double> r, vector<double> m)
    {
        clock_t end__ = clock() + MAX_RUNNING_TIME * CLOCKS_PER_SEC;
        
        setup(x, y, r, m);
        while (clock() < end__) for (int i=0; i<100; i++) update();
        
#if PRINT_FRAMES
        cerr << "total_frames: " << total_frames_ << endl;
#endif
        
        vector<double> res(N*2);
        for (int i = 0; i < N; i++) {
            res[2*i]   = best_x[i];
            res[2*i+1] = best_y[i];
        }
        return res;
    }
    
private:
    
    bool isOverlap(int i, int j) {
        if (i == j) return false;
        double dx = (double)ball_x[i];
        dx -= (double)ball_x[j];
        double dy = (double)ball_y[i];
        dy -= (double)ball_y[j];
        double rsum = input_r[ball_index[i]] + input_r[ball_index[j]];
        return dx*dx + dy*dy < rsum*rsum;
    }
    
    
    bool hasOverlap() {
        for (int i = 0; i < N; i++)
            for (int j = i+1; j < N; j++)
                if (isOverlap(i, j))
                    return true;
        return false;
    }
    
    void expandBalls() {
        for (int i = 0; i < N; i++) {
            ball_x[i] *= 1.1;
            ball_y[i] *= 1.1;
        }
    }
    
    void updateBest() {
        if (hasOverlap()) return;
        
        float cost = currentCost();
        if (best_cost > cost) {
#if PRINT_SCORE_UPDATES
            cerr << "cost: " << cost << endl;
#endif
            best_cost = cost;
            for (int i = 0; i < N; i++) {
                best_x[ball_index[i]] = ball_x[i];
                best_y[ball_index[i]] = ball_y[i];
            }
        }
    }
};
