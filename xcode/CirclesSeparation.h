////////////////////////////////////////////////////////////////////////////////////
// Flags
////////////////////////////////////////////////////////////////////////////////////

#define PROFILE 1
#define PRINT_SIMULATED_FRAMES 1
#define PRINT_SCORE_UPDATES 1
#define PRINT_BEST_PARAMETERS 1
#define PRINT_BEST_PARAMETERS_PERIODICALLY 1

#define ENABLE_RESTART 1
#define ENABLE_SHAKE 1


////////////////////////////////////////////////////////////////////////////////////
// Constant values
////////////////////////////////////////////////////////////////////////////////////

const int MAX_N = 512;
const float MAX_RUNNING_TIME = 9.5;
//const float MAX_RUNNING_TIME = 2.7;
const float INF = 1e100;
const float G = 9.8;
const float E = 0.0;
const float DECAY_PER_FRAME = 0.002;
const float TIME_PER_FRAME = 0.0005;
const float BOUNCE_MARGIN = 0.0001;
const float ANTI_PENETRATION_COEFFICIENT = 0.8;
const float ANTI_PENETRATION_MILD_COEFFICIENT = 0.03;
const float ANTI_PENETRATION_MILD_FRAMES = 300;
const int FRAMES_PER_PERIOD = 2500;
const int MIN_PERIODS = 2;
const int MAX_PERIODS = 11;
const float RESTART_THRESHOLD_VS_PERIOD = 1.001;
const float RESTART_THRESHOLD_VS_ITERATION = 1.01;


////////////////////////////////////////////////////////////////////////////////////
// Headers
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include <xmmintrin.h>
using namespace std;


////////////////////////////////////////////////////////////////////////////////////
// Auxiliary classes
////////////////////////////////////////////////////////////////////////////////////

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

struct Circle
{
    bool is_hover;
    int index;
    float m;
    float inv_r2;
    float inv_m;
};

struct Candidate {
    int a, b;
    Candidate *prev, *next;
};

struct InputStats {
    float total_area;
    float max_r;
    InputStats(): total_area(0.0), max_r(0.0) {}
};

struct Parameter {
    float initial_fill_limit;
};


////////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////////

int N;
Circle c[MAX_N];
float best_cost;
double best_result[MAX_N*2];
double input_r[MAX_N];
double input_m[MAX_N];

float ball_x[MAX_N]__attribute__((aligned(16)));
float ball_y[MAX_N]__attribute__((aligned(16)));
float ball_r[MAX_N]__attribute__((aligned(16)));
float ball_m[MAX_N]__attribute__((aligned(16)));
float ball_vx[MAX_N]__attribute__((aligned(16)));
float ball_vy[MAX_N]__attribute__((aligned(16)));
float ball_ox[MAX_N]__attribute__((aligned(16)));
float ball_oy[MAX_N]__attribute__((aligned(16)));
int ball_index[MAX_N]__attribute__((aligned(16)));;
float ball_gravity[MAX_N]__attribute__((aligned(16)));

unsigned int bounds_x[MAX_N*2];
unsigned int bounds_y[MAX_N*2];
unsigned int BOUNDS_ISLEFT_BIT = (1<<16);
unsigned int BOUNDS_NODEID_MASK = (1<<16) - 1;

int collisions[MAX_N * 8];
int num_collisions;

int overlapping_count[MAX_N][MAX_N];
Candidate candidates[MAX_N][MAX_N];
Candidate* candidates_list_head;

float m_mul_div_sum[MAX_N][MAX_N];
float period_best[MAX_N];


////////////////////////////////////////////////////////////////////////////////////
// Auxiilary functions
////////////////////////////////////////////////////////////////////////////////////

static inline unsigned int randxor() {
    static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
    unsigned int t;
    t=(x^(x<<11));x=y;y=z;z=w; return( w=(w^(w>>19))^(t^(t>>8)) );
}

static inline int myrandom (int i) {
    return randxor() % i;
}

static inline bool greaterMass(const Circle& lhs, const Circle& rhs) {
    return lhs.m > rhs.m;
}

static inline float dot(float* v1, float* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
}

static inline void normalize(float* v) {
    float len = (float)sqrt(v[0] * v[0] + v[1] * v[1]);
    if (len != 0.0f) {
        v[0] /= len;
        v[1] /= len;
    }
}

static inline float overlapping_area(float x1, float y1, float r1, float x2, float y2, float r2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float rsum = r1 + r2;
    float dist_squared = dx*dx + dy*dy;
    if (dist_squared >= rsum*rsum) {
        return 0; // no overlap
    }
    float a = r2, b = r1, c = sqrt(dist_squared);
    if (c + min(r1, r2) <= max(r1, r2)) {
        return (float)M_PI * min(r1, r2) * min(r1, r2); // contained
    }
    float cosA = (b*b + c*c - a*a) / (2.0f * b * c);
    float cosB = (a*a + c*c - b*b) / (2.0f * a * c);
    float fanA = r1 * r1 * acos(cosA);
    float fanB = r2 * r2 * acos(cosB);
    float triangle = c * b * (float)sqrt(1.0 - cosA*cosA);
    return fanA + fanB - triangle;
}


////////////////////////////////////////////////////////////////////////////////////
// Main class
////////////////////////////////////////////////////////////////////////////////////

class CirclesSeparation {
public:
    CirclesSeparation(): hover_circle_(-1), frames_(0), total_frames_(0), periods_(0), frames_in_period_(0), iteration_(0) {}
    
    vector<double> minimumWork(vector<double> x, vector<double> y, vector<double> r, vector<double> m) {
        const clock_t end__ = clock() + MAX_RUNNING_TIME * CLOCKS_PER_SEC;
        
        setup(x, y, r, m);
        while (clock() < end__) for (int i=0; i<100; i++) update();
#if PRINT_SIMULATED_FRAMES
        cerr << "total_frames: " << total_frames_ << endl;
#endif
        printBestParams();
        
        return vector<double>(best_result, best_result + N * 2);
    }
    
    void setup(const vector<double>& x, const vector<double>& y, const vector<double>& r, const vector<double>& m) {
        initializeOnce();
        readInput(x, y, r, m);
        analyzeInput();
        determineParameter();
        precomputeValues();
        start();
    }
    
    void update() {
        ///////////////////////////////////////////////////////
        // apply force
        PROF_START();
        applyForceToCenter();
        applyExternalForce(frames_);
        PROF_END(0);
        
        ///////////////////////////////////////////////////////
        // broad phase (find pairs of balls which may collide)
        PROF_START();
        sortBounds(bounds_x, ball_x, 2*N);
        sortBounds(bounds_y, ball_y, 2*N);
        PROF_END(1);
        
        ///////////////////////////////////////////////////////
        // detect collisions strictly
        PROF_START();
        detectCollisions();
        PROF_END(2);
        
        ///////////////////////////////////////////////////////
        // change velocity
        PROF_START();
        solveConstraints(frames_ < ANTI_PENETRATION_MILD_FRAMES ?
                         ANTI_PENETRATION_MILD_COEFFICIENT : ANTI_PENETRATION_COEFFICIENT);
        decelerateVelocity(1 - DECAY_PER_FRAME);
        PROF_END(3);
        
        ///////////////////////////////////////////////////////
        // change position
        PROF_START();
        integrateVelocity();
        PROF_END(4);
        
        ///////////////////////////////////////////////////////
        // others
        PROF_START();
        frames_++;
        frames_in_period_++;
        total_frames_++;
        if (frames_ % 100 == 0) updateBest();
        if (frames_in_period_ == FRAMES_PER_PERIOD) mayRestart();
        PROF_END(5);
        
        if (total_frames_ % 10000 == 0) {
            PROF_REPORT();
#if PRINT_BEST_PARAMETERS_PERIODICALLY
            printBestParams();
#endif
        }
    }
    
    float currentCost() {
        float cost = 0;
        for (int i = 0; i < N; i++) {
            float d[2];
            d[0] = ball_x[i] - ball_ox[i];
            d[1] = ball_y[i] - ball_oy[i];
            cost += (float)sqrt(d[0]*d[0]+d[1]*d[1]) * ball_m[i];
        }
        return cost;
    }
    
    int currentFrame() {
        return frames_;
    }
    
    int currentIteration() {
        return iteration_;
    }
    
    int currentPeriod() {
        return periods_;
    }
    
    void catch_circle(float x, float y) {
        for (int i = 0; i < N; i++) {
            float d[2];
            d[0] = x - ball_x[i];
            d[1] = y - ball_y[i];
            if (d[0]*d[0] + d[1]*d[1] < ball_r[i]*ball_r[i]) {
                c[i].is_hover = true;
                ball_vx[i] = 0;
                ball_vy[i] = 0;
                hover_circle_ = i;
                break;
            }
        }
    }
    
    void move_circle(float x, float y) {
        if (hover_circle_ >= 0) {
            ball_x[hover_circle_] = x;
            ball_y[hover_circle_] = y;
        }
    }
    
    void release_circle() {
        if (hover_circle_ >= 0) {
            c[hover_circle_].is_hover = false;
            hover_circle_ = -1;
        }
    }
    

private:
    
    int frames_;
    int total_frames_;
    int periods_;
    int frames_in_period_;
    int iteration_;
    float period_best_cost, prev_period_best_cost;
    InputStats input_stats_;
    Parameter param_, best_param_;
    int hover_circle_;
    
    void initializeOnce() {
        srand(10);
        best_cost = INF;
        for (int i = 0; i < MAX_PERIODS; i++) period_best[i] = INF;
    }
    
    void readInput(const vector<double>& x, const vector<double>& y, const vector<double>& r, const vector<double>& m) {
        // read and analyze input data
        N = x.size();
        for (int i = 0; i < N; i++) {
            // save acculate data
            input_r[i] = r[i];
            input_m[i] = m[i];
            
            // properties
            c[i].index = i;
            c[i].m = m[i];
            c[i].is_hover = false;
        }
        
        // order circles. (heavier is prior)
        sort(c, c + N, greaterMass);
        
        for (int i = 0; i < N; i++) {
            ball_x[i] = ball_ox[i] = x[c[i].index];
            ball_y[i] = ball_oy[i] = y[c[i].index];
            ball_m[i] = m[c[i].index];
            ball_r[i] = r[c[i].index];
            ball_index[i] = c[i].index;
        }
    }
    
    void analyzeInput() {
        for (int i = 0; i < N; i++) {
            input_stats_.max_r = max(input_stats_.max_r, ball_r[i]);
            input_stats_.total_area += M_PI * ball_r[i] * ball_r[i];
        }
    }
    
    void determineParameter() {
        
    }
    
    void precomputeValues() {
        // precompute for balls
        for (int i = 0; i < N; i++) {
            c[i].inv_r2 = 1.0 / input_r[c[i].index] / input_r[c[i].index];
            c[i].inv_m = 1.0 / input_m[c[i].index];
            
            ball_gravity[i] = G * TIME_PER_FRAME * min(3.0, 1 + 0.003 * ball_m[i] * c[i].inv_r2);
        }
        
        // precompute for pairs
        for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
            m_mul_div_sum[i][j] = m_mul_div_sum[j][i] = c[i].m * c[j].m / (c[i].m + c[j].m);
        }
    }
    
    void clearState() {
    }
    
    void adjustParameter(int iteration) {
        rand();
        param_.initial_fill_limit = 0.6 + 0.6 * (float)rand() / RAND_MAX;
    }
    
    void start() {
        adjustParameter(iteration_);
        initIteration();
    }
    
    void initIteration() {
        frames_ = 0;
        periods_ = 0;
        period_best_cost = prev_period_best_cost = INF;

        // clear state of balls and collisions
        for (int i = 0; i < N; i++) {
            ball_vx[i] = ball_vy[i] = 0.0f;
            ball_x[i] = ball_ox[i];
            ball_y[i] = ball_oy[i];
        }
        memset(overlapping_count, 0, sizeof(overlapping_count));
        memset(candidates, 0, sizeof(candidates));
        for (int i = 0; i < N; i++) for (int j = i; j < N; j++) {
            candidates[i][j].a = i;
            candidates[i][j].b = j;
        }
        candidates_list_head = NULL;
        for (int i = 0; i < N; i++) {
            bounds_x[2*i] = i | BOUNDS_ISLEFT_BIT;
            bounds_x[2*i+1] = i;
            bounds_y[2*i] = i | BOUNDS_ISLEFT_BIT;
            bounds_y[2*i+1] = i;
        }

        // set initial position for balls
        vector<pair<float, int> > priorities(N);
        for (int i = 0; i < N; i++) {
            float priority =  ball_m[i] * c[i].inv_r2;
            priorities[i].first = priority;
            priorities[i].second = i;
        }
        sort(priorities.rbegin(), priorities.rend());
        
        float filledArea = 0.0;
        vector<int> outer;
        vector<int> inner;
        for (int i = 0; i < N; i++) {
            int index = priorities[i].second;
            if (filledArea > param_.initial_fill_limit) {
                outer.push_back(index);
            } else {
                inner.push_back(index);
                filledArea += ball_r[index] * ball_r[index] * 3.14;
            }
        }
        
        // find weighted center
        float xsum = 0, ysum = 0, msum = 0;
        for (int i = 0; i < N; i++) {
            msum += ball_m[i];
            xsum += ball_m[i] * ball_x[i];
            ysum += ball_m[i] * ball_y[i];
        }
        float center_x = xsum / msum;
        float center_y = ysum / msum;
        
        for (int i = 0; i < outer.size(); i++) {
            int index = outer[i];
            bool canInsert = false;
            if (i < outer.size() / 3) {
                double area_sum = ball_r[index] * ball_r[index] * 3.14f;;
                double search_r = input_stats_.max_r;
                for (int j = 0; j < inner.size(); j++) {
                    int id2 = inner[j];
                    area_sum += overlapping_area(ball_x[index], ball_y[index], search_r,
                                                 ball_x[id2], ball_y[id2], ball_r[id2]);
                }
                if (area_sum < search_r * search_r * 3.14f * 0.9f) {
                    canInsert = true;
                    inner.push_back(index);
                }
            }
            if (!canInsert) {
                float d[2];
                d[0] = ball_x[index] - center_x;
                d[1] = ball_y[index] - center_y;
                normalize(d);
                ball_x[index] += d[0] * (2 + ((float)i / outer.size()) * 3);
                ball_y[index] += d[1] * (2 + ((float)i / outer.size()) * 3);
            }
        }
        sortBounds(bounds_x, ball_x, 2*N);
        sortBounds(bounds_y, ball_y, 2*N);
        
        initPeriod();
    }

    void initPeriod() {
        frames_in_period_ = 0;
        prev_period_best_cost = period_best_cost;
        period_best_cost = INF;
    }
    
    bool shouldRestart() {
        if (periods_ >= MAX_PERIODS) return true;
        if (periods_ < MIN_PERIODS) return false;
        
        bool res = false;
        if (period_best_cost == INF || period_best_cost > period_best[periods_-1] * RESTART_THRESHOLD_VS_ITERATION) {
            res = true;
        }
        if (periods_ >= 2 && period_best_cost > prev_period_best_cost * RESTART_THRESHOLD_VS_PERIOD) {
            res = true;
        }
        period_best[periods_-1] = min(period_best[periods_-1], period_best_cost);
        return res;
    }
    
    void mayRestart() {
        periods_++;
        if (shouldRestart())
            restart();
        else
            initPeriod();
    }
    
    void restart() {
#if ENABLE_RESTART
        iteration_++;
        start();
#endif
    }
    
    void applyForceToCenter() {
        __m128 reg[8];
        for (int i = 0; i < N; i+=4) {
            // SAME AS FOLLOWING CODE
            //
            // float d[2];
            // d[0] = ball_ox[i] - ball_x[i];
            // d[1] = ball_oy[i] - ball_y[i];
            // normalize(d);
            // ball_vx[i] += d[0] * c[i].gravity;
            // ball_vy[i] += d[1] * c[i].gravity;
            
            reg[0] = _mm_load_ps(&ball_ox[i]);
            reg[1] = _mm_load_ps(&ball_x[i]);
            reg[2] = _mm_load_ps(&ball_oy[i]);
            reg[3] = _mm_load_ps(&ball_y[i]);
            reg[4] = _mm_load_ps(&ball_vx[i]);
            reg[5] = _mm_load_ps(&ball_vy[i]);
            reg[6] = _mm_load_ps(&ball_gravity[i]);
            
            reg[0] = _mm_sub_ps(reg[0], reg[1]);
            reg[2] = _mm_sub_ps(reg[2], reg[3]);
            reg[1] = _mm_mul_ps(reg[0], reg[0]);
            reg[3] = _mm_mul_ps(reg[2], reg[2]);
            reg[1] = _mm_add_ps(reg[1], reg[3]);
            reg[1] = _mm_add_ps(reg[1], _mm_set1_ps(1e-10f));
            reg[1] = _mm_rsqrt_ps(reg[1]);
            reg[0] = _mm_mul_ps(reg[0], reg[1]);
            reg[0] = _mm_mul_ps(reg[0], reg[6]);
            reg[2] = _mm_mul_ps(reg[2], reg[1]);
            reg[2] = _mm_mul_ps(reg[2], reg[6]);
            reg[4] = _mm_add_ps(reg[4], reg[0]);
            reg[5] = _mm_add_ps(reg[5], reg[2]);
            
            _mm_store_ps(&ball_vx[i], reg[4]);
            _mm_store_ps(&ball_vy[i], reg[5]);
        }
    }
    
    void applyExternalForce(int frame) {
#if ENABLE_SHAKE
        if (frame != 0 && frame % FRAMES_PER_PERIOD == 0) {
            if (hasOverlap()) {
                expandBalls();
            } else {
                shakeBalls();
            }
        }
        if (periods_ == 0 && frame % 500 == 0)
            shakeBalls();
#endif
    }
    
    void detectCollisions() {
        num_collisions = 0;
        Candidate* node = candidates_list_head;
        while (node) {
            int a = node->a;
            int b = node->b;
            node = node->next;
            float dx = ball_x[a] - ball_x[b];
            float dy = ball_y[a] - ball_y[b];
            float sum = ball_r[a] + ball_r[b] + BOUNCE_MARGIN*2;
            if (dx * dx + dy * dy < sum * sum) {
                collisions[num_collisions++] = (a << 16) | b;
            }
        }
        random_shuffle(collisions, collisions + num_collisions, myrandom);
    }
    
    void solveConstraints(float coefficient) {
        float CD = coefficient / TIME_PER_FRAME;
        if (frames_ < ANTI_PENETRATION_MILD_FRAMES)
            CD = ANTI_PENETRATION_MILD_COEFFICIENT / TIME_PER_FRAME;
        
        for (int k = 0; k < num_collisions; k++) {
            int i = collisions[k] >> 16;
            int j = collisions[k] & ((1<<16)-1);
            if (c[i].is_hover || c[j].is_hover) continue;
            float norm[2];
            norm[0] = ball_x[j] - ball_x[i];
            norm[1] = ball_y[j] - ball_y[i];
            float norm_dist2 = norm[0] * norm[0] + norm[1] * norm[1];
            float minimum_distance = ball_r[i] + ball_r[j] + BOUNCE_MARGIN*2;
            
            float dv[2];
            dv[0] = ball_vx[j] - ball_vx[i];
            dv[1] = ball_vy[j] - ball_vy[i];
            float len = (float)sqrt(norm_dist2);
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
    
    void decelerateVelocity(float rate) {
        __m128 reg[8];
        for (int i = 0; i < N; i += 16) {
            // SAME AS FOLLOWING CODE
            //
            // ball_vx[i] *= 1 - DECAY_PER_FRAME;
            // ball_vy[i] *= 1 - DECAY_PER_FRAME;
            
            int j = i+4, k = j+4, l = k+4;
            reg[0] = _mm_load_ps(&ball_vx[i]);
            reg[1] = _mm_load_ps(&ball_vy[i]);
            reg[2] = _mm_load_ps(&ball_vx[j]);
            reg[3] = _mm_load_ps(&ball_vy[j]);
            reg[4] = _mm_load_ps(&ball_vx[k]);
            reg[5] = _mm_load_ps(&ball_vy[k]);
            reg[6] = _mm_load_ps(&ball_vx[l]);
            reg[7] = _mm_load_ps(&ball_vy[l]);
            
            reg[0] = _mm_mul_ps(reg[0], _mm_set1_ps(rate));
            reg[1] = _mm_mul_ps(reg[1], _mm_set1_ps(rate));
            reg[2] = _mm_mul_ps(reg[2], _mm_set1_ps(rate));
            reg[3] = _mm_mul_ps(reg[3], _mm_set1_ps(rate));
            reg[4] = _mm_mul_ps(reg[4], _mm_set1_ps(rate));
            reg[5] = _mm_mul_ps(reg[5], _mm_set1_ps(rate));
            reg[6] = _mm_mul_ps(reg[6], _mm_set1_ps(rate));
            reg[7] = _mm_mul_ps(reg[7], _mm_set1_ps(rate));
            
            _mm_store_ps(&ball_vx[i], reg[0]);
            _mm_store_ps(&ball_vy[i], reg[1]);
            _mm_store_ps(&ball_vx[j], reg[2]);
            _mm_store_ps(&ball_vy[j], reg[3]);
            _mm_store_ps(&ball_vx[k], reg[4]);
            _mm_store_ps(&ball_vy[k], reg[5]);
            _mm_store_ps(&ball_vx[l], reg[6]);
            _mm_store_ps(&ball_vy[l], reg[7]);
        }
    }
    
    void integrateVelocity() {
        __m128 reg[8];
        for (int i = 0; i < N; i+=8) {
            // SAME AS FOLLOWING CODE
            //
            // ball_x[i] += ball_vx[i] * TIME_PER_FRAME;
            // ball_y[i] += ball_vy[i] * TIME_PER_FRAME;
            
            int j = i+4;
            reg[0] = _mm_load_ps(&ball_vx[i]);
            reg[1] = _mm_load_ps(&ball_vy[i]);
            reg[2] = _mm_load_ps(&ball_x[i]);
            reg[3] = _mm_load_ps(&ball_y[i]);
            reg[4] = _mm_load_ps(&ball_vx[j]);
            reg[5] = _mm_load_ps(&ball_vy[j]);
            reg[6] = _mm_load_ps(&ball_x[j]);
            reg[7] = _mm_load_ps(&ball_y[j]);
            
            reg[0] = _mm_mul_ps(reg[0], _mm_set1_ps(TIME_PER_FRAME));
            reg[1] = _mm_mul_ps(reg[1], _mm_set1_ps(TIME_PER_FRAME));
            reg[2] = _mm_add_ps(reg[2], reg[0]);
            reg[3] = _mm_add_ps(reg[3], reg[1]);
            reg[4] = _mm_mul_ps(reg[4], _mm_set1_ps(TIME_PER_FRAME));
            reg[5] = _mm_mul_ps(reg[5], _mm_set1_ps(TIME_PER_FRAME));
            reg[6] = _mm_add_ps(reg[6], reg[4]);
            reg[7] = _mm_add_ps(reg[7], reg[5]);
            
            _mm_store_ps(&ball_x[i], reg[2]);
            _mm_store_ps(&ball_y[i], reg[3]);
            _mm_store_ps(&ball_x[j], reg[6]);
            _mm_store_ps(&ball_y[j], reg[7]);
        }
    }
    
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
    
    void shakeBalls() {
        static int shake_dir = 0;
        static const float DX[] = {-1, 1, 0, 0};
        static const float DY[] = {0 ,0, -1, 1};
        
        for (int i = 0; i < N; i++) {
            ball_vx[i] += DX[shake_dir];
            ball_vy[i] += DY[shake_dir];
        }
        shake_dir = (shake_dir + 1) % 4;
    }
    
    void updateBest() {
        if (hasOverlap()) return;
        
        float cost = currentCost();
        period_best_cost = min(period_best_cost, cost);
        if (best_cost > cost) {
#if PRINT_SCORE_UPDATES
            cerr << "cost: " << cost << endl;
#endif
            best_cost = cost;
            for (int i = 0; i < N; i++) {
                best_result[ball_index[i] * 2] = ball_x[i];
                best_result[ball_index[i] * 2 + 1] = ball_y[i];
            }
            best_param_ = param_;
        }
    }
    
    void printBestParams() {
#if PRINT_BEST_PARAMETERS
        cerr << N << '\t'
            << input_stats_.total_area << '\t' << input_stats_.max_r << '\t'
            << best_param_.initial_fill_limit << endl;
#endif
    }
    
    void pushPairToList(int a, int b) {
        Candidate* node = &candidates[a][b];
        node->a = a;
        node->b = b;
        if (candidates_list_head) {
            candidates_list_head->prev = node;
        }
        node->next = candidates_list_head;
        node->prev = NULL;
        candidates_list_head = node;
    }
    
    void removePairFromList(int a, int b) {
        Candidate* node = &candidates[a][b];
        if (node->prev) node->prev->next = node->next;
        if (node->next) node->next->prev = node->prev;
        if (node == candidates_list_head) candidates_list_head = node->next;
        node->prev = NULL;
        node->next = NULL;
    }
    
    void sortBounds(unsigned int* bounds, float* pos, int size) {
        for (int i = 0; i < size; i++) {
            int j = i-1;
            if (j < 0) continue;
            while (j >= 0) {
                int id0 = (bounds[j] & BOUNDS_NODEID_MASK);
                int id1 = (bounds[j+1] & BOUNDS_NODEID_MASK);
                float pos0 = pos[id0] + (ball_r[id0] + BOUNCE_MARGIN) * ((bounds[j] & BOUNDS_ISLEFT_BIT) ? -1 : 1);
                float pos1 =  pos[id1] + (ball_r[id1] + BOUNCE_MARGIN) * ((bounds[j+1] & BOUNDS_ISLEFT_BIT) ? -1 : 1);
                if (pos0 > pos1) {
                    int a = min(id0, id1), b = max(id0, id1);
                    if (!(bounds[j] & BOUNDS_ISLEFT_BIT) && (bounds[j+1] & BOUNDS_ISLEFT_BIT)) {
                        overlapping_count[a][b]++;
                        if (overlapping_count[a][b] == 2) pushPairToList(a, b);
                    } else if ((bounds[j] & BOUNDS_ISLEFT_BIT) && !(bounds[j+1] & BOUNDS_ISLEFT_BIT)) {
                        overlapping_count[a][b]--;
                        if (overlapping_count[a][b] == 1) removePairFromList(a, b);
                    }
                    swap(bounds[j], bounds[j+1]);
                } else {
                    break;
                }
                j--;
            }
        }
    }
};

