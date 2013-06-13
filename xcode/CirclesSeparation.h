#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include <xmmintrin.h>
using namespace std;

// flags
#define PROFILE 1
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

unsigned int edges_x[MAX_N*2];
unsigned int edges_y[MAX_N*2];
unsigned int LEFT_BIT = (1<<16);
unsigned int EDGE_MASK = (1<<16) - 1;

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
float ball_gravity[MAX_N]__attribute__((aligned(16)));

float dot(float* v1, float* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
}

void normalize(float* v) {
    float len = (float)sqrt(v[0] * v[0] + v[1] * v[1]);
    if (len != 0.0f) {
        v[0] /= len;
        v[1] /= len;
    }
}

// computing buffers and precomputed values
int pairs[MAX_N * 8];
float m_mul_div_sum[MAX_N][MAX_N];

struct PairNode {
    int a, b;
    PairNode *prev, *next;
};
int overlap_cnt[MAX_N][MAX_N];
PairNode pair_nodes[MAX_N][MAX_N];
PairNode* pair_head;

class CirclesSeparation {
private:
    int hover_circle_;
    int frames_;
    int total_frames_;
    int iteration_;
    InputStats input_stats_;
public:
    CirclesSeparation(): hover_circle_(-1), frames_(0), total_frames_(0), iteration_(0)
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
            
            ball_gravity[i] = G * TIME_PER_FRAME * min(4.0, 1 + 0.008 * ball_m[i] * ball_m[i] * c[i].inv_r2);
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
            float priority =  ball_m[i] * c[i].inv_r2;
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
            d[0] = ball_x[index] - 0.5;
            d[1] = ball_y[index] - 0.5;
            normalize(d);
            ball_x[index] += d[0] * (2 + ((float)i / outer.size()) * 2);
            ball_y[index] += d[1] * (2 + ((float)i / outer.size()) * 2);
        }
        
        memset(overlap_cnt, 0, sizeof(overlap_cnt));
        memset(pair_nodes, 0, sizeof(pair_nodes));
        pair_head = NULL;
        for (int i = 0; i < N; i++) {
            edges_x[2*i] = i | LEFT_BIT;
            edges_x[2*i+1] = i;
            edges_y[2*i] = i | LEFT_BIT;
            edges_y[2*i+1] = i;
        }
        sort_edges(edges_x, ball_x, 2*N);
        sort_edges(edges_y, ball_y, 2*N);
    }

    void sort_edges(unsigned int* edges, float* pos, int size) {
        for (int i = 0; i < size; i++) {
            int j = i-1;
            if (j < 0) continue;
            while (j >= 0) {
                int id0 = (edges[j] & EDGE_MASK);
                int id1 = (edges[j+1] & EDGE_MASK);
                float pos0 = pos[id0] + (ball_r[id0] + BOUNCE_MARGIN) * ((edges[j] & LEFT_BIT) ? -1 : 1);
                float pos1 =  pos[id1] + (ball_r[id1] + BOUNCE_MARGIN) * ((edges[j+1] & LEFT_BIT) ? -1 : 1);
                if (pos0 > pos1) {
                    if (!(edges[j] & LEFT_BIT) && (edges[j+1] & LEFT_BIT)) {
                        int& o_cnt = overlap_cnt[min(id0, id1)][max(id0, id1)];
                        o_cnt++;
                        if (o_cnt == 2) {
                            PairNode* node = &pair_nodes[min(id0, id1)][max(id0, id1)];
                            node->a = min(id0, id1);
                            node->b = max(id0, id1);
                            if (pair_head) {
                                pair_head->prev = node;
                            }
                            node->next = pair_head;
                            node->prev = NULL;
                            pair_head = node;
                        }
                    } else if ((edges[j] & LEFT_BIT) && !(edges[j+1] & LEFT_BIT)) {
                        int& o_cnt = overlap_cnt[min(id0, id1)][max(id0, id1)];
                        o_cnt--;
                        if (o_cnt == 1) {
                            PairNode* node = &pair_nodes[min(id0, id1)][max(id0, id1)];
                            if (node->prev) node->prev->next = node->next;
                            if (node->next) node->next->prev = node->prev;
                            if (node == pair_head) pair_head = node->next;
                            node->prev = NULL;
                            node->next = NULL;
                        }
                    }
                    swap(edges[j], edges[j+1]);
                } else {
                    break;
                }
                j--;
            }
        }
    }
    
    void update() {
        frames_++;
        total_frames_++;
        
        __m128 reg[8];
        
        ///////////////////////////////////////////////////////
        // apply force
        PROF_START();
        for (int i = 0; i < N; i+=4) {
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
        sort_edges(edges_x, ball_x, 2*N);
        sort_edges(edges_y, ball_y, 2*N);
        PROF_END(1);
        
        ///////////////////////////////////////////////////////
        // detect collision
        PROF_START()
        int num_pairs = 0;
        PairNode* node = pair_head;
        while (node) {
            int a = node->a;
            int b = node->b;
            node = node->next;
            float dx = ball_x[a] - ball_x[b];
            float dy = ball_y[a] - ball_y[b];
            float sum = ball_r[a] + ball_r[b] + BOUNCE_MARGIN*2;
            if (dx * dx + dy * dy < sum * sum) {
                pairs[num_pairs++] = (a << 16) | b;
            }
        }
        random_shuffle(pairs, pairs + num_pairs);
        PROF_END(2);
        
        ///////////////////////////////////////////////////////
        // solve constraints
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
        PROF_END(3);
        
        // add friction
        for (int i = 0; i < N; i+=16) {
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

            static const float RATE = 1 - DECAY_PER_FRAME;
            reg[0] = _mm_mul_ps(reg[0], _mm_set1_ps(RATE));
            reg[1] = _mm_mul_ps(reg[1], _mm_set1_ps(RATE));
            reg[2] = _mm_mul_ps(reg[2], _mm_set1_ps(RATE));
            reg[3] = _mm_mul_ps(reg[3], _mm_set1_ps(RATE));
            reg[4] = _mm_mul_ps(reg[4], _mm_set1_ps(RATE));
            reg[5] = _mm_mul_ps(reg[5], _mm_set1_ps(RATE));
            reg[6] = _mm_mul_ps(reg[6], _mm_set1_ps(RATE));
            reg[7] = _mm_mul_ps(reg[7], _mm_set1_ps(RATE));

            _mm_store_ps(&ball_vx[i], reg[0]);
            _mm_store_ps(&ball_vy[i], reg[1]);
            _mm_store_ps(&ball_vx[j], reg[2]);
            _mm_store_ps(&ball_vy[j], reg[3]);
            _mm_store_ps(&ball_vx[k], reg[4]);
            _mm_store_ps(&ball_vy[k], reg[5]);
            _mm_store_ps(&ball_vx[l], reg[6]);
            _mm_store_ps(&ball_vy[l], reg[7]);
        }
        
        
        ///////////////////////////////////////////////////////
        // integrate
        PROF_START();
        for (int i = 0; i < N; i+=8) {
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
            cost += (float)sqrt(d[0]*d[0]+d[1]*d[1]) * ball_m[i];
        }
        return cost;
    }
    
    int currentFrame() {
        return frames_;
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

