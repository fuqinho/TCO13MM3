//
//  CirclesSeparation.h
//  TCO13MM3
//
//  Created by fuqinho on 6/7/13.
//
//

#ifndef __TCO13MM3__CirclesSeparation__
#define __TCO13MM3__CirclesSeparation__


#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
using namespace std;


const double INF = 1e100;

unsigned int randxor() {
    static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
    unsigned int t;
    t=(x^(x<<11));x=y;y=z;z=w; return( w=(w^(w>>19))^(t^(t>>8)) );
}

struct Circle {
    double sx, sy, x, y, r, m;
};

class CirclesSeparation {
    
private:
    int N;
    double best_cost;
    vector<double> best_x, best_y;
    vector<Circle> circles;
    
public:
    
    
    void setup(vector<double> x, vector<double> y, vector<double> r, vector<double> m) {
        N = x.size();
        circles = vector<Circle>(N);
        best_x = vector<double>(N);
        best_y = vector<double>(N);
        best_cost = INF;
        for (int i = 0; i < x.size(); i++) {
            circles[i].x = circles[i].sx = x[i];
            circles[i].y = circles[i].sy = y[i];
            circles[i].r = r[i];
            circles[i].m = m[i];
        }
        
        // HACK
        setRandomPos();
    }
    
    void update() {
        vector<double> vx(circles.size()), vy(circles.size());
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (isOverlap(i, j)) {
                    double dx = circles[i].x - circles[j].x;
                    double dy = circles[i].y - circles[j].y;
                    if (dx == 0.0 && dy == 0.0) dx += 0.1;
                    double len = sqrt(dx*dx + dy*dy);
                    vx[i] += dx/len;
                    vy[i] += dy/len;
                }
            }
        }
        // normalize
        for (int i = 0; i < N; i++) {
            double len = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
            if (len > 0) {
                vx[i] /= len;
                vy[i] /= len;
            }
        }
        // add weight by m
        for (int i = 0; i < N; i++) {
            vx[i] /= circles[i].m;
            vy[i] /= circles[i].m;
        }
        
        // move
        for (int i = 0; i < N; i++) {
            circles[i].x += vx[i] * 0.002;
            circles[i].y += vy[i] * 0.002;
            
            /*
             x[i] = min(x[i], 100.0);
             x[i] = max(x[i], -100.0);
             y[i] = min(y[i], 100.0);
             y[i] = max(y[i], -100.0);
             */
        }
        
        updateBest();
    }
    
    double currentCost() {
        double cost = 0;
        for (int i = 0; i < N; i++) {
            double dx = circles[i].sx - circles[i].x;
            double dy = circles[i].sy - circles[i].y;
            cost += sqrt(dx*dx + dy*dy) * circles[i].m;
        }
        return cost;
    }
    
    const vector<Circle>& current_circles() const {
        return circles;
    }
    
    
    
    vector<double> minimumWork(vector<double> x_, vector<double> y_, vector<double> r, vector<double> m) {
        setup(x_, y_, r, m);
        
        const double MAX_SEC = 9.0;
        clock_t end__ = clock() + MAX_SEC * CLOCKS_PER_SEC;
        
        while (clock() < end__) update();
        
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
        double dx = circles[i].x - circles[j].x;
        double dy = circles[i].y - circles[j].y;
        return sqrt(dx*dx + dy*dy) < circles[i].r + circles[j].r - 1e-10;
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
        cerr << "cost: " << cost << endl;
        if (best_cost > cost) {
            best_cost = cost;
            for (int i = 0; i < N; i++) {
                best_x[i] = circles[i].x;
                best_y[i] = circles[i].y;
            }
        }
    }
    
    void setRandomPos() {
        for (int i = 0; i < N; i++) {
            circles[i].x = (double)rand() / RAND_MAX;
            circles[i].y = (double)rand() / RAND_MAX;
        }
    }
};

#endif /* defined(__TCO13MM3__CirclesSeparation__) */
