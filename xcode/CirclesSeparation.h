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
const double E = 0.05;

unsigned int randxor() {
    static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
    unsigned int t;
    t=(x^(x<<11));x=y;y=z;z=w; return( w=(w^(w>>19))^(t^(t>>8)) );
}

struct Circle {
    int index;
    //double sx, sy, x, y, r, m, vx, vy;
    double r, m;
    Vec o_pos, pos, v;
    bool is_hover;
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
        }
        
        // determine initial position
        sort(c.begin(), c.end(), greaterMass);
        for (int i = 0; i < N; i++) {
            /*
            bool hasCollision = false;
            for (int j = 0; j < i; j++) {
                Vec d = c[i].pos - c[j].pos;
                if (d.isSmallerThan(c[i].r + c[j].r)) {
                    hasCollision = true;
                    break;
                }
            }
            while (hasCollision) {
                hasCollision = false;
                c[i].pos.x = -25 + ((double)rand() / RAND_MAX) * 50;
                c[i].pos.y = -25 + ((double)rand() / RAND_MAX) * 50;
                for (int j = 0; j < i; j++) {
                    Vec d = c[i].pos - c[j].pos;
                    if (d.isSmallerThan(c[i].r + c[j].r)) {
                        hasCollision = true;
                        break;
                    }
                }
            }
             */
            double len_from_center = (Vec(0.5, 0.5) - c[i].pos).length();
            double mass = c[i].m;
            double area = c[i].r * c[i].r;
            double jama = area / mass / len_from_center;
            if (jama > 0.1) {
                cerr << i << " is in the way" << endl;
                Vec d = c[i].pos - Vec(0.5, 0.5);
                d.normalize();
                c[i].pos += d * 2;;
            }
        }
    }
    
    void update(double dt) {
        frames++;
        
        // FORCE to the original position
        for (int i = 0; i < N; i++) {
            Vec d = c[i].o_pos - c[i].pos;
            double len = d.length();
            if (len < 1e-5) continue;
            d.normalize();
            c[i].v += d * G * dt * (0.5 + c[i].m);
        }
        
        // FORCE by the collision
        double CD = 0.5 / dt;
        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {
                if (c[i].is_hover || c[j].is_hover) continue;
                Vec d = c[j].pos - c[i].pos;
                if (d.isSmallerThan(c[i].r + c[j].r + 0.003)) {
                    Vec dv = c[j].v - c[i].v;
                    Vec norm = d;
                    norm.normalize();
                    double D = max(0.0, c[i].r + c[j].r  + 0.003 - d.length());
                    double C = c[i].m * c[j].m / (c[i].m + c[j].m) * ((1 + E) * dot(dv, norm) - CD*D);
                    c[i].v += norm * C / (c[i].m);
                    c[j].v -= norm * C / (c[j].m);
                }
            }
        }
    
        // add friction
        for (int i = 0; i < N; i++) {
            c[i].v *= 0.98;
        }
    
        // move
        for (int i = 0; i < N; i++) {
            c[i].pos += c[i].v * dt;
        }
        
        updateBest();
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
                
    void regular_check() {
        for (int i = N-1; i >= 0; i--) {
            for (int j = 0; j < i; j++) {
                if ((c[i].pos - c[j].pos).isSmallerThan(c[i].r + c[j].r)) {
                    Vec d = c[i].pos - Vec(0.5, 0.5);
                    d.normalize();
                    c[i].pos += d * 2;;
                }
            }
        }
        //cerr << "regular check" << endl;
    }
    
    vector<double> minimumWork(vector<double> x_, vector<double> y_, vector<double> r, vector<double> m) {
        setup(x_, y_, r, m);
        
        const double MAX_SEC = 9.0;
        clock_t lap__ = clock();
        clock_t end__ = clock() + MAX_SEC * CLOCKS_PER_SEC;
        
        while (clock() < end__) {
            update(0.01);
            if (clock() - lap__ > CLOCKS_PER_SEC) {
                regular_check();
                lap__ = clock();
            }
        }
        
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
        cerr << "cost: " << cost << endl;
        if (best_cost > cost) {
            best_cost = cost;
            for (int i = 0; i < N; i++) {
                best_x[c[i].index] = c[i].pos.x;
                best_y[c[i].index] = c[i].pos.y;
            }
        }
    }
};

#endif /* defined(__TCO13MM3__CirclesSeparation__) */
