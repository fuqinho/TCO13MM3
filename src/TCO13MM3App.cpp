#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/TextureFont.h"
#include "CirclesSeparation.h"
#include "TestCases.h"
#include <sstream>

using namespace ci;
using namespace ci::app;
using namespace std;



const int SCREEN_W = 800;
const int SCREEN_H = 600;
const int MAX_N = 600;
double csx[MAX_N];
double csy[MAX_N];
double cx[MAX_N];
double cy[MAX_N];
double cr[MAX_N];

class TCO13MM3App : public AppNative {
public:
    void prepareSettings( Settings *settings );
	void setup();
	void mouseDown( MouseEvent event );
    void mouseUp( MouseEvent event );
    void mouseDrag( MouseEvent event );
	void update();
	void draw();
private:
    Font mFont;
    gl::TextureFontRef mTextureFont;
    
    CirclesSeparation *solution;
    double prj_left;
    double prj_top;
    double prj_ratio;
};

void TCO13MM3App::prepareSettings( Settings *settings )
{
    settings->enableHighDensityDisplay();
    settings->setWindowSize( SCREEN_W, SCREEN_H );
}

void TCO13MM3App::setup()
{
    // create font
    mFont = Font( "Arial", 22 );
    mTextureFont = gl::TextureFont::create( mFont );
    
    solution = new CirclesSeparation;
    solution->setup(vector<double>(XX, XX + sizeof(XX)/sizeof(double)),
                    vector<double>(YY, YY + sizeof(YY)/sizeof(double)),
                    vector<double>(RR, RR + sizeof(RR)/sizeof(double)),
                    vector<double>(MM, MM + sizeof(MM)/sizeof(double)));
}

void TCO13MM3App::mouseDown( MouseEvent event )
{
    solution->catch_circle(Vec(prj_left + event.getX() / prj_ratio,
                               prj_top + event.getY() / prj_ratio));
}

void TCO13MM3App::mouseUp( MouseEvent event )
{
    solution->release_circle();
}

void TCO13MM3App::mouseDrag( MouseEvent event )
{
    solution->move_circle(Vec(prj_left + event.getX() / prj_ratio,
                              prj_top + event.getY() / prj_ratio));
}

void TCO13MM3App::update()
{
    solution->update(0.003);
}

void TCO13MM3App::draw()
{
	// clear out the window with black
	gl::clear( Color( 1.0f, 1.0f, 1.0f ) );
    
    const vector<Circle>& circles = solution->current_circles();
    
    // calculate projection ratio
    double left = -0.5, right = 1.5, top = -0.5, bottom = 1.5;
    for (int i = 0; i < circles.size(); i++) {
        left = min(left, circles[i].pos.x - circles[i].r);
        right = max(right, circles[i].pos.x + circles[i].r);
        top = min(top, circles[i].pos.y - circles[i].r);
        bottom = max(bottom, circles[i].pos.y + circles[i].r);
    }
    double rate_w = SCREEN_W / (right - left);
    double rate_h = SCREEN_H / (bottom - top);
    double rate = min(rate_w, rate_h);
    prj_left = left;
    prj_top = top;
    prj_ratio = rate;
    
    // precompute position on the screen
    for (int i = 0; i < circles.size(); i++) {
        cx[i] = (circles[i].pos.x - left) * rate;
        cy[i] = (circles[i].pos.y - top) * rate;
        cr[i] = circles[i].r * rate;
        csx[i] = (circles[i].o_pos.x - left) * rate;
        csy[i] = (circles[i].o_pos.y - top) * rate;
    }
    
    // draw inner circle
    gl::enableAlphaBlending();
    for (int i = 0; i < circles.size(); i++) {
        gl::color( ColorA ( 0.0f, 0.0f, 1.0, circles[i].m ) );
        gl::drawSolidCircle(Vec2f(cx[i], cy[i]), cr[i]);
    }
    
    // draw circle edge
    gl::color( Color( 0.0f, 0.0f, 0.0f ) );
    for (int i = 0; i < circles.size(); i++)
        gl::drawStrokedCircle(Vec2f(cx[i], cy[i]), cr[i]);
    
    // draw lines
    for (int i = 0; i < circles.size(); i++) {
        gl::lineWidth(8.0 * circles[i].m);
        gl::color( ColorA( 1.0, 0.0, 0.0, circles[i].m ) );
        gl::drawLine(Vec2f(cx[i], cy[i]), Vec2f(csx[i], csy[i]));
    }
    
    // draw current score
    gl::color( Color( 0.0f, 0.0f, 0.0f ) );
    std::stringstream ss;
    ss << "cost: " << solution->currentCost() << "  " << solution->currentFrame();
    mTextureFont->drawString(ss.str(),  Vec2f(0, 20) );
}

CINDER_APP_NATIVE( TCO13MM3App, RendererGl )
