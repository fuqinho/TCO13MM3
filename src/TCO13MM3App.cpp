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
float csx[512];
float csy[512];
float cx[512];
float cy[512];
float cr[512];

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
    float prj_left;
    float prj_top;
    float prj_ratio;
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
    for (int i=0; i<20; i++)
    solution->update();
}

void TCO13MM3App::draw()
{
	// clear out the window with black
	gl::clear( Color( 1.0f, 1.0f, 1.0f ) );
    
    // calculate projection ratio
    float left = -0.5, right = 1.5, top = -0.5, bottom = 1.5;
    for (int i = 0; i < N; i++) {
        left = min(left, c[i].pos.x - c[i].r);
        right = max(right, c[i].pos.x + c[i].r);
        top = min(top, c[i].pos.y - c[i].r);
        bottom = max(bottom, c[i].pos.y + c[i].r);
    }
    float rate_w = SCREEN_W / (right - left);
    float rate_h = SCREEN_H / (bottom - top);
    float rate = min(rate_w, rate_h);
    prj_left = left;
    prj_top = top;
    prj_ratio = rate;
    
    // precompute position on the screen
    for (int i = 0; i < N; i++) {
        cx[i] = (c[i].pos.x - left) * rate;
        cy[i] = (c[i].pos.y - top) * rate;
        cr[i] = c[i].r * rate;
        csx[i] = (c[i].o_pos.x - left) * rate;
        csy[i] = (c[i].o_pos.y - top) * rate;
    }
    
    // draw inner circle
    gl::enableAlphaBlending();
    for (int i = 0; i < N; i++) {
        gl::color( ColorA ( 0.0f, 0.0f, 1.0f, c[i].m ) );
        gl::drawSolidCircle(Vec2f(cx[i], cy[i]), cr[i]);
    }
    
    // draw circle edge
    gl::color( Color( 0.0f, 0.0f, 0.0f ) );
    for (int i = 0; i < N; i++)
        gl::drawStrokedCircle(Vec2f(cx[i], cy[i]), cr[i]);
    
    // draw lines
    for (int i = 0; i < N; i++) {
        gl::lineWidth(8.0 * c[i].m);
        gl::color( ColorA( 1.0f, 0.0f, 0.0f, c[i].m ) );
        gl::drawLine(Vec2f(cx[i], cy[i]), Vec2f(csx[i], csy[i]));
    }
    
    // draw current score
    gl::color( Color( 0.0f, 0.0f, 0.0f ) );
    std::stringstream ss;
    ss << "cost: " << solution->currentCost() << "  " << solution->currentFrame();
    mTextureFont->drawString(ss.str(),  Vec2f(0, 20) );
}

CINDER_APP_NATIVE( TCO13MM3App, RendererGl )
