#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/TextureFont.h"
#include "CirclesSeparation.h"
#include <sstream>

#define READ_TEST_FROM_FILE 1
#define TEST_ID 6

#if READ_TEST_FROM_FILE
const int NUM_TESTS = 200;
int NS[NUM_TESTS];
double* XS[NUM_TESTS];
double* YS[NUM_TESTS];
double* RS[NUM_TESTS];
double* MS[NUM_TESTS];
#else
#include "TestCases.h"
#endif

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

    
#if READ_TEST_FROM_FILE
    stringstream ss;
    ss << "/Users/fuqinho/Dropbox/TCO13MM3/xcode/testcases" << NUM_TESTS << ".txt";
    freopen(ss.str().c_str(), "r", stdin);
    for (int i = 0; i < NUM_TESTS; i++) {
        int n; cin >> n;
        NS[i] = n;
        XS[i] = new double[n];
        YS[i] = new double[n];
        RS[i] = new double[n];
        MS[i] = new double[n];
        for (int j = 0; j < n; j++) cin >> XS[i][j];
        for (int j = 0; j < n; j++) cin >> YS[i][j];
        for (int j = 0; j < n; j++) cin >> RS[i][j];
        for (int j = 0; j < n; j++) cin >> MS[i][j];
    }
    solution->setup(vector<double>(XS[TEST_ID-1], XS[TEST_ID-1] + NS[TEST_ID-1]),
                    vector<double>(YS[TEST_ID-1], YS[TEST_ID-1] + NS[TEST_ID-1]),
                    vector<double>(RS[TEST_ID-1], RS[TEST_ID-1] + NS[TEST_ID-1]),
                    vector<double>(MS[TEST_ID-1], MS[TEST_ID-1] + NS[TEST_ID-1]));
#else
    solution->setup(vector<double>(XX, XX + sizeof(XX)/sizeof(double)),
                    vector<double>(YY, YY + sizeof(YY)/sizeof(double)),
                    vector<double>(RR, RR + sizeof(RR)/sizeof(double)),
                    vector<double>(MM, MM + sizeof(MM)/sizeof(double)));
#endif
}

void TCO13MM3App::mouseDown( MouseEvent event )
{
    solution->catch_circle(prj_left + event.getX() / prj_ratio,
                           prj_top + event.getY() / prj_ratio);
}

void TCO13MM3App::mouseUp( MouseEvent event )
{
    solution->release_circle();
}

void TCO13MM3App::mouseDrag( MouseEvent event )
{
    solution->move_circle(prj_left + event.getX() / prj_ratio,
                          prj_top + event.getY() / prj_ratio);
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
        left = min(left, ball_x[i] - ball_r[i]);
        right = max(right, ball_x[i] + ball_r[i]);
        top = min(top, ball_y[i] - ball_r[i]);
        bottom = max(bottom, ball_y[i] + ball_r[i]);
    }
    float rate_w = SCREEN_W / (right - left);
    float rate_h = SCREEN_H / (bottom - top);
    float rate = min(rate_w, rate_h);
    prj_left = left;
    prj_top = top;
    prj_ratio = rate;
    
    // precompute position on the screen
    for (int i = 0; i < N; i++) {
        cx[i] = (ball_x[i] - left) * rate;
        cy[i] = (ball_y[i] - top) * rate;
        cr[i] = ball_r[i] * rate;
        csx[i] = (ball_ox[i] - left) * rate;
        csy[i] = (ball_oy[i] - top) * rate;
    }
    
    float box_top = (0.0f - top) * rate;
    float box_bottom = (1.0f - top) * rate;
    float box_left = (0.0f - left) * rate;
    float box_right = (1.0f - left) * rate;
    
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
    
    // draw boxes
    gl::color( ColorA ( 2.0f, 1.0f, 2.0f, 0.7f ) );
    gl::lineWidth(2.0f);
    gl::drawLine(Vec2f(box_left, box_top), Vec2f(box_right, box_top));
    gl::drawLine(Vec2f(box_right, box_top), Vec2f(box_right, box_bottom));
    gl::drawLine(Vec2f(box_right, box_bottom), Vec2f(box_left, box_bottom));
    gl::drawLine(Vec2f(box_left, box_bottom), Vec2f(box_left, box_top));
    
    // draw current score
    gl::color( Color( 0.0f, 0.0f, 0.0f ) );
    std::stringstream ss;
    ss << "cost: " << solution->currentCost() << "  " << solution->currentFrame();
    mTextureFont->drawString(ss.str(),  Vec2f(0, 20) );
}

CINDER_APP_NATIVE( TCO13MM3App, RendererGl )
