#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class TCO13MM3App : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void TCO13MM3App::setup()
{
}

void TCO13MM3App::mouseDown( MouseEvent event )
{
}

void TCO13MM3App::update()
{
}

void TCO13MM3App::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 
}

CINDER_APP_NATIVE( TCO13MM3App, RendererGl )
