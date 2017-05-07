#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class Particles_2DApp : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );	
	void update();
	void draw();
};

void Particles_2DApp::setup()
{
}

void Particles_2DApp::mouseDown( MouseEvent event )
{
}

void Particles_2DApp::update()
{
}

void Particles_2DApp::draw()
{
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) ); 
}

CINDER_APP_NATIVE( Particles_2DApp, RendererGl )
