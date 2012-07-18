#include <stdio.h>
#include <string.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <SOIL.h>


#define PVEC2( v ) fprintf( stdout, "" #v "= %g, %g\n", v[0], v[1] )
#define PVEC3( v ) fprintf( stdout, "" #v "= %g, %g, %g\n", v[0], v[1], v[2] )
#define PVEC4( v ) fprintf( stdout, "" #v "= %g, %g, %g, %g\n", v[0], v[1], v[2], v[3] )


#define DEBUG_WINDOW 1
#if DEBUG_WINDOW
	#include <xwl/xwl.h>
	#if __APPLE__
		#include <OpenGL/gl.h>
	#else
		#include <GL/gl.h>
	#endif
	int run = 1;
	void callback( xwl_event_t * e )
	{
		if ( e->type == XWLE_KEYRELEASED )
		{
			if ( e->key == XWLK_ESCAPE )
			{
				run = 0;
			}
			//printf( "\t-> key: %i (%s)\n", e->key, xwl_key_to_string(e->key) );
		}		
	}

#endif

#define VIEWPORT_WIDTH 640.0f
#define VIEWPORT_HEIGHT 480.0f

// object coordinates * (model view matrix) = eye coordinates
// eye coordinates * (projection matrix ) = clip coordinates
// clip coordinates / (perspective divide) = normalized device coordinates
// normalized device coordinates * viewport transform = screen coordinates

// equation from wikipedia: http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
void barycentric_coords( const glm::vec2 & p, const glm::vec2 & a, const glm::vec2 & b, const glm::vec2 & c, float * alpha, float * beta, float * gamma )
{
	float detT = 1.0 / ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	*alpha = ((b.y - c.y)*(p.x - c.x) + (c.x - b.x) * (p.y - c.y)) * detT;
	*beta = ((c.y - a.y)*(p.x - c.x) + (a.x - c.x) * (p.y - c.y)) * detT;

	*gamma = 1.0f - *alpha - *beta;
}

struct Camera
{
	glm::mat4 modelview;
	glm::mat4 projection;
	glm::mat4 modelview_projection;
};

struct Vertex
{
	glm::vec3 position;
	glm::vec4 ndc;
	glm::vec2 screen;
	float r, g, b, a;
	unsigned int color;
};

struct Viewport
{
	glm::vec2 offset;
	glm::vec2 size;
};

// client structures

struct RenderBuffer
{
	unsigned char * pixels;
	unsigned int width;
	unsigned int height;
	unsigned int channels;
	float * zbuffer;
};

struct Triangle
{
	Vertex v[3];
};

Triangle * triangles = 0;
unsigned int num_triangles = 0;


inline unsigned int Color( unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a )
{
	return (_r << 24 | _g << 16 | _b << 8 | _a );
}

inline unsigned int scaleColor( unsigned int a, float t )
{
	float ac[] = { ((a >> 24) & 0xFF)/255.0f, ((a >> 16) & 0xFF)/255.0f, ((a >> 8) & 0xFF)/255.0f, (a & 0xFF)/255.0f };

	return (int(t*ac[0]*255) << 24 ) | (int(t*ac[1]*255) << 16) | (int(t*ac[2]*255) << 8) | (int(t*ac[3]*255));
}




void clearColor( RenderBuffer & rb, unsigned int color )
{
	unsigned char * p = rb.pixels;
	for( int h = 0; h < VIEWPORT_HEIGHT; ++h )
	{
		for( int w = 0; w < VIEWPORT_WIDTH; ++w )
		{
			p[0] = (color >> 24) & 0xFF;
			p[1] = (color >> 16) & 0xFF;
			p[2] = (color >> 8) & 0xFF;
			p[3] = color & 0xFF;
			
			p += rb.channels;
		}
	}
}


void renderTriangle( RenderBuffer & rb, Triangle * t )
{
	// find the bounding rect of the triangle in screen space	
	glm::vec2 mins(rb.width, rb.height), maxs;
	for( int i = 0; i < 3; ++i )
	{
		if ( t->v[i].screen.x < mins.x )
			mins.x = t->v[i].screen.x;
		if ( t->v[i].screen.y < mins.y )
			mins.y = t->v[i].screen.y;
		if ( t->v[i].screen.x > maxs.x )
			maxs.x = t->v[i].screen.x;
		if ( t->v[i].screen.y > maxs.y )
			maxs.y = t->v[i].screen.y;		
	}

	mins.x = floor(mins.x);
	mins.y = floor(mins.y);
	maxs.x = floor(maxs.x);
	maxs.y = floor(maxs.y);	

	//unsigned int max_width = (maxs.x - mins.x);
	//unsigned int max_height = (maxs.y - mins.y);
	//fprintf( stdout, "rect: %i x %i\n", max_width, max_height );
	
	// I invert the direction that the height is traversed to match the same image output as OpenGL
	// Alternatively, I could simply flip the image when I'm done...
//	for( int h = VIEWPORT_HEIGHT; h > 0; --h )	
	for( unsigned int py = mins.y; py < maxs.y; ++py )
	{
		for( unsigned int px = mins.x; px < maxs.x; ++px )
		{
			unsigned char * pixel = &rb.pixels[ (py * rb.width * rb.channels) + (px*rb.channels) ];
			float alpha, beta, gamma;
			glm::vec2 point = glm::vec2(0.5f + px, 0.5f + py);
			barycentric_coords( point, t->v[0].screen, t->v[1].screen, t->v[2].screen, &alpha, &beta, &gamma );

			// determine if this pixel resides in the triangle or not
			if ( alpha > 0 && alpha < 1 && beta > 0 && beta < 1 && gamma > 0 && gamma < 1 )
			{
				// interpolate the color
				unsigned int p = scaleColor(t->v[0].color, alpha) + scaleColor(t->v[1].color, beta) + scaleColor(t->v[2].color, gamma); 
				unsigned char r = (p >> 24) & 0xFF;
				unsigned char g = (p >> 16) & 0xFF;
				unsigned char b = (p >> 8) & 0xFF;
				unsigned char a = p & 0xFF;
				pixel[0] = r;
				pixel[1] = g;
				pixel[2] = b;
				pixel[3] = a;
				//fprintf( stdout, "p = %i\n", (*p >> 24) & 0xFF );
				//pixel[0] = (vertices[0].r * alpha + vertices[1].r * beta + vertices[2].r * gamma) * 255.0f;
				//pixel[1] = (vertices[0].g * alpha + vertices[1].g * beta + vertices[2].g * gamma) * 255.0f;
				//pixel[2] = (vertices[0].b * alpha + vertices[1].b * beta + vertices[2].b * gamma) * 255.0f;
				//pixel[3] = (vertices[0].a * alpha + vertices[1].a * beta + vertices[2].a * gamma) * 255.0f;
				
				//float z = (a.z * alpha + b.z * beta + c.z * gamma);
				
				//p[0] = p[1] = p[2] = z*255.0f;
				//p[3] = 255;
			}
		}
	}

}

void worldSpaceToClipSpace( const glm::mat4 & modelview_projection, glm::vec3 & world, glm::vec4 & out )
{
	out = modelview_projection * glm::vec4(world, 1.0);
}

void normalizedDeviceCoordsToScreen( const glm::vec2 & ndc, glm::vec2 & screenpos, const Viewport & viewport )
{
	//http://stackoverflow.com/questions/8491247/c-opengl-convert-world-coords-to-screen2d-coords
	screenpos[0] = ((ndc[0]+1.0) / 2.0 * viewport.size[0] + viewport.offset[0]);
	screenpos[1] = ((ndc[1]+1.0) / 2.0 * viewport.size[1] + viewport.offset[1]);
}


void renderScene( const Camera & camera, const Viewport & viewport, RenderBuffer & rb )
{
	Triangle * t = triangles;
	// for every object visible in the scene
	for( int i = 0; i < num_triangles; ++i )
	{
		// convert triangle to clip space
		worldSpaceToClipSpace( camera.modelview_projection, t->v[0].position, t->v[0].ndc );
		worldSpaceToClipSpace( camera.modelview_projection, t->v[1].position, t->v[1].ndc );
		worldSpaceToClipSpace( camera.modelview_projection, t->v[2].position, t->v[2].ndc );

		// divide by w, convert to normalized device coordinates
		t->v[0].ndc *= (1/t->v[0].ndc.w);
		t->v[1].ndc *= (1/t->v[1].ndc.w);
		t->v[2].ndc *= (1/t->v[2].ndc.w);

		// convert triangle to screen coordinates
		normalizedDeviceCoordsToScreen( glm::vec2(t->v[0].ndc), t->v[0].screen, viewport );
		normalizedDeviceCoordsToScreen( glm::vec2(t->v[1].ndc), t->v[1].screen, viewport );
		normalizedDeviceCoordsToScreen( glm::vec2(t->v[2].ndc), t->v[2].screen, viewport );

		t++;
	}



	// fill image with a color
	unsigned char * p = rb.pixels;

	clearColor( rb, Color(0,0,0,255) );

	renderTriangle( rb, &triangles[0] );
}

int main( int argc, char ** argv )
{	
	// setup a viewport
	Viewport vp;
	vp.offset = glm::vec2( 0, 0 );
	vp.size = glm::vec2( VIEWPORT_WIDTH, VIEWPORT_HEIGHT );

	// setup camera
	Camera cam1;
	cam1.modelview = glm::translate( glm::mat4(1.0f), glm::vec3( 0, 0, -2 ));
//	cam1.projection = glm::ortho( -VIEWPORT_WIDTH/2.0f, VIEWPORT_WIDTH/2.0f, VIEWPORT_HEIGHT/2.0f, -VIEWPORT_HEIGHT/2.0f, -0.1f, 128.0f );
	cam1.projection = glm::perspective( 60.0f, ((float)VIEWPORT_WIDTH/(float)VIEWPORT_HEIGHT), .1f, 2048.0f );
	cam1.modelview_projection = cam1.projection * cam1.modelview;


	RenderBuffer render_buffer;
	render_buffer.width = (int)VIEWPORT_WIDTH;
	render_buffer.height = (int)VIEWPORT_HEIGHT;
	render_buffer.channels = 4;
	render_buffer.pixels = new unsigned char[ render_buffer.width * render_buffer.height * render_buffer.channels ];
	render_buffer.zbuffer = new float[ render_buffer.width * render_buffer.height ];
//	vertices[0].position = glm::vec3( -50, 0, 0 );
//	vertices[0].r = 1.0f; vertices[0].g = vertices[0].b = 0.0f; vertices[0].a = 1.0f;
//	vertices[1].position = glm::vec3( 50, 0, 0 );
//	vertices[1].r = 0.0f; vertices[1].g = 1.0f; vertices[1].b = 0.0f; vertices[1].a = 1.0f;
//	vertices[2].position = glm::vec3( 0, -50, 0 );
//	vertices[2].r = 0.0f; vertices[2].g = 0.0f; vertices[2].b = 1.0f; vertices[2].a = 1.0f;

	num_triangles = 1;
	triangles = new Triangle[ num_triangles ];
	Triangle t1;
	t1.v[0].position = glm::vec3( -0.4f, 0.2f, -0.2f );
	t1.v[0].color = Color(255, 0, 0, 255);
	t1.v[1].position = glm::vec3( -0.2, -0.6, -0.6f );
	t1.v[1].color = Color(0, 255, 0, 255);
	t1.v[2].position = glm::vec3( 0.2, 0.9, -0.3f );
	t1.v[2].color = Color(0, 0, 255, 255);

	triangles[0] = t1;
	/*
	// three points of a triangle
	Vertex vertices[3];	
	vertices[0].position = glm::vec3( -0.4, 0.2, -0.2f );
	vertices[0].r = 0.0f; 
	vertices[0].g = 0.0f;
	vertices[0].b = 1.0f;
	vertices[0].a = 1.0f;
	vertices[1].position = glm::vec3( -0.2, -0.6, -0.6f );
	vertices[1].r = 0.0f;
	vertices[1].g = 1.0f;
	vertices[1].b = 0.0f;
	vertices[1].a = 1.0f;
	vertices[2].position = glm::vec3( 0.2, 0.9, -0.3f );
	vertices[2].r = 0.0f;
	vertices[2].g = 0.0f;
	vertices[2].b = 1.0f;
	vertices[2].a = 1.0f;

	vertices[0].position = glm::vec3( -0.7, 0.6, -1.2f );
	vertices[0].r = 1.0f; 
	vertices[0].g = 0.0f;
	vertices[0].b = 0.0f;
	vertices[0].a = 1.0f;
	vertices[1].position = glm::vec3( -0.2, -0.3, -0.8f );
	vertices[1].r = 0.0f;
	vertices[1].g = 1.0f;
	vertices[1].b = 0.0f;
	vertices[1].a = 1.0f;
	vertices[2].position = glm::vec3( 0.8, 0.2, -0.2f );
	vertices[2].r = 0.0f;
	vertices[2].g = 0.0f;
	vertices[2].b = 1.0f;
	vertices[2].a = 1.0f;	
*/


	renderScene( cam1, vp, render_buffer );
	/*
	fprintf( stdout, "World Coordinates:\n" );
	PVEC3( vertices[0].position );
	PVEC3( vertices[1].position );
	PVEC3( vertices[2].position );


	// convert triangle to clip space
	glm::vec4 a, b, c;
	worldSpaceToClipSpace( cam1.modelview_projection, vertices[0].position, a );
	worldSpaceToClipSpace( cam1.modelview_projection, vertices[1].position, b );
	worldSpaceToClipSpace( cam1.modelview_projection, vertices[2].position, c );

	fprintf( stdout, "Converted to Clip Coordinates:\n" );
	PVEC4(a);
	PVEC4(b);
	PVEC4(c);	
	
	fprintf( stdout, "Converted to Normalized Device Coordinates:\n" );
	// convert to normalized device coordinates
	a *= (1/a.w);
	b *= (1/b.w);
	c *= (1/c.w);

	PVEC4(a);
	PVEC4(b);
	PVEC4(c);

	fprintf( stdout, "Converted to Screen Coordinates:\n" );
	// convert triangle to screen coordinates
	glm::vec2 wa, wb, wc;
	normalizedDeviceCoordsToScreen( glm::vec2(a), wa, vp );
	normalizedDeviceCoordsToScreen( glm::vec2(b), wb, vp );
	normalizedDeviceCoordsToScreen( glm::vec2(c), wc, vp );

	PVEC2(wa);
	PVEC2(wb);
	PVEC2(wc);
	*/



#if DEBUG_WINDOW
	if ( argc > 1 )
	{
		xwl_startup();

		xwl_windowparams_t wp;
		wp.width = VIEWPORT_WIDTH;
		wp.height = VIEWPORT_HEIGHT;
		wp.flags |= XWL_OPENGL;
		xwl_window_t * window = xwl_create_window( &wp, "title", 0 );
		xwl_set_callback( callback );

		printf( "-> GL_VENDOR: %s\n", glGetString( GL_VENDOR ) );
		printf( "-> GL_RENDERER: %s\n", glGetString( GL_RENDERER ) );
		printf( "-> GL_VERSION: %s\n", glGetString( GL_VERSION ) );

		
		while( run )
		{
			xwl_event_t e;
			memset( &e, 0, sizeof(xwl_event_t) );
			xwl_pollevent( &e );

			glClearColor( 0.25, 0.25, 0.25, 1.0 );
			glClearDepth( 1.0f );
			glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
			glDisable( GL_CULL_FACE );
			glViewport( vp.offset[0], vp.offset[1], vp.size[0], vp.size[1] );

			glMatrixMode( GL_PROJECTION );
			glLoadMatrixf( (const float*)&cam1.projection[0] );

			glMatrixMode( GL_MODELVIEW );
			glLoadMatrixf( (const float*)&cam1.modelview[0] );
/*
			glColor3ub( 255, 255, 255 );
			glBegin( GL_TRIANGLES );
				glColor4fv( &vertices[0].r );
				glVertex3fv( (const float*)&vertices[0].position[0] );
				glColor4fv( &vertices[1].r );
				glVertex3fv( (const float*)&vertices[1].position[0] );
				glColor4fv( &vertices[2].r );			
				glVertex3fv( (const float*)&vertices[2].position[0] );
			glEnd();
*/
			glPointSize( 4.0 );
			glColor3ub( 255, 0, 0 );

			glBegin( GL_POINTS );
				glVertex3i( 0, 0, 0 );
			glEnd();

			xwl_finish();
		}


		xwl_shutdown();
	}
#endif

	int save_result = SOIL_save_image( "output/output.bmp", SOIL_SAVE_TYPE_BMP, render_buffer.width, render_buffer.height, render_buffer.channels, render_buffer.pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}

	delete [] render_buffer.pixels;
	delete [] render_buffer.zbuffer;

	return 0;
}
