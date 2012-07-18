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

struct Vertex
{
	glm::vec3 position;
	float r, g, b, a;
};

struct Triangle
{
	Vertex v[3];	
};

int main( int argc, char ** argv )
{
	const int channels = 4;
	unsigned char * pixels = new unsigned char[ (int)VIEWPORT_WIDTH * (int)VIEWPORT_HEIGHT * channels ];
	

	// three points of a triangle
	
	Vertex vertices[3];
	
//	vertices[0].position = glm::vec3( -50, 0, 0 );
//	vertices[0].r = 1.0f; vertices[0].g = vertices[0].b = 0.0f; vertices[0].a = 1.0f;
//	vertices[1].position = glm::vec3( 50, 0, 0 );
//	vertices[1].r = 0.0f; vertices[1].g = 1.0f; vertices[1].b = 0.0f; vertices[1].a = 1.0f;
//	vertices[2].position = glm::vec3( 0, -50, 0 );
//	vertices[2].r = 0.0f; vertices[2].g = 0.0f; vertices[2].b = 1.0f; vertices[2].a = 1.0f;
	
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

	fprintf( stdout, "World Coordinates:\n" );
	PVEC3( vertices[0].position );
	PVEC3( vertices[1].position );
	PVEC3( vertices[2].position );


	glm::mat4 modelview = glm::mat4(1.0); // world to eye
	modelview = glm::translate( modelview, glm::vec3( 0, 0, -2 ));

	glm::mat4 projection; // eye to clip
//	projection = glm::ortho( -VIEWPORT_WIDTH/2.0f, VIEWPORT_WIDTH/2.0f, VIEWPORT_HEIGHT/2.0f, -VIEWPORT_HEIGHT/2.0f, -0.1f, 128.0f );
//	glm::mat4 rotate = glm::rotate( glm::mat4(1.0), 5.0f, glm::vec3(1,0,0) );
	
	projection = glm::perspective( 60.0f, ((float)VIEWPORT_WIDTH/(float)VIEWPORT_HEIGHT), .1f, 2048.0f );

	// convert triangle to clip space
	glm::mat4 world_to_clip = (projection * modelview);
	glm::vec4 a = world_to_clip * glm::vec4( vertices[0].position, 1.0 );
	glm::vec4 b = world_to_clip * glm::vec4( vertices[1].position, 1.0 );
	glm::vec4 c = world_to_clip * glm::vec4( vertices[2].position, 1.0 );

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

	//http://stackoverflow.com/questions/8491247/c-opengl-convert-world-coords-to-screen2d-coords
	// convert triangle to screen coordinates
	glm::vec2 viewOffset( 0, 0 );
	glm::vec2 viewSize( VIEWPORT_WIDTH, VIEWPORT_HEIGHT );
	
	glm::vec2 wa = glm::vec2(
	wa[0] = (a[0]+1.0) / 2.0 * viewSize[0] + viewOffset[0],
	wa[1] = (a[1]+1.0) / 2.0 * viewSize[1] + viewOffset[1] );

	glm::vec2 wb = glm::vec2(
	wb[0] = (b[0]+1.0) / 2.0 * viewSize[0] + viewOffset[0],
	wb[1] = (b[1]+1.0) / 2.0 * viewSize[1] + viewOffset[1] );

	glm::vec2 wc = glm::vec2(
	wc[0] = (c[0]+1.0) / 2.0 * viewSize[0] + viewOffset[0],
	wc[1] = (c[1]+1.0) / 2.0 * viewSize[1] + viewOffset[1] );

	PVEC2(wa);
	PVEC2(wb);
	PVEC2(wc);

	// fill image with a color
	unsigned char * p = pixels;
	// I invert the direction that the height is traversed to match the same image output as OpenGL
	// Alternatively, I could simply flip the image when I'm done...
//	for( int h = VIEWPORT_HEIGHT; h > 0; --h )
	for( int h = 0; h < VIEWPORT_HEIGHT; ++h )
	{
		for( int w = 0; w < VIEWPORT_WIDTH; ++w )
		{
			float alpha, beta, gamma;
			
			// use the pixel center
			glm::vec2 ndc = glm::vec2(w + 0.5f, h + 0.5f);
			
			if ( 0 ) // normalized device coordinate space
			{
				// convert screen coords back to normalized device coordinates
				ndc.x = ((ndc.x - viewOffset[0]) / (viewSize[0])) * 2.0 - 1.0;
				ndc.y = ((ndc.y - viewOffset[1]) / (viewSize[1])) * 2.0 - 1.0;

				// calculate barycentric coordinates
				barycentric_coords( ndc, glm::vec2(a), glm::vec2(b), glm::vec2(c), &alpha, &beta, &gamma );
			}
			else if ( 1 ) // screen space
			{			
				// calculate barycentric coordinates
				barycentric_coords( ndc, wa, wb, wc, &alpha, &beta, &gamma );
			}

			
			// determine if this pixel resides in the triangle or not
			if ( alpha > 0 && alpha < 1 && beta > 0 && beta < 1 && gamma > 0 && gamma < 1 )
			{
				// interpolate the color
				p[0] = (vertices[0].r * alpha + vertices[1].r * beta + vertices[2].r * gamma) * 255.0f;
				p[1] = (vertices[0].g * alpha + vertices[1].g * beta + vertices[2].g * gamma) * 255.0f;
				p[2] = (vertices[0].b * alpha + vertices[1].b * beta + vertices[2].b * gamma) * 255.0f;
				p[3] = (vertices[0].a * alpha + vertices[1].a * beta + vertices[2].a * gamma) * 255.0f;
				
				//float z = (a.z * alpha + b.z * beta + c.z * gamma);
				
				//p[0] = p[1] = p[2] = z*255.0f;
				//p[3] = 255;
			}
			else
			{
				p[0] = 32;
				p[1] = 32;
				p[2] = 32;
				p[3] = 255;
			}

			p += channels;
		}
	}

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
			glViewport( viewOffset[0], viewOffset[1], viewSize[0], viewSize[1] );

			glMatrixMode( GL_PROJECTION );
			glLoadMatrixf( (const float*)&projection[0] );

			glMatrixMode( GL_MODELVIEW );
			glLoadMatrixf( (const float*)&modelview[0] );

			glColor3ub( 255, 255, 255 );
			glBegin( GL_TRIANGLES );
				glColor4fv( &vertices[0].r );
				glVertex3fv( (const float*)&vertices[0].position[0] );
				glColor4fv( &vertices[1].r );
				glVertex3fv( (const float*)&vertices[1].position[0] );
				glColor4fv( &vertices[2].r );			
				glVertex3fv( (const float*)&vertices[2].position[0] );
			glEnd();

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

	int save_result = SOIL_save_image( "output/output.bmp", SOIL_SAVE_TYPE_BMP, VIEWPORT_WIDTH, VIEWPORT_HEIGHT, channels, pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}

	delete [] pixels;

	return 0;
}
