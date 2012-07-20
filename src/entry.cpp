#include <stdio.h>
#include <string.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <SOIL.h>
#include <float.h>

#define PVEC2( v ) fprintf( stdout, "" #v " = %g, %g\n", v[0], v[1] )
#define PVEC3( v ) fprintf( stdout, "" #v " = %g, %g, %g\n", v[0], v[1], v[2] )
#define PVEC4( v ) fprintf( stdout, "" #v " = %g, %g, %g, %g\n", v[0], v[1], v[2], v[3] )


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
	glm::vec4 eye;
	glm::vec4 clip;
	glm::vec4 ndc;
	glm::vec2 screen;
	glm::vec3 normal;
	float u, v;
	float r, g, b, a;
	unsigned int color;
};

struct Viewport
{
	glm::vec2 offset;
	glm::vec2 size;
};

// client structures

struct Texture
{
	int width;
	int height;
	int channels;
	unsigned char * pixels;
};

struct RenderBuffer
{
	unsigned char * pixels;
	unsigned int width;
	unsigned int height;
	unsigned int channels;
	float * zbuffer;
};

/*
struct VertexShaderBlock
{
	glm::mat4 * modelview;
	glm::mat4 * projection;
	glm::vec3 * vertex_world;	// vertex in world space
	glm::vec4 * vertex_clip;	// vertex in clip space
};
*/

struct Triangle
{
	Vertex v[3];
	Texture * texture;
};

struct Settings
{
	int perspective_correct;
	int backface_culling;
	
	Settings()
	{
		perspective_correct = 1;
		backface_culling = 1;
	}
};

static Settings local_settings;

Settings &settings()
{
	return local_settings;
}

Triangle * triangles = 0;
unsigned int num_triangles = 0;

glm::vec3 light_position_eye;



void loadTexture( const char * filename, Texture & texture )
{
	texture.pixels = 0;
	texture.pixels = SOIL_load_image( filename, &texture.width, &texture.height, &texture.channels, SOIL_LOAD_AUTO );	
//	fprintf( stdout, "loaded image: \"%s\", %i x %i @ %i\n", filename, texture.width, texture.height, texture.channels );
}

void freeTexture( Texture & texture )
{
	if ( texture.pixels != 0 )
	{
		SOIL_free_image_data( texture.pixels );
	}
}

glm::vec4 texture2D( Texture * texture, const glm::vec2 & uv )
{
	int x, y;
	x = int(floor(uv.x * texture->width)) % texture->width;
	y = int(floor(uv.y * texture->height)) % texture->height;
	
//	fprintf( stdout, "u=%g, v=%g, x=%i, y=%i\n", uv.x, uv.y, x, y );
	
	unsigned char * pixel = &texture->pixels[ ((y * texture->width * texture->channels) + (x*texture->channels)) ];
	unsigned char alpha = 255;
	return glm::vec4( pixel[0], pixel[1], pixel[2], alpha );
}

inline unsigned int Color( unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a )
{
	return (_r | _g << 8 | _b << 16 | _a << 24 );
}

inline unsigned int scaleColor( unsigned int a, float t )
{
	float ac[] = { (a & 0xFF), ((a >> 8) & 0xFF), ((a >> 16) & 0xFF), ( (a >> 24) & 0xFF) };

	return (int(t*ac[0]) ) | (int(t*ac[1]) << 8) | (int(t*ac[2]) << 16) | (int(t*ac[3]) << 24);
}

inline void convertColor( unsigned int color, float out[4] )
{
	out[3] = ((color >> 24) & 0xFF) / 255.0;
	out[2] = ((color >> 16) & 0xFF) / 255.0;
	out[1] = ((color >> 8) & 0xFF) / 255.0;
	out[0] = (color & 0xFF) / 255.0;
}

void clearColor( RenderBuffer & rb, unsigned int color )
{
	unsigned char * p = rb.pixels;
	for( int h = 0; h < rb.height; ++h )
	{
		for( int w = 0; w < rb.width; ++w )
		{
			p[0] = (color) & 0xFF;
			p[1] = (color >> 8) & 0xFF;
			p[2] = (color >> 16) & 0xFF;
			p[3] = (color >> 24) & 0xFF;
			
			p += rb.channels;
		}
	}
}

void clearDepth( RenderBuffer & rb, float value )
{
	for( int h = 0; h < rb.height; ++h )
	{
		for( int w = 0; w < rb.width; ++w )
		{
			rb.zbuffer[ (h * rb.width) + w ] = value;
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

	//fminf( fminf( t->v[0].screen.x, t->v[1].screen.x ), t->v[2].screen.x );
	glm::vec3 light_dirA = light_position_eye - glm::vec3(t->v[0].eye);
	glm::vec3 light_dirB = light_position_eye - glm::vec3(t->v[1].eye);
	glm::vec3 light_dirC = light_position_eye - glm::vec3(t->v[2].eye);

	unsigned int idx;
	for( unsigned int py = mins.y; py < maxs.y; ++py )
	{
		for( unsigned int px = mins.x; px < maxs.x; ++px )
		{
			idx = (py * rb.width) + (px);
			unsigned char * pixel = &rb.pixels[ idx * rb.channels ];
			float alpha, beta, gamma;
			glm::vec2 point = glm::vec2(0.5f + px, 0.5f + py);
			barycentric_coords( point, t->v[0].screen, t->v[1].screen, t->v[2].screen, &alpha, &beta, &gamma );

			// determine if this pixel resides in the triangle or not
			if ( alpha > 0 && alpha < 1 && beta > 0 && beta < 1 && gamma > 0 && gamma < 1 )
			{
				// interpolate the value of w (clip space) at this pixel
//				w is not affine in screen space, but 1/w is (identical to z-view)
				float invw = (1 / t->v[0].clip.w) * alpha + (1/t->v[1].clip.w) * beta + (1/t->v[2].clip.w) * gamma;
				
				// package the w values up
				float w[] = { t->v[0].clip.w, t->v[1].clip.w, t->v[2].clip.w };
				
				if ( !settings().perspective_correct )
				{
					invw = 1;
					w[2] = w[1] = w[0] = 1;
				}
				
//				fprintf( stdout, "w: %g %g %g\n", w[0], w[1], w[2] );
				
				// calculate the depth value
				float depth = (t->v[0].ndc.z * alpha + t->v[1].ndc.z * beta + t->v[2].ndc.z * gamma);
				// perform a depth test
				if ( depth < rb.zbuffer[ idx ] )
				{
					rb.zbuffer[ idx ] = depth;

					glm::vec3 vertex_to_light = glm::normalize(light_position_eye - glm::vec3(t->v[0].eye));
					glm::vec3 normal = t->v[0].normal;
//					float ndl = fmax(0.0f, glm::dot( normal, vertex_to_light ) );
					float ndl = 1.0;
					

					
					glm::vec4 texel;
					if ( t->texture )
					{
						float u = (((t->v[0].u/w[0]) * alpha) + ((t->v[1].u/w[1]) * beta) + ((t->v[2].u/w[2]) * gamma)) / invw;
						float v = (((t->v[0].v/w[0]) * alpha) + ((t->v[1].v/w[1]) * beta) + ((t->v[2].v/w[2]) * gamma)) / invw;
						
						// perform texture lookup
						texel = texture2D( t->texture, glm::vec2(u,v) );
					}
					else
					{
						// interpolate the color
						float cA[4];
						float cB[4];
						float cC[4];						
						convertColor( t->v[0].color, cA );
						convertColor( t->v[1].color, cB );
						convertColor( t->v[2].color, cC );
						
						
						float r = ((cA[0]/w[0]) * alpha + (cB[0]/w[1]) * beta + (cC[0]/w[2]) * gamma) / invw;
						float g = ((cA[1]/w[0]) * alpha + (cB[1]/w[1]) * beta + (cC[1]/w[2]) * gamma) / invw;
						float b = ((cA[2]/w[0]) * alpha + (cB[2]/w[1]) * beta + (cC[2]/w[2]) * gamma) / invw;
						float a = ((cA[3]/w[0]) * alpha + (cB[3]/w[1]) * beta + (cC[3]/w[2]) * gamma) / invw;		
						texel = glm::vec4( r*255, g*255, b*255, a*255 );
					}

					pixel[0] = texel[0] * ndl;
					pixel[1] = texel[1] * ndl;
					pixel[2] = texel[2] * ndl;
					pixel[3] = texel[3];
				}
			}
		}
	}

}

void normalizedDeviceCoordsToScreen( const glm::vec2 & ndc, glm::vec2 & screenpos, const Viewport & viewport )
{
	//http://stackoverflow.com/questions/8491247/c-opengl-convert-world-coords-to-screen2d-coords
	screenpos[0] = ((ndc[0]+1.0) / 2.0 * viewport.size[0] + viewport.offset[0]);
	screenpos[1] = ((ndc[1]+1.0) / 2.0 * viewport.size[1] + viewport.offset[1]);
}


void renderScene( const Camera & camera, const Viewport & viewport, RenderBuffer & rb )
{
	clearColor( rb, Color(0,0,0,255) );
	clearDepth( rb, 9 );
	
	glm::vec3 lightposition( 5, 10, 10 );
//	glm::mat3 normal_matrix = glm::inverseTranspose( glm::mat3( camera.modelview ) );
	light_position_eye = glm::vec3(camera.modelview * glm::vec4(lightposition, 1.0) );	
	
	Triangle * t = triangles;
	// for every object visible in the scene
	for( int i = 0; i < num_triangles; ++i )
	{
		// convert verts into eye space (space of the camera, in this regard)
//		PVEC3( t->v[0].position );
		
		t->v[0].eye = camera.modelview * glm::vec4( t->v[0].position, 1.0 );
		t->v[1].eye = camera.modelview * glm::vec4( t->v[1].position, 1.0 );	
		t->v[2].eye = camera.modelview * glm::vec4( t->v[2].position, 1.0 );
//		PVEC4( t->v[0].eye );
		
		// convert triangle to clip space (projection)
		// this places the coordinates in [-1, 1] clip space
		t->v[0].clip = camera.projection * t->v[0].eye;
		t->v[1].clip = camera.projection * t->v[1].eye;
		t->v[2].clip = camera.projection * t->v[2].eye;	
//		PVEC4( t->v[0].clip );
		
		// divide by w, convert to normalized device coordinates [0, 1]
		t->v[0].ndc = t->v[0].clip * (1/t->v[0].clip.w);
		t->v[1].ndc = t->v[1].clip * (1/t->v[1].clip.w);
		t->v[2].ndc = t->v[2].clip * (1/t->v[2].clip.w);
//		PVEC4( t->v[0].ndc );
		
		// convert triangle to screen coordinates; [viewport.offset - viewport.size]
		normalizedDeviceCoordsToScreen( glm::vec2(t->v[0].ndc), t->v[0].screen, viewport );
		normalizedDeviceCoordsToScreen( glm::vec2(t->v[1].ndc), t->v[1].screen, viewport );
		normalizedDeviceCoordsToScreen( glm::vec2(t->v[2].ndc), t->v[2].screen, viewport );
//		PVEC2( t->v[0].screen );		
		
		// perform culling
		glm::vec3 edge1 = glm::vec3(t->v[1].ndc - t->v[0].ndc);
		glm::vec3 edge2 = glm::vec3(t->v[2].ndc - t->v[0].ndc);
		glm::vec3 normal = glm::normalize( glm::cross( edge2, edge1 ) );
		
		if ( !settings().backface_culling )
			normal.z = 1;
		
		// TODO: this test should eventually dot normal and the inverted view vector from the camera
		// currently works fine if camera is viewing down -Z axis.
		if ( normal.z > 0 )
		{		
			renderTriangle( rb, t );
		}
		t++;
	}
}

int main( int argc, char ** argv )
{	
	// setup a viewport
	Viewport vp;
	vp.offset = glm::vec2( 0, 0 );
	vp.size = glm::vec2( VIEWPORT_WIDTH, VIEWPORT_HEIGHT );

	
	fprintf( stdout, "======== Settings ========\n" );
	fprintf( stdout, "\tperpsective_correct: %i\n", settings().perspective_correct );
	fprintf( stdout, "\tbackface_culling: %i\n", settings().backface_culling );
	
	// setup camera
	Camera cam1;
	cam1.modelview = glm::translate( glm::mat4(1.0f), glm::vec3( 0, 0, -1.5 ));
	cam1.projection = glm::perspective( 60.0f, ((float)VIEWPORT_WIDTH/(float)VIEWPORT_HEIGHT), .1f, 2048.0f );
	cam1.modelview_projection = cam1.projection * cam1.modelview;

	
	cam1.modelview = glm::rotate( cam1.modelview, 45.0f, glm::vec3( 1, 0, 0 ) );

	RenderBuffer render_buffer;
	render_buffer.width = (int)VIEWPORT_WIDTH;
	render_buffer.height = (int)VIEWPORT_HEIGHT;
	render_buffer.channels = 4;
	render_buffer.pixels = new unsigned char[ render_buffer.width * render_buffer.height * render_buffer.channels ];
	render_buffer.zbuffer = new float[ render_buffer.width * render_buffer.height ];

	Texture tex;
	
	loadTexture( "assets/checker2.png", tex );	
	
	num_triangles = 2;
	triangles = new Triangle[ num_triangles ];
	Triangle t1;
	t1.v[0].position = glm::vec3( -0.5f, -0.5f, -0.0f );
	t1.v[0].color = Color(255, 0, 0, 255);
	t1.v[0].normal = glm::normalize(glm::vec3( 0, 0, 1 ));
	t1.v[0].u = 0.0; t1.v[0].v = 0.0;	
	t1.v[1].position = glm::vec3( -0.5f, 0.5, -0.0f );
	t1.v[1].color = Color(0, 0, 255, 255);
	t1.v[1].normal = glm::normalize(glm::vec3( 0, 0, 1 ));	
	t1.v[1].u = 0.0; t1.v[1].v = 1.0;
	t1.v[2].position = glm::vec3( 0.5f, 0.5, -0.0 );
	t1.v[2].color = Color(0, 0, 255, 255);
	t1.v[2].normal = glm::normalize(glm::vec3( 0, 0, 1 ));
	t1.v[2].u = 1.0; t1.v[2].v = 1.0;	
	t1.texture = 0;
	t1.texture = &tex;	
	triangles[0] = t1;
	
	Triangle t2;
	t2.v[0].position = t1.v[2].position;
	t2.v[0].color = t1.v[2].color;
	t2.v[0].normal = t1.v[2].normal;
	t2.v[0].u = t1.v[2].u; t2.v[0].v = t1.v[2].v;
	t2.v[1].position = glm::vec3( 0.5, -0.5, 0.0f );
	t2.v[1].color = Color(255, 0, 0, 255);
	t2.v[1].normal = glm::normalize(glm::vec3( 0, 0, 1 ));	
	t2.v[1].u = 1.0; t2.v[1].v = 0.0;
	t2.v[2].position = t1.v[0].position;
	t2.v[2].color = t1.v[0].color;
	t2.v[2].normal = t1.v[0].normal;
	t2.v[2].u = t1.v[0].u; t2.v[2].v = t1.v[0].v;	
	t2.texture = 0;
	t2.texture = &tex;
	triangles[1] = t2;
	
	
	
	
	/*
	t1.v[0].position = glm::vec3( -0.4f, 0.2f, -0.2f );
	t1.v[0].color = Color(0, 0, 255, 255);
	t1.v[0].normal = glm::normalize(glm::vec3( 1, 0, 1 ));
	t1.v[0].u = 0.0; t1.v[0].v = 0.0;	
	t1.v[1].position = glm::vec3( -0.2, -0.6, -1.6f );
	t1.v[1].color = Color(0, 0, 255, 255);
	t1.v[1].normal = glm::normalize(glm::vec3( 1, 0.35, 0.5 ));	
	t1.v[1].u = 0.0; t1.v[1].v = 1.0;
	t1.v[2].position = glm::vec3( 0.5, 0.9, -0.5f );
	t1.v[2].color = Color(0, 255, 0, 255);
	t1.v[2].normal = glm::normalize(glm::vec3( 1, -0.55, 1 ));
	t1.v[2].u = 1.0; t1.v[2].v = 1.0;	
	t1.texture = &tex;
	triangles[0] = t1;
	
	t1.v[0].position = glm::vec3( -0.7, 0.6, -1.2f );
	t1.v[0].color = Color(255, 0, 0, 255);
	t1.v[0].normal = glm::vec3( 0, 0, 1 );
	t1.v[0].u = 0.0; t1.v[0].v = 1.0;
	t1.v[1].position = glm::vec3( -0.2, -0.3, -0.8f );
	t1.v[1].color = Color(0, 255, 0, 255);
	t1.v[1].normal = glm::vec3( 0, 0, 1 );
	t1.v[1].u = 1.0; t1.v[1].v = 1.0;	
	t1.v[2].position = glm::vec3( 0.8, 0.2, -0.2f );
	t1.v[2].color = Color(0, 0, 255, 255);
	t1.v[2].normal = glm::vec3( 0, 0, 1 );
	t1.v[2].u = 1.0; t1.v[2].v = 0.0;	
	t1.texture = &tex;
	triangles[1] = t1;
	*/

	fprintf( stdout, "======== Render Scene ========\n" );
	renderScene( cam1, vp, render_buffer );
	fprintf( stdout, "\ttriangles rendered: %i\n", num_triangles );

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
			
			glBegin( GL_TRIANGLES );
			Triangle * t = triangles;
			for( int i = 0; i < num_triangles; ++i )
			{
				glColor4ubv( (unsigned char*)&t->v[0].color );
				glVertex3fv( (const float*)&t->v[0].position[0] );
				glColor4ubv( (unsigned char*)&t->v[1].color );
				glVertex3fv( (const float*)&t->v[1].position[0] );
				glColor4ubv( (unsigned char*)&t->v[2].color );			
				glVertex3fv( (const float*)&t->v[2].position[0] );
				++t;
			}
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

	int save_result = SOIL_save_image( "output/output.bmp", SOIL_SAVE_TYPE_BMP, render_buffer.width, render_buffer.height, render_buffer.channels, render_buffer.pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}

	delete [] render_buffer.pixels;
	delete [] render_buffer.zbuffer;
	
	freeTexture( tex );

	return 0;
}
