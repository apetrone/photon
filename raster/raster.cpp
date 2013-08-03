#include <stdio.h>
#include <string.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <SOIL.h>
#include <float.h>
#include "common.h"

#define PVEC2( v ) fprintf( stdout, "" #v " = %g, %g\n", v[0], v[1] )
#define PVEC3( v ) fprintf( stdout, "" #v " = %g, %g, %g\n", v[0], v[1], v[2] )
#define PVEC4( v ) fprintf( stdout, "" #v " = %g, %g, %g, %g\n", v[0], v[1], v[2], v[3] )


#define VIEWPORT_WIDTH 640.0f
#define VIEWPORT_HEIGHT 480.0f

// object coordinates * (model view matrix) = eye coordinates
// eye coordinates * (projection matrix ) = clip coordinates
// clip coordinates / (perspective divide) = normalized device coordinates
// normalized device coordinates * viewport transform = screen coordinates


/*
struct VertexShaderBlock
{
	glm::mat4 * modelview;
	glm::mat4 * projection;
	glm::vec3 * vertex_world;	// vertex in world space
	glm::vec4 * vertex_clip;	// vertex in clip space
};
*/



struct Line
{
	Vertex start;
	Vertex end;
};

struct Settings
{
	int perspective_correct;
	int backface_culling;
	int texture_mapping;
	int lighting;
	
	Settings()
	{
		perspective_correct = 1;
		backface_culling = 1;
		texture_mapping = 1;
		lighting = 1;
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
			if ( alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 )
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
					glm::vec3 vertex_to_light = glm::normalize( light_position_eye - glm::vec3( (((t->v[0].eye/w[0]) * alpha) + ((t->v[1].eye/w[1]) * beta) + ((t->v[2].eye/w[2]) * gamma))/invw ) );
					glm::vec3 normal = glm::normalize((((t->v[0].normal_eye/w[0]) * alpha + ((t->v[1].normal_eye/w[1]) * beta) + ((t->v[2].normal_eye/w[2]) * gamma))) / invw );
					float ndl = fmax(0.0f, glm::dot( normal, vertex_to_light ) );
					if ( !settings().lighting )
						ndl = 1.0;
					
					glm::vec4 texel;
					if ( t->texture && settings().texture_mapping )
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

unsigned char * getPixelAt( RenderBuffer & rb, int x, int y )
{
	unsigned int idx = (y * rb.width * rb.channels) + (x*rb.channels);
	return &rb.pixels[ idx ];
}

// http://roguebasin.roguelikedevelopment.org/index.php/Bresenham's_Line_Algorithm
void renderLine( RenderBuffer & rb, Line * l )
{
	unsigned char * pixel;
	float cA[4];
	convertColor( l->start.color, cA );
	
	int x0 = l->start.position.x;
	int x1 = l->end.position.x;
	int y0 = l->start.position.y;
	int y1 = l->end.position.y;
	
	
	int dx = x1 - x0;
	int dy = y1 - y0;
	
	char ix = (dx > 0) - (dx < 0);
	dx = abs(dx) << 1;
	
	char iy = (dy > 0) - (dy < 0);
	dy = abs(dy) << 1;

	
	if ( dx >= dy )
	{
		int error = (dy - dx >> 1);
		
		while( x0 != x1 )
		{
			if ( error >= 0 )
			{
				if ( error || (ix > 0) )
				{
					y0 += iy;
					error -= dx;
				}
			}
			
			x0 += ix;
			error += dy;
			
			pixel = getPixelAt( rb, x0, y0 );
			pixel[0] = (cA[0]*255);
			pixel[1] = (cA[1]*255);
			pixel[2] = (cA[2]*255);		
			pixel[3] = (cA[3]*255);
		}
	}
	else
	{
		int error = (dx - dy >> 1);
		
		while( y0 != y1 )
		{
			if ( error >= 0 )
			{
				if ( error || (iy > 0) )
				{
					x0 += ix;
					error -= dy;
				}
			}
			
			y0 += iy;
			error += dx;
			
			pixel = getPixelAt( rb, x0, y0 );
			pixel[0] = (cA[0]*255);
			pixel[1] = (cA[1]*255);
			pixel[2] = (cA[2]*255);		
			pixel[3] = (cA[3]*255);
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
	clearDepth( rb, 256 );
	
	glm::vec3 lightposition( 0.0, -0.0, 0.25 );
	glm::mat3 normal_matrix = glm::inverseTranspose( glm::mat3( camera.modelview ) );
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
		
		// also transform nomals into eye space
		t->v[0].normal_eye = normal_matrix * t->v[0].normal;
		t->v[1].normal_eye = normal_matrix * t->v[1].normal;
		t->v[2].normal_eye = normal_matrix * t->v[2].normal;
		
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
	
	Line line;
	line.start.color = Color( 0, 255, 0, 255 );
	line.start.position = glm::vec3( 80, 75, 0 );
	line.end.position = glm::vec3( 150, 100, 0 );
	line.end.color = Color( 0, 0, 255, 255 );
	renderLine( rb, &line );
	
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
	fprintf( stdout, "\ttexture_mapping: %i\n", settings().texture_mapping );
	fprintf( stdout, "\tlighting: %i\n", settings().lighting );
	
	// setup camera
	Camera cam1;
	cam1.modelview = glm::translate( glm::mat4(1.0f), glm::vec3( 0, 0, -1.5 ));
	cam1.projection = glm::perspective( 60.0f, ((float)VIEWPORT_WIDTH/(float)VIEWPORT_HEIGHT), .1f, 2048.0f );
	cam1.modelview_projection = cam1.projection * cam1.modelview;

	
	cam1.modelview = glm::rotate( cam1.modelview, 45.0f, glm::vec3( 0, 1, 0 ) );

	RenderBuffer render_buffer;
	render_buffer.width = (int)VIEWPORT_WIDTH;
	render_buffer.height = (int)VIEWPORT_HEIGHT;
	render_buffer.channels = 4;
	render_buffer.pixels = new unsigned char[ render_buffer.width * render_buffer.height * render_buffer.channels ];
	render_buffer.zbuffer = new float[ render_buffer.width * render_buffer.height ];

	Texture tex;
	
	loadTexture( "../assets/checker2.png", tex );	

	num_triangles = 0;
	
	num_triangles = 2;
#if 1
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
#endif
	
	
	
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

	int save_result = SOIL_save_image( "../output/raster.bmp", SOIL_SAVE_TYPE_BMP, render_buffer.width, render_buffer.height, render_buffer.channels, render_buffer.pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}

	delete [] render_buffer.pixels;
	delete [] render_buffer.zbuffer;
	
	freeTexture( tex );

	return 0;
}
