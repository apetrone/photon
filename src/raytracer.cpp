#include <common.hpp>

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
}

// http://wiki.cgsociety.org/index.php/Ray_Sphere_Intersection
// this accepts ray origin and sphere origin which it will calculate the new origin from
bool IntersectRaySphere( const glm::vec3 rayOrigin, const glm::vec3 rayDirection, const glm::vec3 sphereOrigin, float sphereRadius, float * t )
{
	glm::vec3 v = rayOrigin - sphereOrigin;
	float a = glm::dot( rayDirection, rayDirection );
	float b = 2.0f * glm::dot( rayDirection, v );
	float c = glm::dot( v, v ) - (sphereRadius*sphereRadius);
	
	float det = b * b - 4.0f * a * c;

	// miss! no intersections
	if ( det < 0.0f )
	{
		return false;
	}
	
	float detsqrt = sqrt( det );


#if 0
	float t0, t1;	
	float q;
	
	if ( b < 0 )
	{
		q = ( -b - detsqrt ) * 0.5f;
	}
	else
	{
		q = ( -b + detsqrt ) * 0.5f;
	}

	t0 = q / a;
	t1 = c / q;
	
	// make sure t0 is smaller than t1
	if ( t0 > t1 )
	{
		// swap the two
		float temp = t0;
		t0 = t1;
		t1 = temp;
	}
	
	// if t1 < 0, object is in the ray's negative direction (miss)
	if ( t1 < 0 )
	{
		return false;
	}
	
	// if t0  < 0, intersection point is at t1
	fprintf( stdout, "t0: %g, t1: %g\n", t0, t1 );
	if ( t0 < 0 )
	{
		*t = t1;
	}
	else
	{
		*t = t0;
	}
#elif 0
	float rootMin = (-b - detsqrt) / (2.0f * a);
	float rootMax = (-b + detsqrt) / (2.0f * a);	
	
	
	if ( rootMax < 0.0f )
	{
		return false;
	}
	else
	{
		if ( rootMin < 0.0f )
		{
			*t = rootMax;
		}
		else
		{
			*t = rootMin;
		}
	}
#else
		float dn = 2*a;
		
		float td = (-b - detsqrt) / dn;
		if ( td > FLT_EPSILON )
		{
			*t = td;
			return true;
		}
	
		td = (-b + detsqrt) / dn;
	
		if ( td > FLT_EPSILON )
		{
			*t = td;
			return true;
		}
	
#endif

	return true;
}

bool IntersectRayPlane( const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, const glm::vec3 planeOrigin, const glm::vec3 & planeNormal, float * t )
{
	float ldn = glm::dot( rayDirection, planeNormal );
	
	if ( fabsf(ldn) < FLT_EPSILON )
		return false;
	
	float org = glm::dot( (planeOrigin - rayOrigin), planeNormal );

	float planeIntersection = (org / ldn);
	if ( planeIntersection < 0.0f )
		return false;
	
	*t = planeIntersection;
//	fprintf( stdout, "org: %g, ldn: %g, t: %g\n", org, ldn, t );
	
	return true;
} // IntersectRayPlane


typedef bool (*PrimitiveIntersection)( void * obj, const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t );
typedef unsigned int (*PrimitiveColor)( void * obj, const glm::vec3 & point );
typedef glm::vec3 (*PrimitiveNormal)( void * obj, const glm::vec3 & point );

enum 
{
	PRIMITIVE_PLANE = 1,
	PRIMITIVE_SPHERE = 2,
};

struct Primitive
{
	int type;
	PrimitiveIntersection intersection;
	PrimitiveColor colorAtPoint;
	PrimitiveNormal normalAtPoint;
	void * data;
};

struct PlanePrimitive
{
	glm::vec3 origin;
	glm::vec3 normal;
	unsigned int color;
};

struct SpherePrimitive
{
	glm::vec3 origin;	
	float radius;
	unsigned int color;
};


Primitive * _primitives = 0;
unsigned int _num_primitives = 0;
unsigned int _current_primitive = 0;

// ----------------------------------------
void purgePrimitives()
{
	for( unsigned int p = 0; p < _num_primitives; ++p )
	{
		Primitive * primitive = &_primitives[ p ];
		if ( primitive->type != 0 )
		{
			free(primitive->data);
		}
	}
	free(_primitives);
	_num_primitives = 0;
}

void allocPrimitives( unsigned int total_primitives )
{
	_primitives = (Primitive*)malloc( total_primitives * sizeof(Primitive) );
	_num_primitives = total_primitives;
	_current_primitive = 0;
	
	for( unsigned int p = 0; p < _num_primitives; ++p )
	{
		Primitive * primitive = &_primitives[ p ];
		primitive->type = 0;
		primitive->colorAtPoint = 0;
		primitive->normalAtPoint = 0;
		primitive->data = 0;
	}
}











// ----------------------------------------
// Plane
bool PlaneIntersection( void * obj, const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t )
{
	PlanePrimitive * plane = (PlanePrimitive*)obj;
	return IntersectRayPlane( rayOrigin, rayDirection, plane->origin, plane->normal, t );
}

unsigned int PlaneColorAtPoint( void * obj, const glm::vec3 & point )
{
	PlanePrimitive * plane = (PlanePrimitive*)obj;
	return plane->color;
}

glm::vec3 PlaneNormalAtPoint( void * obj, const glm::vec3 & point )
{
	PlanePrimitive * plane = (PlanePrimitive*)obj;
	return plane->normal;
}

// ----------------------------------------
// Sphere
bool SphereIntersection( void * obj, const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t )
{
	SpherePrimitive * sphere = (SpherePrimitive*)obj;
	return IntersectRaySphere( rayOrigin, rayDirection, sphere->origin, sphere->radius, t );
}

unsigned SphereColorAtPoint( void * obj, const glm::vec3 & point )
{
	SpherePrimitive * sphere = (SpherePrimitive*)obj;
	return sphere->color;
}

glm::vec3 SphereNormalAtPoint( void * obj, const glm::vec3 & point )
{
	SpherePrimitive * sphere = (SpherePrimitive*)obj;
	return point - sphere->origin;
}

// ----------------------------------------


void traceScene( const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t, Primitive ** closestPrimitive, Primitive * ignorePrimitive )
{
	*closestPrimitive = 0;
	Primitive * p;
	for( unsigned int i = 0; i < _current_primitive; ++i )
	{
		p = &_primitives[ i ];
		
		if ( p->type == 0 )
			continue;
		
		// ignore these primitives
		if ( ignorePrimitive != 0 && p == ignorePrimitive )
			continue;
		
		if ( p->intersection( p->data, rayOrigin, rayDirection, t ) )
		{
			*closestPrimitive = p;
		}
	}
}


Primitive * nextPrimitive()
{
	Primitive * p = &_primitives[ _current_primitive++ ];
	return p;
}

void addPlane( const glm::vec3 & planeOrigin, const glm::vec3 & planeNormal, unsigned int color )
{
	PlanePrimitive * plane = (PlanePrimitive*)malloc( sizeof(PlanePrimitive) );
	plane->normal = planeNormal;
	plane->origin = planeOrigin;
	plane->color = color;

	Primitive * p = nextPrimitive();
	p->type = PRIMITIVE_PLANE;
	p->data = plane;
	p->intersection = PlaneIntersection;
	p->colorAtPoint = PlaneColorAtPoint;
	p->normalAtPoint = PlaneNormalAtPoint;
}

void addSphere( const glm::vec3 & sphereOrigin, float radius, unsigned int color )
{
	SpherePrimitive * sphere = (SpherePrimitive*)malloc( sizeof(SpherePrimitive) );
	sphere->origin = sphereOrigin;
	sphere->radius = radius;
	sphere->color = color;
	
	Primitive * p = nextPrimitive();
	p->type = PRIMITIVE_SPHERE;
	p->data = sphere;
	p->intersection = SphereIntersection;
	p->colorAtPoint = SphereColorAtPoint;
	p->normalAtPoint = SphereNormalAtPoint;
}

// ----------------------------------------

void raytraceScene( glm::vec3 & eye, RenderBuffer & rb )
{
	glm::vec3 view( 0.0f, 0.0f, -1.0f );
	glm::vec3 screenU( 1, 0, 0 );
	glm::vec3 screenV( 0, 1, 0 );
	
	float aspect = (rb.width/(float)rb.height);
	glm::vec3 lightPosition( 0.0f, 5.0f, 0.0f );

	glm::vec3 sphereOrigin[3] = { glm::vec3( 0, 0, -10.0f ), glm::vec3( -2.3f, 0.2f, -6.0f ), glm::vec3( 1.8f, -0.6f, -3.0f) };
	float sphereRadius[3] = { 1.0f, 1.0f, 1.0f };
	unsigned int sphereColor[3] = { Color( 0, 0, 255, 255 ), Color( 255, 0, 0, 255 ), Color( 255, 255, 255, 255 ) };
	
	glm::vec3 ambient( 0.10f, 0.10f, 0.125f );
	glm::vec3 lightColor( 1.0f, 1.0f, 1.0f );
	
	
	glm::vec3 planeNormal( 0.0f, 1.0f, 0.0f );
	glm::vec3 planeOrigin( 0.0f, -2.0f, 0.0f );
	unsigned int planeColor = Color( 255, 128, 0, 255 );
	
	
	allocPrimitives( 4 );
	addPlane( planeOrigin, planeNormal, planeColor );
	
	addSphere( sphereOrigin[0], sphereRadius[0], sphereColor[0] );
	addSphere( sphereOrigin[1], sphereRadius[1], sphereColor[1] );
	addSphere( sphereOrigin[2], sphereRadius[2], sphereColor[2] );	
	
	glm::vec3 rayOrigin = eye;
	unsigned char * pixels = rb.pixels;
	
	float objectColor[4];
	
	Primitive * prim = 0;
	float t;
	glm::vec3 normal;
	glm::vec3 lightdir;
	
	float backgroundColor[4];
	convertColor( Color( 128, 128, 128, 255 ), backgroundColor );
	
	for( int h = rb.height; h > 0; --h )
//	for( int h = 0; h < rb.height; ++h )
	{
		for( int w = 0; w < rb.width; ++w )
		{
			// we need u and v to be [0,1] inclusive
			glm::vec3 u = screenU * (((float)w / rb.width) * 2.0f - 1.0f) * aspect;
			glm::vec3 v = screenV * (((float)h / rb.height) * 2.0f - 1.0f);		
			glm::vec3 rayDirection = glm::normalize(view + u + v);

			traceScene( rayOrigin, rayDirection, &t, &prim, 0 );
			glm::vec3 intersection = rayOrigin + (t*rayDirection);
			float finalColor[3] = { backgroundColor[0], backgroundColor[1], backgroundColor[2] };
			if ( prim )
			{
				lightdir = glm::normalize( lightPosition - intersection );
				normal = glm::normalize( prim->normalAtPoint( prim->data, intersection ) );
				
				float ndl = 1.0f;
				
				finalColor[0] = ambient[0];
				finalColor[1] = ambient[1];
				finalColor[2] = ambient[2];
				
				Primitive * shadow_primitive = 0;
				float st;
				traceScene( intersection, lightdir, &st, &shadow_primitive, prim );
				if ( shadow_primitive )
				{
					glm::vec3 shadow_int = intersection + lightdir*st;
					// use the blocker's color to color this pixel
					convertColor( shadow_primitive->colorAtPoint( shadow_primitive->data, shadow_int), objectColor );
					
					// use ambient color as shadow
					objectColor[0] = ambient[0];
					objectColor[1] = ambient[1];
					objectColor[2] = ambient[2];
				}
				else
				{
					Primitive * reflected = 0;

					// calculate reflected vector
//					glm::vec3 r = glm::reflect( rayDirection, normal );
//					traceScene( intersection, r, &st, &reflected, prim );
					if ( reflected )
					{
						convertColor( reflected->colorAtPoint( reflected->data, intersection ), objectColor );
					}
					else
					{
						convertColor( prim->colorAtPoint( prim->data, intersection ), objectColor );
					}
					
					ndl = fmax( 0.0f, glm::dot( lightdir, normal ) );
				}
			
				
				finalColor[0] += (ndl * lightColor[0]) * objectColor[0];
				finalColor[1] += (ndl * lightColor[1]) * objectColor[1];
				finalColor[2] += (ndl * lightColor[2]) * objectColor[2];
			}

			// clamp colors
			if ( finalColor[0] > 1 )
				finalColor[0] = 1;
			if ( finalColor[1] > 1 )
				finalColor[1] = 1;
			if ( finalColor[2] > 1 )
				finalColor[2] = 1;
			
			pixels[0] = finalColor[0] * 255.0f;
			pixels[1] = finalColor[1] * 255.0f;
			pixels[2] = finalColor[2] * 255.0f;
			pixels[3] = 255;			
			pixels += rb.channels;
		}
	}
	
	
	purgePrimitives();
}


int main( int argc, char ** argv )
{	
	// setup a viewport
	Viewport vp;
	vp.offset = glm::vec2( 0, 0 );
	vp.size = glm::vec2( 512, 512 );

	// print_settings();
//	fprintf( stdout, "======== Settings ========\n" );
//	fprintf( stdout, "\tperpsective_correct: %i\n", settings().perspective_correct );
//	fprintf( stdout, "\tbackface_culling: %i\n", settings().backface_culling );
//	fprintf( stdout, "\ttexture_mapping: %i\n", settings().texture_mapping );
//	fprintf( stdout, "\tlighting: %i\n", settings().lighting );
	
	// setup camera
//	Camera cam1;
//	cam1.modelview = glm::translate( glm::mat4(1.0f), glm::vec3( 0, 0, -1.5 ));
//	cam1.projection = glm::perspective( 60.0f, ((float)VIEWPORT_WIDTH/(float)VIEWPORT_HEIGHT), .1f, 2048.0f );
//	cam1.modelview_projection = cam1.projection * cam1.modelview;

	RenderBuffer render_buffer;
	render_buffer.width = (int)vp.size.x;
	render_buffer.height = (int)vp.size.y;
	render_buffer.channels = 4;
	render_buffer.pixels = new unsigned char[ render_buffer.width * render_buffer.height * render_buffer.channels ];
	render_buffer.zbuffer = new float[ render_buffer.width * render_buffer.height ];
	

	fprintf( stdout, "======== Render Scene ========\n" );
	
//	clearColor( render_buffer, Color(128, 128, 128, 255) );
	glm::vec3 eye( 0, 0, 1 );
	raytraceScene( eye, render_buffer );
	
	int save_result = SOIL_save_image( "output/raytracer.bmp", SOIL_SAVE_TYPE_BMP, render_buffer.width, render_buffer.height, render_buffer.channels, render_buffer.pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}

	delete [] render_buffer.pixels;
	delete [] render_buffer.zbuffer;

	return 0;
}
