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







// ----------------------------------------
// Material

typedef glm::vec4 (*MaterialColor)( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id );

enum
{
	MATERIAL_LAMBERT = 1,	
};

struct Material
{
	int type;
	void * data;
	
	
	
	MaterialColor colorAtPoint;
};

Material * _materials = 0;
unsigned int _num_materials = 0;
unsigned int _current_material = 0;


struct LambertMaterial
{
	unsigned int color;
};

struct MirrorMaterial
{
	int id;
};


glm::vec4 LambertColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id );
glm::vec4 MirrorColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id );



Material * nextMaterial()
{
	if ( _current_material >= _num_materials )
	{
		fprintf( stdout, "Reached MAX materials! (%i)\n", _num_materials );
		return 0;
	}
	
	Material * m = &_materials[ _current_material++ ];
	return m;
}

void allocMaterials( unsigned int max_materials )
{
	_materials = (Material*)malloc( sizeof(Material) * max_materials );
	_current_material = 0;	
	_num_materials = max_materials;
	
	
	for( unsigned int m = 0; m < _num_materials; ++m )
	{
		Material * material = &_materials[ m ];
		material->type = 0;
		material->colorAtPoint = 0;
		material->data = 0;
	}
}

void purgeMaterials()
{
	for( unsigned int p = 0; p < _num_materials; ++p )
	{
		Material * material = &_materials[ p ];
		if ( material->type != 0 )
		{
			free( material->data );
		}
	}
	free(_materials);
	_num_materials = 0;
}


Material * createLambertMaterial( unsigned int color )
{
	LambertMaterial * mat = (LambertMaterial*)malloc( sizeof(LambertMaterial) );
	mat->color = color;
	
	
	Material * m = nextMaterial();
	if ( m )
	{
		m->data = mat;
		m->colorAtPoint = LambertColorAtPoint;
	}
	
	return m;
}

Material * createMirrorMaterial()
{
	MirrorMaterial * mat = (MirrorMaterial*)malloc( sizeof(MirrorMaterial) );

	Material * m = nextMaterial();
	if ( m )
	{
		m->data = mat;
		m->colorAtPoint = MirrorColorAtPoint;
	}
	return m;
}


// ----------------------------------------
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
	
	Material * material;
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


struct Light
{
	glm::vec3 color;
	glm::vec3 origin;
	glm::vec3 specular;
	int enabled;
	float intensity; // normalized intensity
};


const int MAX_LIGHTS = 2;
struct SceneContext
{
	glm::vec4 background_color;
	int current_trace_depth;
	int max_trace_depth;
	
	glm::vec4 ambient;
	Light lights[ MAX_LIGHTS ];
};

static SceneContext scene;

Primitive * _primitives = 0;
unsigned int _num_primitives = 0;
unsigned int _current_primitive = 0;



void rayTrace( const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t, Primitive ** closestPrimitive, Primitive * ignorePrimitive )
{
	if ( scene.current_trace_depth >= scene.max_trace_depth )
	{
		fprintf( stdout, "Reached max trace depth!\n" );
		return;
	}
	
	scene.current_trace_depth++;
	*closestPrimitive = 0;
	Primitive * p;
	float min_distance = FLT_MAX;
	
	glm::vec3 origin(FLT_EPSILON, FLT_EPSILON, FLT_EPSILON);
	origin += rayOrigin;
	
	for( unsigned int i = 0; i < _current_primitive; ++i )
	{
		p = &_primitives[ i ];
		
		if ( p->type == 0 )
			continue;
		
		// ignore these primitives
		if ( ignorePrimitive != 0 && p == ignorePrimitive )
			continue;
		
		if ( p->intersection( p->data, origin, rayDirection, t ) )
		{
			// find the closest primitive
			// filter against *t greater than FLT_EPSILON to remove some artifacts
			glm::vec3 origin_to_hit = (origin + rayDirection*(*t));
			origin_to_hit = origin_to_hit - origin;
			float distance = glm::length( origin_to_hit );
			if ( distance < min_distance && (*t) > FLT_EPSILON )
			{
				min_distance = distance;
				*closestPrimitive = p;
			}
		}
	}
	scene.current_trace_depth--;
}

Primitive * rayTracePrimitive( const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t, bool find_closest )
{
	// if find_closest is false, it will immediately return the first result it finds...
	Primitive * closestPrimitive = 0;
	Primitive * p;
	float min_distance = FLT_MAX;
	float closest_t = 0;
	float e = (FLT_EPSILON*500.0);
	glm::vec3 origin = rayOrigin + (rayDirection*e);
	for( unsigned int i = 0; i < _current_primitive; ++i )
	{
		p = &_primitives[ i ];
		
		if ( p->type == 0 )
			continue;
		
		float local_t;
		if ( p->intersection( p->data, origin, rayDirection, &local_t ) )
		{
			// find the closest primitive
			// filter against *t greater than FLT_EPSILON to remove some artifacts
			glm::vec3 origin_to_hit = (origin + rayDirection*local_t);
			origin_to_hit = origin_to_hit - origin;
			float distance = glm::length( origin_to_hit );
			if ( (distance < min_distance) && (local_t > FLT_EPSILON) )
			{
				min_distance = distance;
				closestPrimitive = p;
				closest_t = local_t;
			}
			
			if ( !find_closest )
				break;
		}
	}
	
	*t = closest_t;
	
	return closestPrimitive;
}

glm::vec4 rayTraceColor( const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection )
{
	if ( scene.current_trace_depth >= scene.max_trace_depth )
	{
		fprintf( stdout, "Reached max trace depth!\n" );
		return scene.background_color;
	}
	
	scene.current_trace_depth++;
	float t;
	Primitive * closestPrimitive = rayTracePrimitive( rayOrigin, rayDirection, &t, true );
	
	glm::vec4 fcolor(0,0,0,1);
	if ( closestPrimitive )
	{
		fcolor = scene.ambient;
		Light * light;
		for( int light_num = 0; light_num < MAX_LIGHTS; ++light_num )
		{
			light = &scene.lights[ light_num ];
			if ( !light->enabled )
				continue;
			glm::vec3 hit = (rayOrigin + rayDirection*t);
			glm::vec3 normal = closestPrimitive->normalAtPoint( closestPrimitive->data, hit );
			glm::vec4 color = closestPrimitive->material->colorAtPoint( closestPrimitive, hit, normal, rayDirection );			
			fcolor.r += color.r;
			fcolor.g += color.g;
			fcolor.b += color.b;
		}

		fcolor.a = 1.0f;
	}
	else
	{
		fcolor = scene.background_color;
	}

	scene.current_trace_depth--;
	return fcolor;	
}

glm::vec4 phongLighting( glm::vec4 diffuseColor, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id )
{
	glm::vec4 color;
	Light * light;
	float specularLevel = 0.75f;
	float specularPower = 75.0f;	
	for( int light_num = 0; light_num < MAX_LIGHTS; ++light_num )
	{
		light = &scene.lights[ light_num ];
		if ( !light->enabled )
			continue;
		
		glm::vec3 light_direction = glm::normalize(light->origin - point);
		
		float nDotL = glm::clamp( glm::dot(glm::normalize(normal), light_direction), 0.0f, 1.0f );
		
		// only send a shadow probe ray if the normal is facing the light
		Primitive * shadow_prim = 0;
		if ( nDotL > 0.0 )
		{
			float t;
			shadow_prim = rayTracePrimitive( point, light_direction, &t, false );
		}
		
		glm::vec3 R = glm::reflect( light_direction, normal );
		float i = glm::pow( glm::clamp( glm::dot(id,R), 0.0f, 1.0f ), specularPower );
		glm::vec3 light_color = light->color;
		
		color.r += light->intensity * (nDotL * light_color[0] * diffuseColor[0]) + i*light->specular[0]*specularLevel;
		color.g += light->intensity * (nDotL * light_color[1] * diffuseColor[1]) + i*light->specular[1]*specularLevel;
		color.b += light->intensity * (nDotL * light_color[2] * diffuseColor[2]) + i*light->specular[2]*specularLevel;
		
		// if the shadow probe found geometry; multiply to create a shadow
		if ( shadow_prim )
		{
			color.r *= 0.5;
			color.g *= 0.5;
			color.b *= 0.5;
		}
	}
	
	
	return color;
}


glm::vec4 LambertColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id )
{
	LambertMaterial * m = (LambertMaterial*)primitive->material->data;

	glm::vec4 diffuseColor;
	convertColor( m->color, &diffuseColor[0] );
	glm::vec4 color;
	
	color = phongLighting( diffuseColor, point, normal, id );
		
	color.a = 1;
	return color;
}

glm::vec4 MirrorColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id )
{
//	MirrorMaterial * m = (MirrorMaterial*)primitive->material->data;

	
	float t;
	Primitive * p = 0;
	glm::vec3 r = glm::normalize( glm::reflect(id, normal) );
	
	
	p = rayTracePrimitive( point, r, &t, true );
//	glm::vec4 diffuseColor = rayTraceColor( point, r );
//	return phongLighting( diffuseColor, point, normal, id );
	
	
//	return rayTraceColor( point, r );
	//rayTrace( point, r, &t, &p, primitive );

	if ( p )
	{
//		glm::vec3 reflect_hit = point + (t * r);
//		return p->material->colorAtPoint( p, reflect_hit, p->normalAtPoint( p->data, reflect_hit ), r );
		glm::vec4 diffuseColor = rayTraceColor( point, r );
		return diffuseColor;
//		return phongLighting( diffuseColor, point, normal, id );
		
	}

//	return glm::vec4(0,0,0,1);
	return scene.background_color;
}

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
		primitive->normalAtPoint = 0;
		primitive->data = 0;
		primitive->material = 0;
	}
}











// ----------------------------------------
// Plane
bool PlaneIntersection( void * obj, const glm::vec3 & rayOrigin, const glm::vec3 & rayDirection, float * t )
{
	PlanePrimitive * plane = (PlanePrimitive*)obj;
	return IntersectRayPlane( rayOrigin, rayDirection, plane->origin, plane->normal, t );
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

glm::vec3 SphereNormalAtPoint( void * obj, const glm::vec3 & point )
{
	SpherePrimitive * sphere = (SpherePrimitive*)obj;
	return point - sphere->origin;
}

// ----------------------------------------




Primitive * nextPrimitive()
{
	if ( _current_primitive >= _num_primitives )
	{
		fprintf( stdout, "Reached MAX primitives! (%i)\n", _num_primitives );
		return 0;
	}
	
	Primitive * p = &_primitives[ _current_primitive++ ];
	return p;
}

void addPlane( const glm::vec3 & planeOrigin, const glm::vec3 & planeNormal, Material * material )
{
	PlanePrimitive * plane = (PlanePrimitive*)malloc( sizeof(PlanePrimitive) );
	plane->normal = planeNormal;
	plane->origin = planeOrigin;


	Primitive * p = nextPrimitive();
	if ( p )
	{
		p->type = PRIMITIVE_PLANE;
		p->data = plane;
		p->intersection = PlaneIntersection;
		p->normalAtPoint = PlaneNormalAtPoint;
		p->material = material;
	}
}

void addSphere( const glm::vec3 & sphereOrigin, float radius, Material * material )
{
	SpherePrimitive * sphere = (SpherePrimitive*)malloc( sizeof(SpherePrimitive) );
	sphere->origin = sphereOrigin;
	sphere->radius = radius;
	
	Primitive * p = nextPrimitive();
	if ( p )
	{
		p->type = PRIMITIVE_SPHERE;
		p->data = sphere;
		p->intersection = SphereIntersection;
		p->normalAtPoint = SphereNormalAtPoint;
		p->material = material;
	}
}

// ----------------------------------------

void render_scene( glm::vec3 & eye, RenderBuffer & rb )
{
	glm::vec3 view( 0.0f, 0.0f, -1.0f );
	glm::vec3 screenU( 1, 0, 0 );
	glm::vec3 screenV( 0, 1, 0 );
	
	float aspect = (rb.width/(float)rb.height);

		
	allocPrimitives( 4 );
	allocMaterials( 16 );
	
	
	Material * redMaterial = createLambertMaterial( Color( 255, 0, 0, 255 ) );
	Material * greenMaterial = createLambertMaterial( Color( 0, 255, 0, 255 ) );
	Material * blueMaterial = createLambertMaterial( Color( 0, 0, 255, 255 ) );
	Material * whiteMaterial = createLambertMaterial( Color( 255, 255, 255, 255 ) );
	Material * yellowMaterial = createLambertMaterial( Color( 255, 255, 0, 255 ) );	
	Material * mirrorMaterial = createMirrorMaterial();
	
	addPlane( glm::vec3( 0.0f, -2.0f, 0.0f ), glm::vec3( 0.0f, 1.0f, 0.0f ), greenMaterial );
	
	addSphere( glm::vec3( 0, 0, -10.0f ), 1.0f, blueMaterial );
	addSphere( glm::vec3( -2.0f, 0.2f, -4.0f ), 1.0f, redMaterial );
	addSphere( glm::vec3( 1.8f, -0.6f, -2.0f), 1.0f, mirrorMaterial );
	
	glm::vec3 rayOrigin = eye;
	unsigned char * pixels = rb.pixels;
	
	float objectColor[4];
	
	Primitive * prim = 0;
	float t;
	glm::vec3 normal;
	glm::vec3 lightdir;
	
	for( int h = rb.height; h > 0; --h )
//	for( int h = 0; h < rb.height; ++h )
	{
		for( int w = 0; w < rb.width; ++w )
		{
			// we need u and v to be [0,1] inclusive
			glm::vec3 u = screenU * (((float)w / rb.width) * 2.0f - 1.0f) * aspect;
			glm::vec3 v = screenV * (((float)h / rb.height) * 2.0f - 1.0f);		
			glm::vec3 rayDirection = glm::normalize(view + u + v);

#if 0
			rayTrace( rayOrigin, rayDirection, &t, &prim, 0 );
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
				rayTrace( intersection, lightdir, &st, &shadow_primitive, prim );
				if ( shadow_primitive )
				{
					glm::vec3 shadow_int = intersection + lightdir*st;
					
					// use ambient color as shadow
					objectColor[0] = ambient[0];
					objectColor[1] = ambient[1];
					objectColor[2] = ambient[2];
				}
				else
				{
					convertColor( prim->material->colorAtPoint( prim, intersection, normal, rayDirection ), objectColor );					
					ndl = glm::clamp( glm::dot( lightdir, normal ), 0.0f, 1.0f );
				}
			
				// diffuse = (ndl * lightColor) * objectColor;
//				finalColor[0] += (ndl * lightColor[0]) * objectColor[0];
//				finalColor[1] += (ndl * lightColor[1]) * objectColor[1];
//				finalColor[2] += (ndl * lightColor[2]) * objectColor[2];
				
				// phong = diffuse * (L dot N) + specular * (V dot R) * N
				glm::vec3 R = glm::reflect( lightdir, normal );				
				float i = glm::pow( glm::clamp( glm::dot(rayDirection,R), 0.0f, 1.0f ), 80.0f);
				
				float n = 1.0f;
				finalColor[0] += (ndl * lightColor[0]) * objectColor[0] + i*lightSpecular[0]*n;
				finalColor[1] += (ndl * lightColor[1]) * objectColor[1] + i*lightSpecular[1]*n;
				finalColor[2] += (ndl * lightColor[2]) * objectColor[2] + i*lightSpecular[2]*n;
				
			}

			// clamp colors
			if ( finalColor[0] > 1 )
				finalColor[0] = 1;
			if ( finalColor[1] > 1 )
				finalColor[1] = 1;
			if ( finalColor[2] > 1 )
				finalColor[2] = 1;
#endif
			
			
			glm::vec4 color = rayTraceColor( rayOrigin, rayDirection );

			if ( color[0] > 1 )
				color[0] = 1;
			if ( color[1] > 1 )
				color[1] = 1;
			if ( color[2] > 1 )
				color[2] = 1;
			
			pixels[0] = color[0] * 255.0f;
			pixels[1] = color[1] * 255.0f;
			pixels[2] = color[2] * 255.0f;			

//			pixels[0] = finalColor[0] * 255.0f;
//			pixels[1] = finalColor[1] * 255.0f;
//			pixels[2] = finalColor[2] * 255.0f;
			pixels[3] = 255;			
			pixels += rb.channels;
		}
	}
	
	
	purgePrimitives();
	purgeMaterials();
}


int main( int argc, char ** argv )
{	
	// setup a viewport
	Viewport vp;
	vp.offset = glm::vec2( 0, 0 );
	vp.size = glm::vec2( 512, 512 );

	scene.background_color = glm::vec4( 0.5f, 0.5f, 0.5f, 1.0f );
	scene.current_trace_depth = 0;
	scene.max_trace_depth = 8;
	
	// setup lights in the scene
	scene.ambient = glm::vec4( 0.10f, 0.10f, 0.125f, 1.0f );
	Light * light;
	
	light = &scene.lights[0];
	light->enabled = 1;
	light->intensity = 0.5f;
	light->color = glm::vec3( 1, 1, 1 );
	light->origin = glm::vec3( 0, 5, 0 );
	light->specular = glm::vec3( 1.0f, 1.0f, 1.0f );
	
	light = &scene.lights[1];
	light->enabled = 1;
	light->intensity = 0.5f;
	light->color = glm::vec3( 1.0f, 0.0f, 0.0f );
	light->origin = glm::vec3( 5, 5, -5 );
	light->specular = glm::vec3( 1.0f, 0.0f, 0.0f );
	
	
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
	glm::vec3 eye( 0, 0, 1 );
	render_scene( eye, render_buffer );
	
	int save_result = SOIL_save_image( "output/raytracer.bmp", SOIL_SAVE_TYPE_BMP, render_buffer.width, render_buffer.height, render_buffer.channels, render_buffer.pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}

	delete [] render_buffer.pixels;
	delete [] render_buffer.zbuffer;

	return 0;
}
