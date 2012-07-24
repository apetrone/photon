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


struct PhongMaterial
{
	unsigned int color;
	float specular_power;
	float specular_level;
};

struct MirrorMaterial
{
	unsigned int color;
	float reflectance;
	float specular_power;
	float specular_level;	
};


struct RefractiveMaterial
{
	unsigned int color;
	float refractive_index;
	float refractance;
	float specular_power;
	float specular_level;	
};


glm::vec4 LambertColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id );
glm::vec4 MirrorColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id );
glm::vec4 RefractiveColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id );


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


Material * createPhongMaterial( unsigned int color, float specularPower, float specularLevel )
{
	PhongMaterial * mat = (PhongMaterial*)malloc( sizeof(PhongMaterial) );
	mat->color = color;
	mat->specular_power = specularPower;
	mat->specular_level = specularLevel;
	
	
	Material * m = nextMaterial();
	if ( m )
	{
		m->data = mat;
		m->colorAtPoint = LambertColorAtPoint;
	}
	
	return m;
}

Material * createMirrorMaterial( unsigned int color, float reflectance, float specularPower, float specularLevel )
{
	MirrorMaterial * mat = (MirrorMaterial*)malloc( sizeof(MirrorMaterial) );
	mat->color = color;
	mat->reflectance = reflectance;
	mat->specular_power = specularPower;
	mat->specular_level = specularLevel;

	Material * m = nextMaterial();
	if ( m )
	{
		m->data = mat;
		m->colorAtPoint = MirrorColorAtPoint;
	}
	return m;
}

Material * createRefractiveMaterial( unsigned int color, float refractive_index )
{
	RefractiveMaterial * mat = (RefractiveMaterial*)malloc( sizeof(RefractiveMaterial) );
	mat->color = color;
	mat->refractive_index = refractive_index;

	Material * m = nextMaterial();
	if ( m )
	{
		m->data = mat;
		m->colorAtPoint = RefractiveColorAtPoint;
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

glm::vec4 phongLighting( glm::vec4 diffuseColor, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id, float specularPower, float specularLevel )
{
	glm::vec4 color;
	Light * light;	
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
	PhongMaterial * m = (PhongMaterial*)primitive->material->data;

	glm::vec4 diffuseColor;
	convertColor( m->color, &diffuseColor[0] );
	glm::vec4 color;
	

	color = phongLighting( diffuseColor, point, normal, id, m->specular_power, m->specular_level );
		
	color.a = 1;
	return color;
}

glm::vec4 MirrorColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id )
{
	MirrorMaterial * m = (MirrorMaterial*)primitive->material->data;
	
	float t;
	Primitive * p = 0;
	glm::vec3 r = glm::normalize( glm::reflect(id, normal) );
	
	glm::vec4 matColor;

	convertColor( m->color, &matColor[0] );

	// probe and see if the reflected vector hits anything
	p = rayTracePrimitive( point, r, &t, true );
	glm::vec4 diffuseColor;

	float rdn = 0;
	// it does, we'll add that object's diffuse color
	if ( p )
	{
		diffuseColor = rayTraceColor( point, r );

		// scale the reflected color based on the angle between reflection and surface normal
		rdn = m->reflectance * fmax( 0.0f, glm::dot(r, glm::normalize(normal)) );
	}

	glm::vec4 myColor = phongLighting( matColor, point, normal, id, m->specular_power, m->specular_level );

	glm::vec4 outColor = myColor + (rdn * diffuseColor);


	outColor.a = 1;
	return outColor;
}


glm::vec4 RefractiveColorAtPoint( struct Primitive * primitive, const glm::vec3 & point, const glm::vec3 & normal, const glm::vec3 & id )
{
	RefractiveMaterial * m = (RefractiveMaterial*)primitive->material->data;
	
	float t;
	Primitive * p = 0;
	glm::vec3 r = glm::normalize( glm::reflect(id, normal) );
	
	glm::vec4 matColor;

	convertColor( m->color, &matColor[0] );

	// probe and see if the reflected vector hits anything
	p = rayTracePrimitive( point, r, &t, true );
	glm::vec4 diffuseColor;

	float rdn = 0;
	// it does, we'll add that object's diffuse color
	if ( p )
	{
		diffuseColor = rayTraceColor( point, r );

		// scale the reflected color based on the angle between reflection and surface normal
		rdn = 0.5f * fmax( 0.0f, glm::dot(r, glm::normalize(normal)) );
	}

	glm::vec4 myColor = phongLighting( matColor, point, normal, id, 75.0f, 0.75f );

	glm::vec4 outColor = myColor + (rdn * diffuseColor);


	outColor.a = 1;
	return outColor;
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

		
	allocPrimitives( 8 );
	allocMaterials( 16 );
	
	
	Material * redMaterial = createPhongMaterial( Color( 255, 0, 0, 255 ), 75.0f, 0.75f );
	Material * greenMaterial = createPhongMaterial( Color( 0, 255, 0, 255 ), 75.0f, 0.75f );
	Material * blueMaterial = createPhongMaterial( Color( 0, 0, 255, 255 ), 75.0f, 0.75f );
	Material * whiteMaterial = createPhongMaterial( Color( 255, 255, 255, 255 ), 75.0f, 0.75f );
	Material * yellowMaterial = createPhongMaterial( Color( 255, 255, 0, 255 ), 75.0f, 0.75f );	
	Material * mirrorMaterial = createMirrorMaterial( Color(0, 0, 0, 255 ), 0.5f, 100.0f, 1.0f );
	Material * refraction = createRefractiveMaterial( Color( 0, 0, 0, 255), 1.30 );
	
	addPlane( glm::vec3( 0.0f, -2.0f, 0.0f ), glm::vec3( 0.0f, 1.0f, 0.0f ), greenMaterial );
	
	addSphere( glm::vec3( 0, 0, -10.0f ), 1.0f, blueMaterial );
	addSphere( glm::vec3( -2.0f, 0.2f, -4.0f ), 1.0f, redMaterial );
	addSphere( glm::vec3( 1.8f, -0.6f, -2.0f), 1.0f, mirrorMaterial );
	addSphere( glm::vec3( -2.75f, -1.0f, -1.0f), 1.0f, yellowMaterial );
	
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
