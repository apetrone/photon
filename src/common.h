/*
	Copyright (c) 2013, <Adam Petrone>
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:
		* Redistributions of source code must retain the above copyright
			notice, this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright
			notice, this list of conditions and the following disclaimer in the
			documentation and/or other materials provided with the distribution.
		* Neither the name of the <organization> nor the
			names of its contributors may be used to endorse or promote products
			derived from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
	DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
	ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#pragma once

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
	glm::vec3 normal_eye;
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

struct Triangle
{
	Vertex v[3];
	Texture * texture;
};


void barycentric_coords( const glm::vec2 & p, const glm::vec2 & a, const glm::vec2 & b, const glm::vec2 & c, float * alpha, float * beta, float * gamma );


inline void convertColor( unsigned int color, float out[4] )
{
	out[3] = ((color >> 24) & 0xFF) / 255.0;
	out[2] = ((color >> 16) & 0xFF) / 255.0;
	out[1] = ((color >> 8) & 0xFF) / 255.0;
	out[0] = (color & 0xFF) / 255.0;
}

inline unsigned int Color( unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a )
{
	return (_r | _g << 8 | _b << 16 | _a << 24 );
}

inline unsigned int scaleColor( unsigned int a, float t )
{
	float c[4];
	convertColor( a, c );
	
	c[0] *= t;
	c[1] *= t;
	c[2] *= t;
	c[3] = 1.0;
	
//	float ac[] = { (a & 0xFF), ((a >> 8) & 0xFF), ((a >> 16) & 0xFF), ( (a >> 24) & 0xFF) };
	
	return (int(255*c[0]) ) | (int(255*c[1]) << 8) | (int(255*c[2]) << 16) | (int(255) << 24);
}

void freeTexture( Texture & texture );
void loadTexture( const char * filename, Texture & texture );
void writeTexture( const char * filename, Texture & texture );

