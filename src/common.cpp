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
#include "common.h"


// equation from wikipedia: http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
void barycentric_coords( const glm::vec2 & p, const glm::vec2 & a, const glm::vec2 & b, const glm::vec2 & c, float * alpha, float * beta, float * gamma )
{
	float detT = 1.0 / ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	*alpha = ((b.y - c.y)*(p.x - c.x) + (c.x - b.x) * (p.y - c.y)) * detT;
	*beta = ((c.y - a.y)*(p.x - c.x) + (a.x - c.x) * (p.y - c.y)) * detT;

	*gamma = 1.0f - *alpha - *beta;
}

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


void writeTexture( const char * filename, Texture & texture )
{
	int save_result = SOIL_save_image( filename, SOIL_SAVE_TYPE_TGA, texture.width, texture.height, texture.channels, texture.pixels );
	if ( save_result != 1 )
	{
		fprintf( stdout, "saving to output failed!\n" );
	}	
}