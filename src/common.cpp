#include <common.hpp>


// equation from wikipedia: http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
void barycentric_coords( const glm::vec2 & p, const glm::vec2 & a, const glm::vec2 & b, const glm::vec2 & c, float * alpha, float * beta, float * gamma )
{
	float detT = 1.0 / ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	*alpha = ((b.y - c.y)*(p.x - c.x) + (c.x - b.x) * (p.y - c.y)) * detT;
	*beta = ((c.y - a.y)*(p.x - c.x) + (a.x - c.x) * (p.y - c.y)) * detT;

	*gamma = 1.0f - *alpha - *beta;
}

