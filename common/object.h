#ifndef OBJECT_H
#define OBJECT_H

#include "App.h"

#include <vector>
#include <glm\glm.hpp>
class object
{
public:
	object();
	~object();

	void readFile(const char *filename);
	void render();
protected:
	GLuint                  mVBuffer;
	GLuint                  mIBuffer;
	vector<glm::vec3>       vertex;
	vector<glm::uvec3>      indices;
};


#endif