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
	void readFileIndice4(const char *filename);
	void readFileIndice5(const char *filename);
	void render();
	void renderQuad();
	void renderPolygon();
protected:
	GLuint                  mVBuffer;
	GLuint                  mIBuffer;
	vector<vec3>            vertex;
	vector<uint>            indices;
	uint                    mIndiceNum;
};


#endif