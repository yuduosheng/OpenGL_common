#include "object.h"

object::object()
{
}
object::~object()
{
	vertex.clear();
	indices.clear();
}
void object::readFile(const char *filename)
{
	ifstream fileStream(filename, ios::in);
	glm::vec3 vData;
	glm::uvec3 iData;
	if (fileStream.is_open())
	{
		string head;
		while (fileStream >> head)
		{
			if (head == "v")
			{
				fileStream >> vData.x >> vData.y >> vData.z;
				vertex.push_back(vData);
				cout << vData.x <<" "<< vData.y << " " << vData.z <<endl;
			}
			else
				if (head == "f")
				{
				    fileStream >> iData.x >> iData.y >> iData.z;
					iData.x -= 1;
					iData.y -= 1;
					iData.z -= 1;
					indices.push_back(iData);
					cout << iData.x << " " << iData.y << " " << iData.z << endl;
				}
		}

		fileStream.close();
	}
}

void object::render()
{
	glGenBuffers(1, &mVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, mVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		vertex.size() * sizeof(glm::vec3),
		&vertex[0],
		GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &mIBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mIBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		indices.size() * sizeof(glm::uvec3),
		&indices[0],
		GL_STATIC_DRAW);

	glDrawElements(GL_TRIANGLES, indices.size() * 3.0f, GL_UNSIGNED_INT, 0);
}