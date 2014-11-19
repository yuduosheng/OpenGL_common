#include "object.h"

object::object()
{
	mIndiceNum = 0;
}
object::~object()
{
	vertex.clear();
	indices.clear();
}
void object::readFile(const char *filename)
{//read obj file which the number of vertex of a primitive is 3
	ifstream fileStream(filename, ios::in);
	glm::vec3 vData;
	glm::uint iData;
	if (fileStream.is_open())
	{
		string head;
		while (fileStream >> head)
		{
			if (head == "v")
			{
				fileStream >> vData.x >> vData.y >> vData.z;
				vertex.push_back(vData);
				//cout << vData.x <<" "<< vData.y << " " << vData.z <<endl;
			}
			else
				if (head == "f")
				{
				    mIndiceNum++;
				    fileStream >> iData;
					iData -= 1;
					indices.push_back(iData);
					fileStream >> iData;
					iData -= 1;
					indices.push_back(iData);
					fileStream >> iData;
					iData -= 1;
					indices.push_back(iData);
				}
		}

		fileStream.close();
	}
}
void object::readFileIndice4(const char *filename)
{//read obj file which the number of vertex of a primitive is 4
	ifstream fileStream(filename, ios::in);
	glm::vec3 vData;
	glm::uint iData;
	if (fileStream.is_open())
	{
		string head;
		while (fileStream >> head)
		{
			if (head == "v")
			{
				fileStream >> vData.x >> vData.y >> vData.z;
				vertex.push_back(vData);
				//cout << vData.x << " " << vData.y << " " << vData.z << endl;
			}
			else
				if (head == "f")
				{
				    mIndiceNum++;
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
					fileStream >> iData;
					iData -= 1;
					indices.push_back(iData);
				}
		}

		fileStream.close();
	}
}
void object::readFileIndice5(const char *filename)
{//read obj file which the number of vertex of a primitive is 5
	ifstream fileStream(filename, ios::in);
	glm::vec3 vData;
	glm::uint iData;
	if (fileStream.is_open())
	{
		string head;
		while (fileStream >> head)
		{
			if (head == "v")
			{
				fileStream >> vData.x >> vData.y >> vData.z;
				vertex.push_back(vData);
				//cout << vData.x << " " << vData.y << " " << vData.z << endl;
			}
			else
				if (head == "f")
				{
				    mIndiceNum++;
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
				    fileStream >> iData;
				    iData -= 1;
				    indices.push_back(iData);
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
		indices.size() * sizeof(uint),
		&indices[0],
		GL_STATIC_DRAW);

	glDrawElements(GL_TRIANGLES, mIndiceNum * 3.0f, GL_UNSIGNED_INT, 0);
}
void object::renderQuad()
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
		indices.size() * sizeof(uint),
		&indices[0],
		GL_STATIC_DRAW);

	glDrawElements(GL_QUADS, mIndiceNum * 4.0f, GL_UNSIGNED_INT, 0);
}
void object::renderPolygon()
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
		indices.size() * sizeof(uint),
		&indices[0],
		GL_STATIC_DRAW);

	glDrawElements(GL_POLYGON, mIndiceNum * 5.0f, GL_UNSIGNED_INT, 0);
}