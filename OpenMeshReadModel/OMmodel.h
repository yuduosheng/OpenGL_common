#include "App.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

#define	FlatShading 0x08
#define	SmoothShading 0x10

class OMmodel
{
private:
	MyMesh                  mesh;
	GLint                   meshVetexNum = 0;
	GLint                   meshFaceNum = 0;
	GLint                   meshHalfEdgeNum = 0;
	GLint                   meshBoundryEdgeNum = 0;
	GLint                   mSM = FlatShading;
	vector<OpenMesh::Vec3f> meshVertexBuffer;//vertex buffer
	vector<OpenMesh::Vec3f> meshVertexNormalBuffer;//vertex normal buffer
	vector<OpenMesh::Vec3f> meshFaceNormalBuffer;//face buffer normal


	GLuint                  meshVBuffer;
	GLuint                  meshNBuffer;
	GLuint                  meshFNormal;
public:
	OMmodel();
	~OMmodel();
	void PrintMeshStatus()
	{
		cout << "The number of vertices: " << meshVetexNum << endl;
		cout << "The number of faces: " << meshFaceNum << endl;
		cout << "The number of halfedges: " << meshHalfEdgeNum << endl;
		cout << "The number of boundry-edges: " << meshBoundryEdgeNum << endl;
	}
	void SetSM()
	{ 
		if (mSM == FlatShading)
			mSM = SmoothShading;
		else
			mSM = FlatShading;
	};
	bool OpenMeshReadFile(const char * filename);
	void RenderModel();
};