#include "App.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

#define	FlatShading 0x08
#define	SmoothShading 0x10
#define	Valence 0x01
#define	MeanCurvature 0x02
#define	GaussianCurvature 0x04
class OMmodel
{
private:
	MyMesh                  mesh;
	GLint                   meshVetexNum = 0;
	GLint                   meshFaceNum = 0;
	GLint                   meshHalfEdgeNum = 0;
	GLint                   meshBoundryEdgeNum = 0;
	GLint                   mSM = FlatShading;
	GLint                   mCM = Valence;
	vector<OpenMesh::Vec3f> meshVertexBuffer;//vertex buffer
	vector<OpenMesh::Vec3f> meshVertexNormalBuffer;//vertex normal buffer
	vector<OpenMesh::Vec3f> meshFaceNormalBuffer;//face buffer normal
	vector<OpenMesh::Vec3f> meshVertexColorBuffer;//vertex color
	vector<OpenMesh::Vec3f> meshCurColorBuffer;//mean curvatrue color
	vector<OpenMesh::Vec3f> meshGCurColorBuffer;//Gaussian curvatrue color
	GLuint                  meshVBuffer;
	GLuint                  meshNBuffer;
	GLuint                  meshFNormal;
	GLuint                  meshVColor;
public:
	OMmodel();
	~OMmodel();
	void PrintMeshStatus()
	{
		std::cout << "The number of vertices: " << meshVetexNum << endl;
		std::cout << "The number of faces: " << meshFaceNum << endl;
		std::cout << "The number of halfedges: " << meshHalfEdgeNum << endl;
		std::cout << "The number of boundry-edges: " << meshBoundryEdgeNum << endl;
	}
	void SetSM()
	{ 
		if (mSM == FlatShading)
			mSM = SmoothShading;
		else
			mSM = FlatShading;
	};
	void SetCM()
	{
		if (mCM == Valence)
			mCM = MeanCurvature;
		else
			if (mCM == MeanCurvature)
				mCM = GaussianCurvature;
			else
				mCM = Valence;
	};
	void OpenMeshReadFile(const char * filename);
	void RenderModel();
	void RenderModelWithColor();
};