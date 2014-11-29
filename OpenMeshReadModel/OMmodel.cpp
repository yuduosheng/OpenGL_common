#include "OMmodel.h"

OMmodel::OMmodel()
{
}
OMmodel::~OMmodel()
{
	glDeleteBuffers(1, &meshVBuffer);
	glDeleteBuffers(1, &meshNBuffer);
	glDeleteBuffers(1, &meshFNormal);
	glDeleteBuffers(1, &meshVColor);
}
bool OMmodel::OpenMeshReadFile(const char * filename)
{
	// request vertex normals, so the mesh reader can use normal information
	// if available
	mesh.request_vertex_normals();
	// request face normals,
	mesh.request_face_normals();
	//request vertex color
	mesh.request_vertex_colors();
	// assure we have vertex normals
	if (!mesh.has_vertex_normals())
	{
		std::cerr << "ERROR: Standard vertex property 'Normals' not available!\n";
		return 1;
	}

	OpenMesh::IO::Options opt;
	// read mesh from file
	if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
	{
		std::cerr << "Error: Cannot read mesh from " << filename << std::endl;
		return 1;
	}

	// If the file did not provide vertex normals, then calculate them
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// let the mesh update the normals
		mesh.update_normals();
	}

	// Get the face-vertex circulator of face _fh
	// get the vertex nomal
	int valence = 0;
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(*fv_it); vv_it.is_valid(); ++vv_it)
			{
				++valence;
			}
			if (valence <= 4)
			{
				meshVertexColorBuffer.push_back(OpenMesh::Vec3f(0.0f, 0.0f, 1.0f));
			}
			if (valence >= 5 && valence <= 7)
			{
				meshVertexColorBuffer.push_back(OpenMesh::Vec3f(0.0f, 1.0f, 0.0f));
			}
			if (valence >= 8)
			{
				meshVertexColorBuffer.push_back(OpenMesh::Vec3f(1.0f, 0.0f, 0.0f));
			}
			//cout<< valence <<endl;
			valence = 0;
			meshVertexBuffer.push_back(mesh.point(*fv_it));
			meshVertexNormalBuffer.push_back(mesh.normal(*fv_it));
			meshFaceNormalBuffer.push_back(mesh.normal(*f_it));
		}
	}
	// don't need the normals anymore? Remove them!
	mesh.release_vertex_normals();
	// dispose the face normals, as we don't need them anymore
	mesh.release_face_normals();
	//release color
	mesh.release_vertex_colors();


	mesh.request_vertex_status();
	meshVetexNum = mesh.n_vertices();
	mesh.request_face_status();
	meshFaceNum = mesh.n_faces();
	mesh.request_halfedge_status();
	meshHalfEdgeNum = mesh.n_halfedges();

	// iterate over all halfedges
	for (MyMesh::HalfedgeIter h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it)
	{
		if (!mesh.face_handle(*h_it).is_valid())
		{
			++meshBoundryEdgeNum;
		}
	}

	mesh.release_face_status();
	mesh.release_vertex_status();
	mesh.release_halfedge_status();
}
void OMmodel::RenderModel()
{
	glGenBuffers(1, &meshVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		meshVertexBuffer.size() * sizeof(OpenMesh::Vec3f),
		&(meshVertexBuffer[0]),
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	if (mSM == SmoothShading)
	{
		glGenBuffers(1, &meshNBuffer);
		glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
		glBufferData(GL_ARRAY_BUFFER, meshVertexNormalBuffer.size() * sizeof(OpenMesh::Vec3f), &meshVertexNormalBuffer[0], GL_STATIC_DRAW);

		// 2rd attribute buffer : normals
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, meshNBuffer);
		glVertexAttribPointer(
			1,                                // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			NULL                         // array buffer offset
			);
	}
	if (mSM == FlatShading)
	{
		glGenBuffers(1, &meshFNormal);
		glBindBuffer(GL_ARRAY_BUFFER, meshFNormal);
		glBufferData(GL_ARRAY_BUFFER, meshFaceNormalBuffer.size() * sizeof(OpenMesh::Vec3f), &meshFaceNormalBuffer[0], GL_STATIC_DRAW);

		// 2rd attribute buffer : normals
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, meshFNormal);
		glVertexAttribPointer(
			1,                                // attribute
			3,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			NULL                         // array buffer offset
			);
	}

	glDrawArrays(GL_TRIANGLES, 0, meshVertexBuffer.size());
}
void OMmodel::RenderModelWithColor()
{
	glGenBuffers(1, &meshVBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glBufferData(GL_ARRAY_BUFFER,
		meshVertexBuffer.size() * sizeof(OpenMesh::Vec3f),
		&(meshVertexBuffer[0]),
		GL_STATIC_DRAW);
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVBuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);


	glGenBuffers(1, &meshVColor);
	glBindBuffer(GL_ARRAY_BUFFER, meshVColor);
	glBufferData(GL_ARRAY_BUFFER, meshVertexColorBuffer.size() * sizeof(OpenMesh::Vec3f), &(meshVertexColorBuffer[0]), GL_STATIC_DRAW);
	// 2rd attribute buffer : vertices
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, meshVColor);
	glVertexAttribPointer(
		1,                                // attribute
		3,                                // size
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		0,                                // stride
		NULL                         // array buffer offset
		);

	glDrawArrays(GL_TRIANGLES, 0, meshVertexBuffer.size());
}