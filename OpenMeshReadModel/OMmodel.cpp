#include "OMmodel.h"

ofstream debug("debug.txt");
ofstream debug2("debug2.txt");
float cot(float d)
{
	return (1 / tan(d));
}
float Acos(OpenMesh::Vec3f vec1, OpenMesh::Vec3f vec2, OpenMesh::Vec3f vec3)
{
	float c = vec3.norm();
	float b = vec2.norm();
	float a = vec1.norm();
	float costheta = (a*a + b*b - c*c) / (2 * a * b);
	return acos(costheta);
}
OpenMesh::Vec3f RGBtoHSV(OpenMesh::Vec3f colorRGB){

	float Max = max(colorRGB[0], max(colorRGB[1], colorRGB[2]));
	float Min = min(colorRGB[0], min(colorRGB[1], colorRGB[2]));
	OpenMesh::Vec3f hsv;
	// h
	if (Max == Min){
		hsv[0] = 0;
	}
	else if (Max == colorRGB[0]){
		hsv[0] = (int)(60 * (colorRGB[1] - colorRGB[2]) / (Max - Min) + 360) % 360;
	}
	else if (Max == colorRGB[1]){
		hsv[0] = (60 * (colorRGB[2] - colorRGB[0]) / (Max - Min)) + 120;
	}
	else if (Max == colorRGB[2]){
		hsv[0] = (60 * (colorRGB[0] - colorRGB[1]) / (Max - Min)) + 240;
	}

	// s
	if (Max == 0){
		hsv[1] = 0;
	}
	else{
		hsv[1] = (255 * ((Max - Min) / Max));
	}

	// v
	hsv[2] = Max;

	return hsv;
}
OpenMesh::Vec3f HSVtoRGB(OpenMesh::Vec3f hsv){
	float f;
	int i;
	float p, q, t;
	OpenMesh::Vec3f rgb;

	i = (int)floor(hsv[0] / 60.0f) % 6;
	f = (float)(hsv[0] / 60.0f) - (float)floor(hsv[0] / 60.0f);
	p = (int)round(hsv[2] * (1.0f - (hsv[1] / 255.0f)));
	q = (int)round(hsv[2] * (1.0f - (hsv[1] / 255.0f) * f));
	t = (int)round(hsv[2] * (1.0f - (hsv[1] / 255.0f) * (1.0f - f)));

	switch (i){
	case 0: rgb[0] = hsv[2];       rgb[1] = t;            rgb[2] = p;         break;
	case 1: rgb[0] = q;            rgb[1] = hsv[2];       rgb[2] = p;         break;
	case 2: rgb[0] = p;            rgb[1] = hsv[2];       rgb[2] = t;         break;
	case 3: rgb[0] = p;            rgb[1] = q;            rgb[2] = hsv[2];    break;
	case 4: rgb[0] = t;            rgb[1] = p;            rgb[2] = hsv[2];    break;
	case 5: rgb[0] = hsv[2];       rgb[1] = p;            rgb[2] = q;         break;
	}

	return rgb;
}
OpenMesh::Vec3f interporlationColor(float max, float min, float now)
{
	/*
	debug << now <<"xxxx"<< max << min <<endl;

	OpenMesh::Vec3f green = OpenMesh::Vec3f(0.0f, 1.0f, 0.0f);
	OpenMesh::Vec3f red = OpenMesh::Vec3f(1.0f, 0.0f, 0.0f);
	OpenMesh::Vec3f blue = OpenMesh::Vec3f(0.0f, 0.0f, 1.0f);

	if (now >= min && now <= max / 3)
		return blue;
	else
		if (now >= max / 3 && now <= max*2 / 3)
			return green;
		else
			return red;
	*/
	float factor = (now - min) / (max - min);

	return HSVtoRGB(OpenMesh::Vec3f(factor*360, 1.0f, 1.0f));
}
OMmodel::OMmodel()
{
}
OMmodel::~OMmodel()
{
	if (debug.good())
		debug.close();
	if (debug2.good())
		debug2.close();
	glDeleteBuffers(1, &meshVBuffer);
	glDeleteBuffers(1, &meshNBuffer);
	glDeleteBuffers(1, &meshFNormal);
	glDeleteBuffers(1, &meshVColor);
}
void OMmodel::OpenMeshReadFile(const char * filename)
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
	}

	OpenMesh::IO::Options opt;
	// read mesh from file
	if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
	{
		std::cerr << "Error: Cannot read mesh from " << filename << std::endl;
	}

	// If the file did not provide vertex normals, then calculate them
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// let the mesh update the normals
		mesh.update_normals();
	}
	// this vertex property stores the computed  mean curvature
	OpenMesh::VPropHandleT<double> HCurvature;
	mesh.add_property(HCurvature);
	//store valence
	OpenMesh::VPropHandleT<int> valence;
	mesh.add_property(valence);
	//store gaussian curvature
	OpenMesh::VPropHandleT<double> GCurvature;
	mesh.add_property(GCurvature);
	//caculate curvature and valence
	vector<MyMesh::Point> oneRing;
	float maxCur = 0;
	float minCur = 0;
	float maxGCur = 0;
	float minGCur = 0;
    int val = 0;

	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		oneRing.clear();
		for (MyMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
		{
			++val;
			oneRing.push_back(mesh.point(*vv_it));
		}

		OpenMesh::Vec3f H = OpenMesh::Vec3f(0.0f, 0.0f, 0.0f);
		float A = 0;
		float Gcurvature = 0;
		float Hcurvature = 0;
		float theta = 0;
		float alpha = 0;
		float beta = 0;
		float arccosAlpha = 0;
		float arccosBeta = 0;
		float arccosTheta = 0;
		float cotAlpha = 0;
		float cotBeta = 0;
		for (int i = 0; i < val; ++i)
		{
			OpenMesh::Vec3f Vvtovi, Vvtovi1, Vvi_1tovi, Vvi_1tov, Vvi1tovi, Vvi1tov, Vvitovi1;
			MyMesh::Point PositionV = mesh.point(*v_it);
			Vvtovi = oneRing[i] - PositionV;
			OpenMesh::Vec3f nVvtovi = Vvtovi.normalized();
			if (i == 0)
			{
				Vvi_1tovi = oneRing[0] - oneRing[val - 1];
				Vvi_1tov = PositionV - oneRing[val - 1];
				Vvi_1tovi.normalize();
				Vvi_1tov.normalize();
			}
			else
			{
				Vvi_1tovi = oneRing[i] - oneRing[i - 1];
				Vvi_1tov = PositionV - oneRing[i - 1];
				Vvi_1tovi.normalize();
				Vvi_1tov.normalize();
			}
			if (i == val - 1)
			{
				Vvtovi1 = oneRing[0] - PositionV;
				Vvi1tovi = oneRing[i] - oneRing[0];
				Vvi1tov = PositionV - oneRing[0];
				Vvtovi1.normalize();
				Vvi1tovi.normalize();
				Vvi1tov.normalize();
			}
			else
			{
				Vvtovi1 = oneRing[i + 1] - PositionV;
				Vvi1tovi = oneRing[i] - oneRing[i + 1];
				Vvi1tov = PositionV - oneRing[i + 1];
				Vvtovi1.normalize();
				Vvi1tovi.normalize();
				Vvi1tov.normalize();
			}

			alpha = dot(Vvi_1tovi, Vvi_1tov);
			beta = dot(Vvi1tov, Vvi1tovi);
			arccosAlpha = acos(alpha);
			arccosBeta = acos(beta);
			cotAlpha = cot(arccosAlpha);
			cotBeta = cot(arccosBeta);
			A += ((cotAlpha + cotBeta) * Vvtovi.sqrnorm());
			H += ((cotAlpha + cotBeta) * Vvtovi);

			//theta += Acos(Vvtovi, Vvtovi1, Vvi1tovi);
			theta = dot(nVvtovi, Vvtovi1);
			//debug2 << nVvtovi.sqrnorm()<<" "<<acos(theta) << endl;
			arccosTheta += acos(theta);

		}
		//debug << arccosTheta << " " << H << " " << A << endl;;
		A = A / 8.0f;
		H = 2.0f * (H / A);
		Hcurvature = 0.5f * H.norm();
		Gcurvature = fabs((2.0f * M_PI) - arccosTheta)/ A;
		maxCur = max(Hcurvature, maxCur);
		minCur = min(Hcurvature, minCur);
		maxGCur = max(Gcurvature, maxGCur);
		minGCur = min(Gcurvature, minGCur);
		mesh.property(valence, *v_it) = val;
		mesh.property(HCurvature, *v_it) = Hcurvature;
		mesh.property(GCurvature, *v_it) = Gcurvature;

		debug << val << " " << Hcurvature << " " << Gcurvature << endl;

		val = 0;
	}
	std::cout << 2.0f * M_PI << endl;
	std::cout << maxCur << "   " << minCur <<endl;
	std::cout << maxGCur << "   " << minGCur << endl;

	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			int value = mesh.property(valence, *fv_it);
			if (value <= 4)
			{
				meshVertexColorBuffer.push_back(OpenMesh::Vec3f(0.0f, 0.0f, 1.0f));
			}
			if (value >= 5 && value <= 7)
			{
				meshVertexColorBuffer.push_back(OpenMesh::Vec3f(0.0f, 1.0f, 0.0f));
			}
			if (value >= 8)
			{
				meshVertexColorBuffer.push_back(OpenMesh::Vec3f(1.0f, 0.0f, 0.0f));
			}

			//cout << mesh.property(curvature, *fv_it) << endl;
			meshCurColorBuffer.push_back(interporlationColor(maxCur/20, minCur, mesh.property(HCurvature, *fv_it)));
			meshGCurColorBuffer.push_back(interporlationColor(maxGCur/20, 0, mesh.property(GCurvature, *fv_it)));
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
	if (mCM == Valence)
	    glBufferData(GL_ARRAY_BUFFER, meshVertexColorBuffer.size() * sizeof(OpenMesh::Vec3f), &(meshVertexColorBuffer[0]), GL_STATIC_DRAW);
	if (mCM == MeanCurvature)
	    glBufferData(GL_ARRAY_BUFFER, meshCurColorBuffer.size() * sizeof(OpenMesh::Vec3f), &(meshCurColorBuffer[0]), GL_STATIC_DRAW);
	if (mCM == GaussianCurvature)
		glBufferData(GL_ARRAY_BUFFER, meshGCurColorBuffer.size() * sizeof(OpenMesh::Vec3f), &(meshGCurColorBuffer[0]), GL_STATIC_DRAW);
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

