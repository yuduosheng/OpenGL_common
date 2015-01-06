#include "deformableobject.h"

DeformableObject::DeformableObject()
{
    d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
    d16 = (1.0f - nu) * d15;
    d17 = nu * d15;
    d18 = Y / 2 / (1.0f + nu);

	D = glm::vec3(d16, d17, d18); //Isotropic elasticity matrix D

	ReadModelFromFile("bunny.1");
	total_tetrahedra = tetrahedra.size();

	total_points = X.size();
	mass.resize(total_points);

	//copy positions to buffer 
	A_row.resize(total_points);
	K_row.resize(total_points);
	b.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);
	F0.resize(total_points);
	residual.resize(total_points);
	update.resize(total_points);
	prev_1.resize(total_points);

	//fill in V
	memset(&(V[0].x), 0, total_points*sizeof(glm::vec3));

	CalculateK();
	ClearStiffnessAssembly();
	RecalcMassMatrix();
	InitializePlastic();
}
DeformableObject::~DeformableObject()
{
	OnShutdown();
}

void DeformableObject::CalculateK() {

	for (size_t k = 0; k<tetrahedra.size(); k++) {

		glm::vec3 x0 = Xi[tetrahedra[k].indices[0]];
		glm::vec3 x1 = Xi[tetrahedra[k].indices[1]];
		glm::vec3 x2 = Xi[tetrahedra[k].indices[2]];
		glm::vec3 x3 = Xi[tetrahedra[k].indices[3]];

		//For this check page no.: 344-346 of Kenny Erleben's book Physics based Animation
		//Eq. 10.30(a-c)
		glm::vec3 e10 = x1 - x0;
		glm::vec3 e20 = x2 - x0;
		glm::vec3 e30 = x3 - x0;
		glm::vec3 e12 = x1 - x2;
		glm::vec3 e23 = x2 - x3;
		glm::vec3 e31 = x3 - x1;
		totalLength += glm::length(e10);
		totalLength += glm::length(e20);
		totalLength += glm::length(e30);
		totalLength += glm::length(e12);
		totalLength += glm::length(e23);
		totalLength += glm::length(e31);
		tetrahedra[k].e1 = e10;
		tetrahedra[k].e2 = e20;
		tetrahedra[k].e3 = e30;

		tetrahedra[k].volume = GetTetraVolume(e10, e20, e30);

		//Eq. 10.32
		glm::mat3 E = glm::mat3(e10.x, e10.y, e10.z,
			e20.x, e20.y, e20.z,
			e30.x, e30.y, e30.z);
		float detE = glm::determinant(E);
		float invDetE = 1.0f / detE;

		//Eq. 10.40 (a) & Eq. 10.42 (a)
		//Shape function derivatives wrt x,y,z
		// d/dx N0
		float invE10 = (e20.z*e30.y - e20.y*e30.z)*invDetE;
		float invE20 = (e10.y*e30.z - e10.z*e30.y)*invDetE;
		float invE30 = (e10.z*e20.y - e10.y*e20.z)*invDetE;
		float invE00 = -invE10 - invE20 - invE30;

		//Eq. 10.40 (b) & Eq. 10.42 (b)
		// d/dy N0
		float invE11 = (e20.x*e30.z - e20.z*e30.x)*invDetE;
		float invE21 = (e10.z*e30.x - e10.x*e30.z)*invDetE;
		float invE31 = (e10.x*e20.z - e10.z*e20.x)*invDetE;
		float invE01 = -invE11 - invE21 - invE31;

		//Eq. 10.40 (c) & Eq. 10.42 (c)
		// d/dz N0
		float invE12 = (e20.y*e30.x - e20.x*e30.y)*invDetE;
		float invE22 = (e10.x*e30.y - e10.y*e30.x)*invDetE;
		float invE32 = (e10.y*e20.x - e10.x*e20.y)*invDetE;
		float invE02 = -invE12 - invE22 - invE32;

		glm::vec3 B[4];
		//Eq. 10.43 
		//Bn ~ [bn cn dn]^T
		// bn = d/dx N0 = [ invE00 invE10 invE20 invE30 ]
		// cn = d/dy N0 = [ invE01 invE11 invE21 invE31 ]
		// dn = d/dz N0 = [ invE02 invE12 invE22 invE32 ]
		tetrahedra[k].B[0] = glm::vec3(invE00, invE01, invE02);
		tetrahedra[k].B[1] = glm::vec3(invE10, invE11, invE12);
		tetrahedra[k].B[2] = glm::vec3(invE20, invE21, invE22);
		tetrahedra[k].B[3] = glm::vec3(invE30, invE31, invE32);

		for (int i = 0; i<4; i++) {
			for (int j = 0; j<4; j++) {
				glm::mat3& Ke = tetrahedra[k].Ke[i][j];
				float d19 = tetrahedra[k].B[i].x;
				float d20 = tetrahedra[k].B[i].y;
				float d21 = tetrahedra[k].B[i].z;
				float d22 = tetrahedra[k].B[j].x;
				float d23 = tetrahedra[k].B[j].y;
				float d24 = tetrahedra[k].B[j].z;
				Ke[0][0] = d16 * d19 * d22 + d18 * (d20 * d23 + d21 * d24);
				Ke[0][1] = d17 * d19 * d23 + d18 * (d20 * d22);
				Ke[0][2] = d17 * d19 * d24 + d18 * (d21 * d22);

				Ke[1][0] = d17 * d20 * d22 + d18 * (d19 * d23);
				Ke[1][1] = d16 * d20 * d23 + d18 * (d19 * d22 + d21 * d24);
				Ke[1][2] = d17 * d20 * d24 + d18 * (d21 * d23);

				Ke[2][0] = d17 * d21 * d22 + d18 * (d19 * d24);
				Ke[2][1] = d17 * d21 * d23 + d18 * (d20 * d24);
				Ke[2][2] = d16 * d21 * d24 + d18 * (d20 * d23 + d19 * d22);

				Ke *= tetrahedra[k].volume;
			}
		}
	}
}
void DeformableObject::ReadModelFromFile(const char *filename)
{
	string fileName = filename;
	string suffix = ".node";
	string suffix2 = ".ele";
	fileName += suffix;
	ifstream nodeFile(fileName.c_str(), ios::in);
	if (!nodeFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}

	nodeFile >> total_points;
	X.resize(total_points);
	Xi.resize(total_points);
	IsFixed.resize(total_points);
	int noUse;
	nodeFile >> noUse >> noUse >> noUse;
	int ind = 0;
	glm::vec3 nodePosition;
	for (int i = 0; i < total_points; ++i)
	{
		nodeFile >> noUse;
		nodeFile >> nodePosition.x >> nodePosition.y >> nodePosition.z;
		//debug << nodePosition.x << " " << nodePosition.y << " " << nodePosition.z << endl;
		X[i] = nodePosition;
		Xi[i] = X[i];
		//Make the first few points fixed
		if (Xi[i].y < 0.05)
		{
			IsFixed[i] = true;
			ind = i;
		}
		else
		{
			IsFixed[i] = false;
		}
	}
	//move to floor
	float k = X[ind].y;
	for (size_t i = 0; i<total_points; i++) {
		X[i].y -= k;
	}

	nodeFile.close();

	fileName = filename;
	fileName += suffix2;

	ifstream eleFile(fileName.c_str(), ios::in);
	if (!eleFile.is_open())
	{
		cout << fileName << " open fail." << endl;
		return;
	}

	eleFile >> total_tetrahedra;
	eleFile >> noUse >> noUse;

	for (int i = 0; i < total_tetrahedra; ++i)
	{
		int p0, p1, p2, p3;
		eleFile >> noUse;
		eleFile >> p0 >> p1 >> p2 >> p3;
		//debug << p0 << " " << p1 << " " << p2 << " " << p3 <<endl;
		AddTetrahedron(p0, p1, p2, p3);
	}
	eleFile.close();
}


void DeformableObject::RecalcMassMatrix() {
	//This is a lumped mass matrix
	//Based on Eq. 10.106 and pseudocode in Fig. 10.9 on page 358
	for (size_t i = 0; i<total_points; i++) {
		if (IsFixed[i])
			mass[i] = std::numeric_limits<float>::max();
		else
			mass[i] = 1.0f / total_points;
	}

	for (int i = 0; i<total_tetrahedra; i++) {
		float m = (density*tetrahedra[i].volume)* 0.25f;
		mass[tetrahedra[i].indices[0]] += m;
		mass[tetrahedra[i].indices[1]] += m;
		mass[tetrahedra[i].indices[2]] += m;
		mass[tetrahedra[i].indices[3]] += m;
	}
}

void DeformableObject::InitializePlastic() {
	for (size_t i = 0; i<tetrahedra.size(); i++) {
		for (int j = 0; j<6; j++)
			tetrahedra[i].plastic[j] = 0;
	}
}

void DeformableObject::ClearStiffnessAssembly() {
	for (size_t k = 0; k<total_points; k++) {
		F0[k].x = 0.0f;
		F0[k].y = 0.0f;
		F0[k].z = 0.0f;

		for (matrix_iterator Kij = K_row[k].begin(); Kij != K_row[k].end(); ++Kij)
			Kij->second = glm::mat3(0);
	}
}


void DeformableObject::OnShutdown() {
	X.clear();
	Xi.clear();
	V.clear();
	mass.clear();
	F.clear();
	IsFixed.clear();
	tetrahedra.clear();
	K_row.clear();
	A_row.clear();
	F0.clear();
	b.clear();
	residual.clear();
	prev_1.clear();
	update.clear();
}

void DeformableObject::ComputeForces() {
	size_t i = 0;
	for (i = 0; i<total_points; i++) {
		F[i] = glm::vec3(0);

		//add gravity force only for non-fixed points
		F[i] += gravity*mass[i];
	}
}


glm::mat3 DeformableObject::ortho_normalize(glm::mat3 A) {
	glm::vec3 row0(A[0][0], A[0][1], A[0][2]);
	glm::vec3 row1(A[1][0], A[1][1], A[1][2]);
	glm::vec3 row2(A[2][0], A[2][1], A[2][2]);

	float L0 = glm::length(row0);
	if (L0)
		row0 /= L0;

	row1 -= row0 * glm::dot(row0, row1);
	float L1 = glm::length(row1);
	if (L1)
		row1 /= L1;

	row2 = glm::cross(row0, row1);

	return glm::mat3(row0,
		row1,
		row2);
}

void DeformableObject::UpdateOrientation() {
	for (int k = 0; k<total_tetrahedra; k++) {
		//Based on description on page 362-364 
		float div6V = 1.0f / tetrahedra[k].volume*6.0f;

		int i0 = tetrahedra[k].indices[0];
		int i1 = tetrahedra[k].indices[1];
		int i2 = tetrahedra[k].indices[2];
		int i3 = tetrahedra[k].indices[3];

		glm::vec3 p0 = X[i0];
		glm::vec3 p1 = X[i1];
		glm::vec3 p2 = X[i2];
		glm::vec3 p3 = X[i3];

		glm::vec3 e1 = p1 - p0;
		glm::vec3 e2 = p2 - p0;
		glm::vec3 e3 = p3 - p0;

		//Eq. 10.129 on page 363 & Eq. 10.131 page 364
		//n1,n2,n3 approximate Einv
		glm::vec3 n1 = glm::cross(e2, e3) * div6V;
		glm::vec3 n2 = glm::cross(e3, e1) * div6V;
		glm::vec3 n3 = glm::cross(e1, e2) * div6V;

		//Now get the rotation matrix from the initial undeformed (model/material coordinates)
		//We get the precomputed edge values
		e1 = tetrahedra[k].e1;
		e2 = tetrahedra[k].e2;
		e3 = tetrahedra[k].e3;

		//Based on Eq. 10.133		

		tetrahedra[k].Re[0][0] = e1.x * n1.x + e2.x * n2.x + e3.x * n3.x;
		tetrahedra[k].Re[0][1] = e1.x * n1.y + e2.x * n2.y + e3.x * n3.y;
		tetrahedra[k].Re[0][2] = e1.x * n1.z + e2.x * n2.z + e3.x * n3.z;

		tetrahedra[k].Re[1][0] = e1.y * n1.x + e2.y * n2.x + e3.y * n3.x;
		tetrahedra[k].Re[1][1] = e1.y * n1.y + e2.y * n2.y + e3.y * n3.y;
		tetrahedra[k].Re[1][2] = e1.y * n1.z + e2.y * n2.z + e3.y * n3.z;

		tetrahedra[k].Re[2][0] = e1.z * n1.x + e2.z * n2.x + e3.z * n3.x;
		tetrahedra[k].Re[2][1] = e1.z * n1.y + e2.z * n2.y + e3.z * n3.y;
		tetrahedra[k].Re[2][2] = e1.z * n1.z + e2.z * n2.z + e3.z * n3.z;

		tetrahedra[k].Re = ortho_normalize(tetrahedra[k].Re);

	}
}
void DeformableObject::ResetOrientation() {
	for (int k = 0; k<total_tetrahedra; k++) {
		tetrahedra[k].Re = I;
	}
}

void DeformableObject::StiffnessAssembly() {

	for (int k = 0; k<total_tetrahedra; k++) {
		glm::mat3 Re = tetrahedra[k].Re;
		glm::mat3 ReT = glm::transpose(Re);


		for (int i = 0; i < 4; ++i) {
			//Based on pseudocode given in Fig. 10.11 on page 361
			glm::vec3 f = glm::vec3(0.0f, 0.0f, 0.0f);
			for (int j = 0; j < 4; ++j) {
				glm::mat3 tmpKe = tetrahedra[k].Ke[i][j];
				glm::vec3 x0 = Xi[tetrahedra[k].indices[j]];
				glm::vec3 prod = glm::vec3(tmpKe[0][0] * x0.x + tmpKe[0][1] * x0.y + tmpKe[0][2] * x0.z, //tmpKe*x0;
					tmpKe[1][0] * x0.x + tmpKe[1][1] * x0.y + tmpKe[1][2] * x0.z,
					tmpKe[2][0] * x0.x + tmpKe[2][1] * x0.y + tmpKe[2][2] * x0.z);
				f += prod;
				if (j >= i) {
					//Based on pseudocode given in Fig. 10.12 on page 361
					glm::mat3 tmp = Re*tmpKe*ReT;
					int index = tetrahedra[k].indices[i];

					K_row[index][tetrahedra[k].indices[j]] += (tmp);

					if (j > i) {
						index = tetrahedra[k].indices[j];
						K_row[index][tetrahedra[k].indices[i]] += (glm::transpose(tmp));
					}
				}

			}
			int idx = tetrahedra[k].indices[i];
			F0[idx] -= Re*f;
		}
	}
}

void DeformableObject::AddPlasticityForce(float dt) {
	for (int k = 0; k<total_tetrahedra; k++) {
		float e_total[6];
		float e_elastic[6];
		for (int i = 0; i<6; ++i)
			e_elastic[i] = e_total[i] = 0;

		//--- Compute total strain: e_total  = Be (Re^{-1} x - x0)
		for (unsigned int j = 0; j<4; ++j) {

			glm::vec3 x_j = X[tetrahedra[k].indices[j]];
			glm::vec3 x0_j = Xi[tetrahedra[k].indices[j]];
			glm::mat3 ReT = glm::transpose(tetrahedra[k].Re);
			glm::vec3 prod = glm::vec3(ReT[0][0] * x_j.x + ReT[0][1] * x_j.y + ReT[0][2] * x_j.z, //tmpKe*x0;
				ReT[1][0] * x_j.x + ReT[1][1] * x_j.y + ReT[1][2] * x_j.z,
				ReT[2][0] * x_j.x + ReT[2][1] * x_j.y + ReT[2][2] * x_j.z);

			glm::vec3 tmp = prod - x0_j;

			//B contains Jacobian of shape funcs. B=SN
			float bj = tetrahedra[k].B[j].x;
			float cj = tetrahedra[k].B[j].y;
			float dj = tetrahedra[k].B[j].z;

			e_total[0] += bj*tmp.x;
			e_total[1] += cj*tmp.y;
			e_total[2] += dj*tmp.z;
			e_total[3] += cj*tmp.x + bj*tmp.y;
			e_total[4] += dj*tmp.x + bj*tmp.z;
			e_total[5] += dj*tmp.y + cj*tmp.z;
		}

		//--- Compute elastic strain
		for (int i = 0; i<6; ++i)
			e_elastic[i] = e_total[i] - tetrahedra[k].plastic[i];

		//--- if elastic strain exceeds c_yield then it is added to plastic strain by c_creep
		float norm_elastic = 0;
		for (int i = 0; i<6; ++i)
			norm_elastic += e_elastic[i] * e_elastic[i];
		norm_elastic = sqrt(norm_elastic);
		if (norm_elastic > yield) {
			float amount = dt*glm::min(creep, (1.0f / dt));  //--- make sure creep do not exceed 1/dt
			for (int i = 0; i<6; ++i)
				tetrahedra[k].plastic[i] += amount*e_elastic[i];
		}

		//--- if plastic strain exceeds c_max then it is clamped to maximum magnitude
		float norm_plastic = 0;
		for (int i = 0; i<6; ++i)
			norm_plastic += tetrahedra[k].plastic[i] * tetrahedra[k].plastic[i];
		norm_plastic = sqrt(norm_plastic);

		if (norm_plastic > m_max) {
			float scale = m_max / norm_plastic;
			for (int i = 0; i<6; ++i)
				tetrahedra[k].plastic[i] *= scale;
		}

		for (size_t n = 0; n<4; ++n) {
			float* e_plastic = tetrahedra[k].plastic;
			//bn, cn and dn are the shape function derivative wrt. x,y and z axis
			//These were calculated in CalculateK function

			//Eq. 10.140(a) & (b) on page 365
			float bn = tetrahedra[k].B[n].x;
			float cn = tetrahedra[k].B[n].y;
			float dn = tetrahedra[k].B[n].z;
			float D0 = D.x;
			float D1 = D.y;
			float D2 = D.z;
			glm::vec3 f = glm::vec3(0);

			float  bnD0 = bn*D0;
			float  bnD1 = bn*D1;
			float  bnD2 = bn*D2;
			float  cnD0 = cn*D0;
			float  cnD1 = cn*D1;
			float  cnD2 = cn*D2;
			float  dnD0 = dn*D0;
			float  dnD1 = dn*D1;
			float  dnD2 = dn*D2;

			//Eq. 10.141 on page 365
			f.x = bnD0*e_plastic[0] + bnD1*e_plastic[1] + bnD1*e_plastic[2] + cnD2*e_plastic[3] + dnD2*e_plastic[4];
			f.y = cnD1*e_plastic[0] + cnD0*e_plastic[1] + cnD1*e_plastic[2] + bnD2*e_plastic[3] + +dnD2*e_plastic[5];
			f.z = dnD1*e_plastic[0] + dnD1*e_plastic[1] + dnD0*e_plastic[2] + bnD2*e_plastic[4] + cnD2*e_plastic[5];

			f *= tetrahedra[k].volume;
			int idx = tetrahedra[k].indices[n];
			F[idx] += tetrahedra[k].Re*f;
		}
	}
}
void DeformableObject::DynamicsAssembly(float dt) {
	float dt2 = dt*dt;

	for (size_t k = 0; k<total_points; k++) {

		float m_i = mass[k];
		b[k] = glm::vec3(0.0, 0.0, 0.0);
		//printf("Idx: %3d\n",k);

		matrix_map tmp = K_row[k];
		matrix_iterator Kbegin = tmp.begin();
		matrix_iterator Kend = tmp.end();
		for (matrix_iterator K = Kbegin; K != Kend; ++K)
		{
			unsigned int j = K->first;
			glm::mat3 K_ij = K->second;
			glm::vec3 x_j = X[j];
			glm::mat3& A_ij = A_row[k][j];

			A_ij = K_ij * (dt2);
			glm::vec3 prod = glm::vec3(K_ij[0][0] * x_j.x + K_ij[0][1] * x_j.y + K_ij[0][2] * x_j.z,
				K_ij[1][0] * x_j.x + K_ij[1][1] * x_j.y + K_ij[1][2] * x_j.z,
				K_ij[2][0] * x_j.x + K_ij[2][1] * x_j.y + K_ij[2][2] * x_j.z);

			b[k] -= prod;//K_ij * x_j;

			if (k == j)
			{
				float c_i = mass_damping*m_i;
				float tmp = m_i + dt*c_i;
				A_ij[0][0] += tmp;
				A_ij[1][1] += tmp;
				A_ij[2][2] += tmp;
			}
		}

		b[k] -= F0[k];
		b[k] += F[k];
		b[k] *= dt;
		b[k] += V[k] * m_i;
	}
}


void DeformableObject::ConjugateGradientSolver(float dt) {


	for (size_t k = 0; k<total_points; k++) {
		if (IsFixed[k])
			continue;
		residual[k] = b[k];

		matrix_iterator Abegin = A_row[k].begin();
		matrix_iterator Aend = A_row[k].end();
		for (matrix_iterator A = Abegin; A != Aend; ++A)
		{
			unsigned int j = A->first;
			glm::mat3& A_ij = A->second;
			float v_jx = V[j].x;
			float v_jy = V[j].y;
			float v_jz = V[j].z;
			glm::vec3 prod = glm::vec3(A_ij[0][0] * v_jx + A_ij[0][1] * v_jy + A_ij[0][2] * v_jz, //A_ij * prev_1[j]
				A_ij[1][0] * v_jx + A_ij[1][1] * v_jy + A_ij[1][2] * v_jz,
				A_ij[2][0] * v_jx + A_ij[2][1] * v_jy + A_ij[2][2] * v_jz);
			residual[k] -= prod;//  A_ij * v_j;

		}
		prev_1[k] = residual[k];
	}

	for (int i = 0; i<i_max; i++) {
		float d = 0;
		float d2 = 0;

		for (size_t k = 0; k<total_points; k++) {

			if (IsFixed[k])
				continue;

			update[k].x = 0; update[k].y = 0; update[k].z = 0;

			matrix_iterator Abegin = A_row[k].begin();
			matrix_iterator Aend = A_row[k].end();
			for (matrix_iterator A = Abegin; A != Aend; ++A) {
				unsigned int j = A->first;
				glm::mat3& A_ij = A->second;
				float prevx = prev_1[j].x;
				float prevy = prev_1[j].y;
				float prevz = prev_1[j].z;
				glm::vec3 prod = glm::vec3(A_ij[0][0] * prevx + A_ij[0][1] * prevy + A_ij[0][2] * prevz, //A_ij * prev_1[j]
					A_ij[1][0] * prevx + A_ij[1][1] * prevy + A_ij[1][2] * prevz,
					A_ij[2][0] * prevx + A_ij[2][1] * prevy + A_ij[2][2] * prevz);
				update[k] += prod;//A_ij*prev_1[j];

			}
			d += glm::dot(residual[k], residual[k]);
			d2 += glm::dot(prev_1[k], update[k]);
		}

		if (fabs(d2)<tiny)
			d2 = tiny;

		float d3 = d / d2;
		float d1 = 0;


		for (size_t k = 0; k<total_points; k++) {
			if (IsFixed[k])
				continue;

			V[k] += prev_1[k] * d3;
			residual[k] -= update[k] * d3;
			d1 += glm::dot(residual[k], residual[k]);
		}

		if (i >= i_max && d1 < tolerence)
			break;

		if (fabs(d)<tiny)
			d = tiny;

		float d4 = d1 / d;

		for (size_t k = 0; k<total_points; k++) {
			if (IsFixed[k])
				continue;
			prev_1[k] = residual[k] + prev_1[k] * d4;
		}
	}
}

void DeformableObject::UpdatePosition(float dt) {
	for (size_t k = 0; k<total_points; k++) {
		if (IsFixed[k])
			continue;
		X[k] += float(dt)*V[k];
	}
}

void DeformableObject::GroundCollision()
{
	for (size_t i = 0; i<total_points; i++) {
		if (X[i].y<0) //collision with ground
			X[i].y = 0;
	}
}

void DeformableObject::StepPhysics(float dt) {

	ComputeForces();

	ClearStiffnessAssembly();

	if (bUseStiffnessWarping)
		UpdateOrientation();
	else
		ResetOrientation();

	StiffnessAssembly();

	AddPlasticityForce(dt);

	DynamicsAssembly(dt);

	ConjugateGradientSolver(dt);
	//EigenSolve();

	UpdatePosition(dt);

	GroundCollision();
}

void DeformableObject::renderModel()
{
	glColor3f(0.75, 0.75, 0.75);
	glBegin(GL_LINES);

	for (int i = 0; i<total_tetrahedra; i++) {
		int i0 = tetrahedra[i].indices[0];
		int i1 = tetrahedra[i].indices[1];
		int i2 = tetrahedra[i].indices[2];
		int i3 = tetrahedra[i].indices[3];
		glm::vec3 p1 = X[i0];
		glm::vec3 p2 = X[i1];
		glm::vec3 p3 = X[i2];
		glm::vec3 p4 = X[i3];

		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();


	//draw points	
	glBegin(GL_POINTS);
	for (int i = 0; i<total_points; i++) {
		glm::vec3 p = X[i];
		int is = (i == selected_index);
		glColor3f((float)!is, (float)is, (float)is);
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
}