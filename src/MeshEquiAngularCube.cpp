//Author: Xavier Corbillon
//IMT Atlantique
#include "MeshEquiAngularCube.hpp"
#include <iostream>
//standard includes
#include <cmath>

using namespace IMT;
constexpr float PI = 3.141592653589793238462643383279502884L;
constexpr float U_BORDER = 2.25e-4;
constexpr float V_BORDER = 3.75e-4;

MeshEquiAngularCube::MeshEquiAngularCube(GLfloat scale, size_t numTriangles): Mesh()
{
  // Figure out how many quads we have per edge.  There
  // is a minimum of 1.
  size_t numQuads = numTriangles / 2;
  size_t numQuadsPerFace = numQuads / 6;
  size_t numQuadsPerEdge = static_cast<size_t> (
    std::sqrt(numQuadsPerFace));
  if (numQuadsPerEdge < 1) { numQuadsPerEdge = 1; }

  // Construct a white square with the specified number of
  // quads as the +Z face of the cube.  We'll copy this
  // and we'll
  // adjust the coordinates by rotation to match each face.
  std::vector<GLfloat> faceBufferData;
  std::vector<GLfloat> tmpUVBufferData;
  for (size_t i = 0; i < numQuadsPerEdge; i++) {
    for (size_t j = 0; j < numQuadsPerEdge; j++) {

      // Send the two triangles that make up this quad, where the
      // quad covers the appropriate fraction of the face from
      // -scale to scale in X and Y.
      GLfloat Z = scale;
      GLfloat minX = -scale + i * (2 * scale) / numQuadsPerEdge;
      GLfloat maxX = -scale + (i+1) * (2 * scale) / numQuadsPerEdge;
      GLfloat minY = -scale + j * (2 * scale) / numQuadsPerEdge;
      GLfloat maxY = -scale + (j + 1) * (2 * scale) / numQuadsPerEdge;
      {
        faceBufferData.push_back(minX);
        faceBufferData.push_back(maxY);
        faceBufferData.push_back(Z);

        faceBufferData.push_back(minX);
        faceBufferData.push_back(minY);
        faceBufferData.push_back(Z);

        faceBufferData.push_back(maxX);
        faceBufferData.push_back(minY);
        faceBufferData.push_back(Z);

        faceBufferData.push_back(maxX);
        faceBufferData.push_back(maxY);
        faceBufferData.push_back(Z);

        faceBufferData.push_back(minX);
        faceBufferData.push_back(maxY);
        faceBufferData.push_back(Z);

        faceBufferData.push_back(maxX);
        faceBufferData.push_back(minY);
        faceBufferData.push_back(Z);
      }
      // auto tmptmpUVBufferData = GetUVs(float(i), float(j), float(numQuadsPerEdge));
      // tmpUVBufferData.insert(tmpUVBufferData.end(),
      //   tmptmpUVBufferData.begin(), tmptmpUVBufferData.end());
    }
  }
  // Make a copy of the vertices for each face, then modulate
  // the color by the face color and rotate the coordinates to
  // put them on the correct cube face.

  std::vector<GLfloat> tmpVertexBufferData;
  

  // right face
  // +Z is blue and is in the same location as the original
  // faces.
  {
    // X = X, Y = Y, Z = Z
    auto faceId = 4;
    std::array<GLfloat, 3> scales = { 1.0f, 1.0f, 1.0f };
    std::array<size_t, 3> indices = { 0, 1, 2 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceId, faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
    // out.push_back((float)1/3);
    // out.push_back(0.0);
    // out.push_back((float)1/3);
    // out.push_back(0.5);
    // out.push_back(0.0);
    // out.push_back(0.5);
    // out.push_back(0.0);
    // out.push_back(0.0);
    // out.push_back((float)1/3);
    // out.push_back(0.0);
    // out.push_back(0.0);
    // out.push_back(0.5);
    
  }

  // left face
  // -Z is cyan and is in the opposite size from the
  // original face (mirror all 3).
  {
    // X = -X, Y = -Y, Z = -Z
    auto faceId = 5;
    std::array<GLfloat, 3> scales = { -1.0f, -1.0f, -1.0f };
    std::array<size_t, 3> indices = { 0, 1, 2 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceId, faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
    // out.push_back(1.0);
    // out.push_back(0.5);
    // out.push_back(1.0);
    // out.push_back(0.0);
    // out.push_back((float)2/3);
    // out.push_back(0.0);
    // out.push_back((float)2/3);
    // out.push_back(0.5);
    // out.push_back(1.0);
    // out.push_back(0.5);
    // out.push_back((float)2/3);
    // out.push_back(0.0);
   
  }

  // front face
  // +X is red and is rotated -90 degrees from the original
  // around Y.
  {
    // X = Z, Y = Y, Z = -X
    auto faceId = 0;
    std::array<GLfloat, 3> scales = { 1.0f, 1.0f, -1.0f };
    std::array<size_t, 3> indices = { 2, 1, 0 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceId, faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
    // out.push_back((float)2/3);
    // out.push_back(1.0);
    // out.push_back((float)1/3);
    // out.push_back(1.0);
    // out.push_back((float)1/3);
    // out.push_back(0.5);
    // out.push_back((float)2/3);
    // out.push_back(0.5);
    // out.push_back((float)2/3);
    // out.push_back(1.0);
    // out.push_back((float)1/3);
    // out.push_back(0.5);
    
  }

  // back face
  // -X is magenta and is rotated 90 degrees from the original
  // around Y.
  {
    // X = -Z, Y = Y, Z = X
    auto faceId = 1;
    std::array<GLfloat, 3> scales = { -1.0f, 1.0f, 1.0f };
    std::array<size_t, 3> indices = { 2, 1, 0 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceId, faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
    // out.push_back((float)2/3);
    // out.push_back(0.0);
    // out.push_back((float)2/3);
    // out.push_back(0.5);
    // out.push_back((float)1/3);
    // out.push_back(0.5);
    // out.push_back((float)1/3);
    // out.push_back(0.0);
    // out.push_back((float)2/3);
    // out.push_back(0.0);
    // out.push_back((float)1/3);
    // out.push_back(0.5);
    
  }

  // top face
  // +Y is green and is rotated -90 degrees from the original
  // around X.
  {
    // X = X, Y = Z, Z = -Y
    auto faceId = 2;
    std::array<GLfloat, 3> scales = { 1.0f, 1.0f, -1.0f };
    std::array<size_t, 3> indices = { 0, 2, 1 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceId, faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
    // out.push_back(1.0);
    // out.push_back(0.5);
    // out.push_back(1.0);
    // out.push_back(1.0);
    // out.push_back((float)2/3);
    // out.push_back(1.0);
    // out.push_back((float)2/3);
    // out.push_back(0.5);
    // out.push_back(1.0);
    // out.push_back(0.5);
    // out.push_back((float)2/3);
    // out.push_back(1.0);
    
  }

  // bottom face
  // -Y is yellow and is rotated 90 degrees from the original
  // around X.
  {

    // X = X, Y = -Z, Z = Y
    auto faceId = 3;
    std::array<GLfloat, 3> scales = { 1.0f, -1.0f, 1.0f };
    std::array<size_t, 3> indices = { 0, 2, 1 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceId, faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
    // out.push_back(0.0);
    // out.push_back(1.0);
    // out.push_back(0.0);
    // out.push_back(0.5);
    // out.push_back((float)1/3);
    // out.push_back(0.5);
    // out.push_back((float)1/3);
    // out.push_back(1.0);
    // out.push_back(0.0);
    // out.push_back(1.0);
    // out.push_back((float)1/3);
    // out.push_back(0.5);

  }

  AppendVertexBufferData(GetVertexs(tmpVertexBufferData));
  AppendUvBufferData(VertexToUVs(tmpVertexBufferData));
}


std::vector<GLfloat> MeshEquiAngularCube::VertexRotate(
    GLfloat const &faceId,
    std::vector<GLfloat> const &inVec,
    std::array<size_t, 3> const &indices,
    std::array<GLfloat, 3> const &scales)
{
  std::vector<GLfloat> out;
  size_t elements = inVec.size() / 3;
  if (elements * 3 != inVec.size()) {
    // We don't have an even multiple of 3 elements, so bail.
    return out;
  }
  // add faceId
  out.resize(inVec.size()+elements);
  for (size_t i = 0; i < elements; i++) {
    out[4 * i] = faceId;
    for (size_t p = 0; p < 3; p++) {
      out[4 * i + p + 1] = inVec[3*i + indices[p]] * scales[p];
    }
  }
  return out;
}

std::vector<GLfloat> MeshEquiAngularCube::GetVertexs( std::vector<GLfloat> const& inputVertexs)
{
  std::vector<GLfloat> out;
  for (size_t i = 0; i < inputVertexs.size(); i += 4)
  {
	auto& faceID = inputVertexs[i];
    auto& x = inputVertexs[i+1];
    auto& y = inputVertexs[i+2];
    auto& z = inputVertexs[i+3];
    out.push_back(x);
    out.push_back(y);
    out.push_back(z);	
  }
  return out;
}


std::vector<GLfloat> MeshEquiAngularCube::VertexToUVs( std::vector<GLfloat> const& inputVertexs)
{
  std::vector<GLfloat> out;
  //Generate CubeMap UV map
  for(size_t i = 0; i < inputVertexs.size(); i += 4){
    
    auto& faceId = inputVertexs[i];
    auto& x = inputVertexs[i+1];
    auto& y = inputVertexs[i+2];
    auto& z = inputVertexs[i+3];
    auto u = 0.0f, v = 0.0f;
	
	if (faceId == 1)
	{
		// back face
		u = y / abs(x);
		v = -z / abs(x);
	}
	else if (faceId == 0)
	{
		// front face
		u = -z / abs(x);
		v = -y / abs(x);
	}
	else if (faceId == 5)
	{
		// left face
		u = -x / abs(z);
		v = -y / abs(z);
	}
	else if (faceId == 4)
	{
		// right face
		u = x / abs(z);
		v = -y / abs(z);
	}
	else if (faceId == 2)
	{
		// top face
		u = x / abs(y);
		v = z / abs(y);
	}
	else
	{
		// bottom face
		u = x / abs(y);
		v = -z / abs(y);
	}

	u = 4.0 / PI * std::atan(u);
	v = 4.0 / PI * std::atan(v);

	u = (u + 1.0) / 2.0;
	v = (v + 1.0) / 2.0;
	if (faceId == 1)
	{
		// back face
		u = u / 3.0 + 1.0 / 3.0;
		v = v / 2.0 + 0.5;
		v = (v <= 0.5 + V_BORDER ? (0.5 + V_BORDER) : v);
	}
	else if (faceId == 0)
	{
		// front face
		u = -u / 3.0 + 2.0 / 3.0;
		v = v / 2.0;
		v = (v >= 0.5 - V_BORDER ? (0.5 - V_BORDER) : v);
	}
	else if (faceId == 5)
	{
		// left face
		u = -u / 3.0 + 1.0 / 3.0;
		v = v / 2.0;
		v = (v >= 0.5 - V_BORDER ? (0.5 - V_BORDER) : v);
	}
	else if (faceId == 4)
	{
		// right face
		u = -u / 3.0 + 1.0;
		v = v / 2.0;
		v = (v >= 0.5 - V_BORDER ? (0.5 - V_BORDER) : v);
	}
	else if (faceId == 2)
	{
		// top face
		u = u / 3.0 + 2.0 / 3.0;
		v = -v / 2.0 + 1.0;
		v = (v <= 0.5 + V_BORDER ? (0.5 + V_BORDER) : v);
	}
	else
	{
		// bottom face
		u = -u / 3.0 + 1.0 / 3.0;
		v = v / 2.0 + 0.5;
		v = (v <= 0.5 + V_BORDER ? (0.5 + V_BORDER) : v);
	}

	u = (u >= 1 - U_BORDER ? (1 - U_BORDER) : u);
	u = (u <= 0 + U_BORDER ? (0 + U_BORDER) : u);
	v = (v >= 1 - V_BORDER ? (1 - V_BORDER) : v);
	v = (v <= 0 + V_BORDER ? (0 + V_BORDER) : v);
	
    out.push_back(u);
    out.push_back(v);
    
  }
  return out;
}

void MeshEquiAngularCube::InitImpl(void)
{
  // UV
  glBindBuffer(GL_ARRAY_BUFFER, GetUvBufferId());
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);

  // VBO
  glBindBuffer(GL_ARRAY_BUFFER, GetVertexBufferId());
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
}
