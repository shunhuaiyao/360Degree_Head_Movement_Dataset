//Author: Xavier Corbillon
//IMT Atlantique
#include "MeshCubeEquiUV.hpp"
#include <iostream>

//standard includes
#include <cmath>

using namespace IMT;
constexpr float PI = 3.141592653589793238462643383279502884L;

MeshCubeEquiUV::MeshCubeEquiUV(GLfloat scale, size_t numTriangles): Mesh()
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
      { //regular square
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
    std::array<GLfloat, 3> scales = { 1.0f, 1.0f, 1.0f };
    std::array<size_t, 3> indices = { 0, 1, 2 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
  }

  // left face
  // -Z is cyan and is in the opposite size from the
  // original face (mirror all 3).
  {
    // X = -X, Y = -Y, Z = -Z
    std::array<GLfloat, 3> scales = { -1.0f, -1.0f, -1.0f };
    std::array<size_t, 3> indices = { 0, 1, 2 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
  }

  // front face
  // +X is red and is rotated -90 degrees from the original
  // around Y.
  {
    // X = Z, Y = Y, Z = -X
    std::array<GLfloat, 3> scales = { 1.0f, 1.0f, -1.0f };
    std::array<size_t, 3> indices = { 2, 1, 0 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
  }

  // back face
  // -X is magenta and is rotated 90 degrees from the original
  // around Y.
  {
    // X = -Z, Y = Y, Z = X
    std::array<GLfloat, 3> scales = { -1.0f, 1.0f, 1.0f };
    std::array<size_t, 3> indices = { 2, 1, 0 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
  }

  // top face
  // +Y is green and is rotated -90 degrees from the original
  // around X.
  {
    // X = X, Y = Z, Z = -Y
    std::array<GLfloat, 3> scales = { 1.0f, 1.0f, -1.0f };
    std::array<size_t, 3> indices = { 0, 2, 1 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
  }

  // bottom face
  // -Y is yellow and is rotated 90 degrees from the original
  // around X.
  {

    // X = X, Y = -Z, Z = Y
    std::array<GLfloat, 3> scales = { 1.0f, -1.0f, 1.0f };
    std::array<size_t, 3> indices = { 0, 2, 1 };
    std::vector<GLfloat> myFaceBufferData =
      VertexRotate(faceBufferData, indices, scales);

    tmpVertexBufferData.insert(tmpVertexBufferData.end(),
      myFaceBufferData.begin(), myFaceBufferData.end());
  }
  AppendVertexBufferData(tmpVertexBufferData);
  AppendUvBufferData(VertexToUVs(tmpVertexBufferData));

  // std::vector<GLfloat> tmpVertexBufferData;
  // VertexData* vertices;
  // GLuint* indices;

  // int slices = 10;
  // int layercuts = 10;
  // float x, y, z, u, v, d; //d=diameter of circle at height y
  // float theta; //theta used for angle about y axis
  // float phi = 0.0f; //phi used for angle to determine y and d.
  // float thetad = 2 * PI / slices; //theta delta
  // float phid = PI / (layercuts + 1); //phi delta
  // int NumVertices = (slices * layercuts) + 2;
  // int NumIndices = ((slices + 2) * 2) + ((layercuts - 1) * (slices + 1) * 2);
  // vertices = new VertexData[NumVertices];
  // indices = new GLuint[NumIndices];
  // unsigned int indexPerLayer = (slices * 2) + 2; //indices for quad layer
  // unsigned int indexPerPointLayer = slices + 2; //indices for triangle layer (at poles)
  // unsigned int layerStartVertex; //beginning vert of current layercut
  // unsigned int layerStartIndex; //beginning index of current layercut

  // vertices[0] = 

  // for (unsigned int i = 0; i < layercuts; i++) {
  //   y = std::cos(phi);
  //   d = std::sin(phi);
  //   theta = 0.0f;
  //   for (unsigned int j = 0; j < slices; j++) {
  //     x = std::sin(theta) * d;
  //     z = std::cos(theta) * d;
  //     //u = 0.5f + (std::atan2(z, x) / (2 * PI));
  //     //v = 0.5f - (std::asin(y) / PI);
  //     //vertices[(i * slices) + j + 1] = { { x / 2, y / 2, z / 2 },{ u, v },{ 1.0f, 1.0f, 1.0f },{ x, y, z } };
  //     tmpVertexBufferData.insert(tmpVertexBufferData.end(), x);
  //     tmpVertexBufferData.insert(tmpVertexBufferData.end(), y);
  //     tmpVertexBufferData.insert(tmpVertexBufferData.end(), z);
  //     //std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
  //     theta += thetad;
  //   }
  //   phi += phid;
  // }

  // tmpVertex[0] = 0.0f;
  // tmpVertex[1] = -1.0f;
  // tmpVertex[2] = 0.0f;
  // tmpVertexBufferData.insert(tmpVertexBufferData.end(), tmpVertex, tmpVertex+3);

  // std::cout << "tmpVertexBufferData size: " << tmpVertexBufferData.size() << std::endl;
  // AppendVertexBufferData(tmpVertexBufferData);
  // AppendUvBufferData(VertexToUVs(tmpVertexBufferData));

  // std::vector<GLfloat> vertices;
  // int lats = 1000;
  // int longs = 1000;
  // for(int i = 0; i <= lats; i++) {
  //   double lat0 = PI * (-0.5 + (double) (i - 1) / lats);
  //   double z0  = sin(lat0);
  //   double zr0 =  cos(lat0);

  //   double lat1 = PI * (-0.5 + (double) i / lats);
  //   double z1 = sin(lat1);
  //   double zr1 = cos(lat1);

  //   for(int j = 0; j <= longs; j++) {
  //     double lng = 2 * PI * (double) (j - 1) / longs;
  //     double x = cos(lng);
  //     double y = sin(lng);
  //     vertices.push_back(x * zr0);
  //     vertices.push_back(y * zr0);
  //     vertices.push_back(z0);
         
  //     vertices.push_back(x * zr1);
  //     vertices.push_back(y * zr1);
  //     vertices.push_back(z1);
  //   }
  // }
  // AppendVertexBufferData(vertices);
  // AppendUvBufferData(VertexToUVs(vertices));
}


std::vector<GLfloat> MeshCubeEquiUV::VertexRotate(
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
  out.resize(inVec.size());
  for (size_t i = 0; i < elements; i++) {
    for (size_t p = 0; p < 3; p++) {
      out[3 * i + p] = inVec[3*i + indices[p]] * scales[p];
    }
  }
  return out;
}

std::vector<GLfloat> MeshCubeEquiUV::VertexToUVs( std::vector<GLfloat> const& inputVertexs)
{
  std::vector<GLfloat> out;
  //Generate Equirectangular UV map
  for(size_t i = 0; i < inputVertexs.size(); i += 3)
  {
    auto& x = inputVertexs[i];
    auto& y = inputVertexs[i+1];
    auto& z = inputVertexs[i+2];
    auto rho = std::sqrt(x*x+y*y+z*z);
    auto theta = std::atan2(z, x) + PI;
    auto phi = std::acos(y / rho);
    auto u = (theta/(2.f*PI));
    auto v = (phi/PI);
    out.push_back(u);
    out.push_back(v);
  }
  //Check if two vertix are at opposite side of the texture
  for(size_t i = 0; i < out.size(); i += 6)
  {
      if (std::abs(out[i] - out[i+2]) > 0.7)
      {
        if (out[i] < 0.3)
        {
          out[i] += 1.f;
        }
        else
        {
          out[i+2] += 1.f;
        }
      }
      if (std::abs(out[i] - out[i+4]) > 0.7)
      {
        if (out[i] < 0.3)
        {
          out[i] += 1.f;
        }
        else
        {
          out[i+4] += 1.f;
        }
      }
      if (std::abs(out[i+2] - out[i+4]) > 0.7)
      {
        if (out[i+2] < 0.3)
        {
          out[i+2] += 1.f;
        }
        else
        {
          out[i+4] += 1.f;
        }
      }
  }

  // for(size_t i = 0; i < inputVertexs.size(); i += 3){
  //   auto& x = inputVertexs[i];
  //   auto& y = inputVertexs[i+1];
  //   auto& z = inputVertexs[i+2];

  //   auto rho = std::sqrt(x*x + y*y + z*z);
  //   auto theta = PI - std::atan2(z, x);
  //   auto phi = std::acos(y/rho);
    
  //   auto u = -(theta / (2*PI));
  //   auto v = (phi / PI);
  //   out.push_back(u);
  //   out.push_back(v);
  // }
  
  // //Check if two vertix are at opposite side of the texture
  // for(size_t i = 0; i < out.size(); i += 6)
  // {
  //     if (std::abs(out[i] - out[i+2]) > 0.7)
  //     {
  //       if (out[i] < 0.3)
  //       {
  //         out[i] += 1.f;
  //       }
  //       else
  //       {
  //         out[i+2] += 1.f;
  //       }
  //     }
  //     if (std::abs(out[i] - out[i+4]) > 0.7)
  //     {
  //       if (out[i] < 0.3)
  //       {
  //         out[i] += 1.f;
  //       }
  //       else
  //       {
  //         out[i+4] += 1.f;
  //       }
  //     }
  //     if (std::abs(out[i+2] - out[i+4]) > 0.7)
  //     {
  //       if (out[i+2] < 0.3)
  //       {
  //         out[i+2] += 1.f;
  //       }
  //       else
  //       {
  //         out[i+4] += 1.f;
  //       }
  //     }
  // }

  return out;
}

void MeshCubeEquiUV::InitImpl(void)
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
