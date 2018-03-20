//Author: Xavier Corbillon
//IMT Atlantique
#include "MeshSegmentedSphere.hpp"
#include <iostream>

//standard includes
#include <cmath>

using namespace IMT;
constexpr float PI = 3.14159265358979323846L;
constexpr float PI_2 = 1.57079632679489661923L;
constexpr float square = 1040.0L;
constexpr float U_BORDER = 1.6e-4;
constexpr float V_BORDER = 9.6e-4;

MeshSegmentedSphere::MeshSegmentedSphere(GLfloat scale, size_t numTriangles): Mesh()
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
  // int slices = 360*5;
  // int layercuts = 179*5;
  // float x, y, z, d; //d=diameter of circle at height y
  // float theta; //theta used for angle about y axis
  // float phi = 0.0f; //phi used for angle to determine y and d.
  // float thetad = 2 * PI / slices; //theta delta
  // float phid = PI / (layercuts + 1); //phi delta

  // // triangle meshes from north pole to the first layercut
  // GLfloat north_pole[3] = {0.0f, 1.0f, 0.0f};
  // phi = phid;
  // y = std::cos(phi);
  // d = std::sin(phi);
  // theta = 0.0f;
  // for (unsigned int i = 0; i < slices; i++) {
  //   auto x1 = std::sin(theta) * d;
  //   auto z1 = std::cos(theta) * d;
  //   auto x2 = std::sin(theta+thetad) * d;
  //   auto z2 = std::cos(theta+thetad) * d;
  //   if (i == slices-1)
  //   {
  //     x2 = std::sin(0.0f) * d;
  //     z2 = std::cos(0.0f) * d;
  //   }
  //   GLfloat tmpVertex[6] = {x1, y, z1, x2, y, z2};
  //   tmpVertexBufferData.insert(tmpVertexBufferData.end(), north_pole, north_pole+3);
  //   tmpVertexBufferData.insert(tmpVertexBufferData.end(), tmpVertex, tmpVertex+6);
  //   theta += thetad;
  // }

  // // triangle meshes between the first and the last layercuts 
  // phi = phid;
  // for (unsigned int i = 0; i < layercuts-1; i++) {
  //   auto y1 = std::cos(phi);
  //   auto d1 = std::sin(phi);
  //   auto y2 = std::cos(phi+phid);
  //   auto d2 = std::sin(phi+phid);
  //   theta = 0.0f;
  //   for (unsigned int j = 0; j < slices; j++) {
  //     auto x1 = std::sin(theta) * d1;
  //     auto z1 = std::cos(theta) * d1;
  //     auto x2 = std::sin(theta) * d2;
  //     auto z2 = std::cos(theta) * d2;
  //     auto x3 = std::sin(theta+thetad) * d2;
  //     auto z3 = std::cos(theta+thetad) * d2;
  //     if (i == slices-1)
  //     {
  //       x3 = std::sin(0.0f) * d2;
  //       z3 = std::cos(0.0f) * d2;
  //     }
  //     // u = 0.5f + (std::atan2(z, x) / (2 * PI));
  //     // v = 0.5f - (std::asin(y) / PI);
  //     // vertices[(i * slices) + j + 1] = { { x / 2, y / 2, z / 2 },{ u, v },{ 1.0f, 1.0f, 1.0f },{ x, y, z } };
  //     GLfloat tmpVertex[9] = {x1, y1, z1, x2, y2, z2, x3, y2, z3};
  //     tmpVertexBufferData.insert(tmpVertexBufferData.end(), tmpVertex, tmpVertex+9);
  //     theta += thetad;
  //   }
  //   phi += phid;
  // }
  // phi = layercuts * phid;
  // for (unsigned int i = 0; i < layercuts-1; i++) {
  //   auto y1 = std::cos(phi);
  //   auto d1 = std::sin(phi);
  //   auto y2 = std::cos(phi-phid);
  //   auto d2 = std::sin(phi-phid);
  //   theta = 0.0f;
  //   for (unsigned int j = 0; j < slices; j++) {
  //     auto x1 = std::sin(theta) * d1;
  //     auto z1 = std::cos(theta) * d1;
  //     auto x2 = std::sin(theta) * d2;
  //     auto z2 = std::cos(theta) * d2;
  //     auto x3 = std::sin(theta-thetad) * d2;
  //     auto z3 = std::cos(theta-thetad) * d2;
  //     if (i == slices-1)
  //     {
  //       x3 = std::sin(0.0f) * d2;
  //       z3 = std::cos(0.0f) * d2;
  //     }
  //     // u = 0.5f + (std::atan2(z, x) / (2 * PI));
  //     // v = 0.5f - (std::asin(y) / PI);
  //     // vertices[(i * slices) + j + 1] = { { x / 2, y / 2, z / 2 },{ u, v },{ 1.0f, 1.0f, 1.0f },{ x, y, z } };
  //     GLfloat tmpVertex[9] = {x1, y1, z1, x2, y2, z2, x3, y2, z3};
  //     tmpVertexBufferData.insert(tmpVertexBufferData.end(), tmpVertex, tmpVertex+9);
  //     theta -= thetad;
  //   }
  //   phi -= phid;
  // }

  // // triangle meshes from south pole to the last layercut
  // GLfloat south_pole[3] = {0.0f, -1.0f, 0.0f};
  // phi = layercuts * phid;
  // y = std::cos(phi);
  // d = std::sin(phi);
  // theta = 0.0f;
  // for (unsigned int i = 0; i < slices; i++) {
  //   auto x1 = std::sin(theta) * d;
  //   auto z1 = std::cos(theta) * d;
  //   auto x2 = std::sin(theta+thetad) * d;
  //   auto z2 = std::cos(theta+thetad) * d;
  //   if (i == slices-1)
  //   {
  //     x2 = std::sin(0.0f) * d;
  //     z2 = std::cos(0.0f) * d;
  //   }
  //   GLfloat tmpVertex[6] = {x1, y, z1, x2, y, z2};
  //   tmpVertexBufferData.insert(tmpVertexBufferData.end(), tmpVertex, tmpVertex+6);
  //   tmpVertexBufferData.insert(tmpVertexBufferData.end(), south_pole, south_pole+3);
  //   theta += thetad;
  // }

  // // for (int i = 0; i < tmpVertexBufferData.size(); i+=3)
  // // {
  // //   std::cout << "(x, y, z): (" << tmpVertexBufferData[i] << ", " << tmpVertexBufferData[i+1] << ", " << tmpVertexBufferData[i+2] << ")" << std::endl;
  // // }
  // AppendVertexBufferData(tmpVertexBufferData);
  // AppendUvBufferData(VertexToUVs(tmpVertexBufferData));
}


std::vector<GLfloat> MeshSegmentedSphere::VertexRotate(
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

std::vector<GLfloat> MeshSegmentedSphere::VertexToUVs( std::vector<GLfloat> const& inputVertexs)
{
  std::vector<GLfloat> out;

  float max_u = -10000.0;
  float min_u = 10000.0;
  float max_v = -10000.0;
  float min_v = 10000.0;

  //Generate SegmentedSphere UV map
  for(size_t i = 0; i < inputVertexs.size(); i += 3)
  {
    auto& x = inputVertexs[i];
    auto& y = inputVertexs[i+1];
    auto& z = inputVertexs[i+2];
    auto rho = std::sqrt(x*x+y*y+z*z);
    auto theta = std::atan2(z, x)+PI;
    auto phi = std::asin(y / rho);
    auto u = 0.0f, v = 0.0f;
    auto faceid = 0;

    if (y > rho*std::sqrt(2) / 2)
    {
      faceid = 0;
      u = ((square * std::sin(theta) * (PI_2 - phi) / PI_2 + square / 2)) / 6240.0;
      v = ((square / 2 * (1 + std::cos(theta) * 2 * (PI_2 - phi) / PI_2))) / 1040.0;
      // u = (square * std::sin(theta) * (PI_2 - phi) / PI_2 + square / 2);
      // u = (u < 1? 0.5 : (u-0.5)) / 6241.0;
      // v = ((square / 2 * (1 + std::cos(theta) * 2 * (PI_2 - phi) / PI_2)) + 0.5) / 1041.0;
      
    }
    else if (y < -rho*std::sqrt(2) / 2)
    {
      faceid = 1;
      u = ((square * std::sin(theta) * (PI_2 + phi) / PI_2 + square / 2) + square) / 6240.0;
      v = ((square / 2 * (1 - std::cos(theta) * 2 * (PI_2 + phi) / PI_2))) / 1040.0;
      // u = (square * std::sin(theta) * (PI_2 + phi) / PI_2 + square / 2);
      // u = ((u < 1? 0.5 : (u-0.5)) + 1040.0) / 6241.0;
      // v = ((square / 2 * (1 - std::cos(theta) * 2 * (PI_2 + phi) / PI_2)) + 0.5) / 1041.0;
      
    }
    else if (z >= 0 && x < 0)
    {
      faceid = 2;
      u = ((PI + theta)*square / PI_2 - 2*square) / 6240.0;
      v = (((PI_2 - phi) * square / PI_2 - square / 2)) / 1040.0;
      // u = ((PI + theta)*square / PI_2);
      // u = ((u < 1? 0.5 : (u-0.5)) + 2*1040.0) / 6241.0;
      // v = (((PI_2 - phi) * square / PI_2 - square / 2) + 0.5) / 1041.0;
    }
    else if (z > 0 && x >= 0)
    {
      faceid = 3;
      u = ((PI + theta)*square / PI_2 - 2*square) / 6240.0;
      v = (((PI_2 - phi) * square / PI_2 - square / 2)) / 1040.0;
      // u = ((PI + theta)*square / PI_2 - square);
      // u = ((u < 1? 0.5 : (u-0.5)) + 3*1040.0) / 6241.0;
      // v = (((PI_2 - phi) * square / PI_2 - square / 2) + 0.5) / 1041.0;
      u = (u <= (0.333537f+U_BORDER)? (0.333537f+U_BORDER) : u);
    }
    else if (z <= 0 && x > 0)
    {
      faceid = 4;
      u = ((PI + theta)*square / PI_2 + 2*square) / 6240.0;
      v = (((PI_2 - phi) * square / PI_2 - square / 2)) / 1040.0;
      // u = ((PI + theta)*square / PI_2 - 2 * square);
      // u = ((u < 1? 0.5 : (u-0.5)) + 4*1040.0) / 6241.0;
      // v = (((PI_2 - phi) * square / PI_2 - square / 2) + 0.5) / 1041.0;
      u = (u >= (1.0f-U_BORDER)? (1.0f-U_BORDER) : u);
    }
    else if (z < 0 && x <= 0)
    {
      faceid = 5;
      u = ((PI + theta)*square / PI_2 + 2*square) / 6240.0;
      v = (((PI_2 - phi) * square / PI_2 - square / 2)) / 1040.0;
      // u = ((PI + theta)*square / PI_2 - 3 * square);
      // u = ((u < 1? 0.5 : (u-0.5)) + 5*1040.0) / 6241.0;
      // v = (((PI_2 - phi) * square / PI_2 - square / 2) + 0.5) / 1041.0;
      
    }
    
    u = (u >= (1-U_BORDER)? (1-U_BORDER) : u);
    u = (u <= U_BORDER? U_BORDER : u);
    v = (v >= (1-V_BORDER)? (1-V_BORDER) : v);
    v = (v <= V_BORDER? V_BORDER : v);

    if (faceid == 6)
    {
      if (max_u < u)
      {
          max_u = u;
      }

      if (min_u > u)
      {
          min_u = u;
      }

      if (max_v < v)
      {
          max_v = v;
      }

      if (min_v > v)
      {
          min_v = v;
      }
      std::cout << "u:[" << min_u << ", " << max_u << "] v:[" << min_v << ", " << max_v << "]" << std::endl;
      u = 0.0f;
      v = 0.0f;
      
    }

    out.push_back(u);
    out.push_back(v);
  }

  return out;
}

void MeshSegmentedSphere::InitImpl(void)
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
