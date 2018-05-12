// Author: Xavier Corbillon
// IMT Atlantique
//
// Description:
// Implementation of a MeshEquiAngularCube: a mesh that represent a cube
#pragma once

//Internal includes
#include "Mesh.hpp"

namespace IMT {

class MeshEquiAngularCube: public Mesh
{
public:
  MeshEquiAngularCube(GLfloat scale, size_t numTriangles = 6*2*15*15);

  virtual ~MeshEquiAngularCube() = default;

private:
  MeshEquiAngularCube(const MeshEquiAngularCube&) = delete;
  MeshEquiAngularCube& operator=(const MeshEquiAngularCube&) = delete;

  // Swizzle each triple of coordinates by the specified
  // index and then multiply by the specified scale.  This
  // lets us implement a poor-man's rotation matrix, where
  // we pick which element (0-2) and which polarity (-1 or
  // 1) to use.
  std::vector<GLfloat> VertexRotate(
      GLfloat const &faceId,
      std::vector<GLfloat> const &inVec,
      std::array<size_t, 3> const &indices,
      std::array<GLfloat, 3> const &scales);

  std::vector<GLfloat> GetVertexs( std::vector<GLfloat> const& inputVertexs);
  std::vector<GLfloat> VertexToUVs( std::vector<GLfloat> const& inputVertexs);

  // // return UV map for each vertex of the quad of normalized index (i,j) for face f.
  // std::vector<GLfloat> GetUVs(float i, float j, float numQuadsPerEdge);

  // std::vector<GLfloat> TransposeUVs( std::vector<GLfloat> const& inputUVs,
  //     size_t faceId );

  virtual void InitImpl(void) override;
};
}
