#pragma once

#include "DP/DP.hpp"

#include "unsupported/Eigen/MatrixFunctions"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace parameterization_playground {

// circle mapping
class CM
{
public:
    CM();
    ~CM();

    VertexData<Vector2> solve(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        ManifoldSurfaceMesh& tar_mesh,
        VertexPositionGeometry& tar_geom);

private:
};

} // namespace parameterization_playground
