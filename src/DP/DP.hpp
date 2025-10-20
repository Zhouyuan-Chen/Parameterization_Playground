#pragma once

#include "BFF/BFF.hpp"

#include "unsupported/Eigen/MatrixFunctions"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace parameterization_playground {

class DP
{
public:
    DP();
    ~DP();

    VertexData<Vector2> solve(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        std::vector<Vector2>& deformed_shape);

    VertexData<Vector2> get_bff_uv() { return bff_uv; }
    VertexData<double> get_bff_u() { return bff_u; }
    VertexData<double> get_dp_u() { return dp_u; }
    VertexData<Vector2> get_dp_uv() { return dp_uv; }
    FaceData<Vector3> get_dp_euRT_X() { return dp_euRT_X; }
    FaceData<Vector3> get_dp_euRT_Y() { return dp_euRT_Y; }
    VertexData<double> get_dp_div_ueRT_X() { return div_ueRT_X; }
    VertexData<double> get_dp_div_ueRT_Y() { return div_ueRT_Y; }
    VertexData<double> get_dp_lap_u() { return dp_lap_u; }

    FaceData<double> get_test() { return test; }

private:
    VertexData<Vector2> bff_uv;
    VertexData<double> bff_u;
    VertexData<double> dp_u;
    VertexData<Vector2> dp_uv;
    FaceData<Vector3> dp_euRT_X;
    FaceData<Vector3> dp_euRT_Y;
    VertexData<double> div_ueRT_X;
    VertexData<double> div_ueRT_Y;
    VertexData<double> dp_lap_u;

    FaceData<double> test;
};

} // namespace parameterization_playground
