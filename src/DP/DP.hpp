#pragma once

#include "BFF/BFF.hpp"


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
    FaceData<Vector3> get_dp_wR() { return dp_wR; }
    VertexData<double> get_dp_div_wR_x() { return div_wR_x; }
    VertexData<double> get_dp_div_wR_y() { return div_wR_y; }

private:
    VertexData<Vector2> bff_uv;
    VertexData<double> bff_u;
    VertexData<double> dp_u;
    VertexData<Vector2> dp_uv;
    FaceData<Vector3> dp_wR;
    VertexData<double> div_wR_x;
    VertexData<double> div_wR_y;
};

} // namespace parameterization_playground
