#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <Eigen/SparseCholesky>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace parameterization_playground {

class BFF
{
public:
    BFF();
    ~BFF();

    enum ConstraintType { BND_SCALE, BND_CURVATURE, BND_FIT_SHAPE };

    VertexData<Vector2> solve(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& constraint,
        ConstraintType type);

    VertexData<Vector2>
    solve(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, std::vector<Vector2>& shape);

    VertexData<Vector2> solve(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& constraint,
        std::vector<Vector2>& shape,
        ConstraintType type);

    VertexData<double> get_u() { return result_u; }
    VertexData<double> get_k() { return curvature; }

private:
    VertexData<double> result_u;
    VertexData<double> curvature;
};

} // namespace parameterization_playground
