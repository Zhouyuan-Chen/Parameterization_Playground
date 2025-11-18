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

class RectangularParameterization
{
public:
    RectangularParameterization();

    void FlattenMesh(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom);
    void SolveOptimizationProblem(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom);
    void hessObjective(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& u,
        VertexData<double>& v,
        FaceData<double>& theta,
        Eigen::SparseMatrix<double>& H);
    void gradObjective(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& u,
        VertexData<double>& v,
        FaceData<double>& theta,
        Eigen::VectorXd& g);
    void BuildHessian(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& u,
        VertexData<double>& v,
        FaceData<double>& theta,
        EdgeData<double>& lambda,
        Eigen::SparseMatrix<double>& H);
    void BuildJacobian(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& u,
        VertexData<double>& v,
        FaceData<double>& theta,
        EdgeData<double>& lambda,
        Eigen::SparseMatrix<double>& J);
    void BuildSystem(
        ManifoldSurfaceMesh& mesh,
        VertexPositionGeometry& geom,
        VertexData<double>& u,
        VertexData<double>& v,
        FaceData<double>& theta,
        EdgeData<double>& lambda,
        Eigen::VectorXd& F);


    VertexData<double> ret_u;
    VertexData<double> ret_v;
    FaceData<double> ret_theta;
    EdgeData<double> ret_lambda;
};

} // namespace parameterization_playground