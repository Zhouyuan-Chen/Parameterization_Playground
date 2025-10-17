#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "DP/DP.hpp"

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main(int argc, char** argv)
{
    polyscope::init();

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " mesh.obj" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    // Read in a polygon mesh
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);

    double bnd_num = 0;
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            bnd_num++;
        }
    }
    // Disk
    // std::vector<Vector2> shape;
    // for (int i = 0; i < bnd_num; i++) {
    //     double t = (i * 1.0) / (bnd_num * 1.0);
    //     shape.push_back(Vector2({cos(t * 2 * M_PI), sin(t * 2 * M_PI)}));
    // }

    std::vector<Vector2> shape(
        {Vector2({0, 0}),
         Vector2({1, 0}),
         Vector2({2, 0}),
         Vector2({2, 1}),
         Vector2({2, 2}),
         Vector2({1, 2}),
         Vector2({1, 1}),
         Vector2({0, 1})});

    parameterization_playground::DP dp;
    auto dp_uv = dp.solve(*mesh, *geom, shape);

    FaceData<double> test(*mesh);
    double testi = 0;

    for (Face f : mesh->faces()) {
        test[f] = testi++;
    }

    polyscope::registerSurfaceMesh("3D domain", geom->vertexPositions, mesh->getFaceVertexList());
    polyscope::getSurfaceMesh("3D domain")->addVertexParameterizationQuantity("uv", dp_uv);
    polyscope::getSurfaceMesh("3D domain")->addVertexScalarQuantity("u", dp.get_dp_u());
    polyscope::getSurfaceMesh("3D domain")
        ->addVertexScalarQuantity("wR_div_x", dp.get_dp_div_ueRT_X());
    polyscope::getSurfaceMesh("3D domain")
        ->addVertexScalarQuantity("wR_div_y", dp.get_dp_div_ueRT_Y());
    polyscope::getSurfaceMesh("3D domain")->addFaceVectorQuantity("euRT_X", dp.get_dp_euRT_X());
    polyscope::getSurfaceMesh("3D domain")->addFaceVectorQuantity("euRT_Y", dp.get_dp_euRT_Y());


    polyscope::registerSurfaceMesh2D("2D domain", dp_uv, mesh->getFaceVertexList());
    polyscope::getSurfaceMesh("2D domain")->addVertexParameterizationQuantity("uv", dp_uv);
    polyscope::getSurfaceMesh("2D domain")->addVertexScalarQuantity("u", dp.get_dp_u());

    polyscope::show();

    return 0;
}