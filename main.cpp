#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "BFF/BFF.hpp"

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main(int argc, char** argv)
{
    polyscope::init();

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " mesh.obj" << " type" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string type = argv[2];
    // Read in a polygon mesh
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geom;
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);

    parameterization_playground::BFF bff;

    VertexData<double> scales(*mesh, 0);
    VertexData<double> curvatures(*mesh, 0);
    double bnd_num = 0;
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            bnd_num++;
        }
    }
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            curvatures[v] = 2 * PI / bnd_num;
        }
    }

    VertexData<Vector2> uv;
    if (type == "1") {
        uv = bff.solve(*mesh, *geom, scales, parameterization_playground::BFF::BND_SCALE);
    } else {
        uv = bff.solve(*mesh, *geom, curvatures, parameterization_playground::BFF::BND_CURVATURE);
    }

    polyscope::registerSurfaceMesh("3D domain", geom->vertexPositions, mesh->getFaceVertexList());

    polyscope::registerSurfaceMesh2D("2D domain", uv, mesh->getFaceVertexList());

    polyscope::getSurfaceMesh("3D domain")->addVertexParameterizationQuantity("uv", uv);
    polyscope::getSurfaceMesh("2D domain")->addVertexParameterizationQuantity("uv", uv);

    polyscope::show();

    return 0;
}