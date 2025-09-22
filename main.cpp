#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

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

    polyscope::registerSurfaceMesh("input mesh", geom->vertexPositions, mesh->getFaceVertexList());

    polyscope::show();

    return 0;
}