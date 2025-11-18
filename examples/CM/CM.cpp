#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "igl/AABB.h"

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

    std::string filename1 = argv[1];
    std::string filename2 = argv[2];
    // Read in a polygon mesh
    std::unique_ptr<ManifoldSurfaceMesh> mesh, tar_mesh;
    std::unique_ptr<VertexPositionGeometry> geom, tar_geom;
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename1);
    std::tie(tar_mesh, tar_geom) = readManifoldSurfaceMesh(filename2);

    double bnd_num = 0;
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            bnd_num++;
        }
    }
    // Disk
    std::vector<Vector2> shape;
    for (int i = 0; i < bnd_num; i++) {
        double t = (i * 1.0) / (bnd_num * 1.0);
        shape.push_back(Vector2({cos(t * 2 * M_PI), sin(t * 2 * M_PI)}));
    }

    parameterization_playground::DP dp;
    auto dp_uv = dp.solve(*mesh, *geom, shape);
    parameterization_playground::DP dp_target;
    auto dp_uv_target = dp_target.solve(*tar_mesh, *tar_geom, shape);

    dp_uv = dp.get_bff_uv();
    dp_uv_target = dp_target.get_bff_uv();
    // polyscope::registerSurfaceMesh("3D domain", geom->vertexPositions,
    // mesh->getFaceVertexList()); polyscope::getSurfaceMesh("3D
    // domain")->addVertexParameterizationQuantity("uv_dp", dp_uv); polyscope::getSurfaceMesh("3D
    // domain")
    //     ->addVertexParameterizationQuantity("uv_bff", dp.get_bff_uv());
    // polyscope::getSurfaceMesh("3D domain")->addVertexScalarQuantity("u", dp.get_dp_u());
    // polyscope::getSurfaceMesh("3D domain")
    //     ->addVertexScalarQuantity("wR_div_x", dp.get_dp_div_ueRT_X());
    // polyscope::getSurfaceMesh("3D domain")
    //     ->addVertexScalarQuantity("wR_div_y", dp.get_dp_div_ueRT_Y());
    // polyscope::getSurfaceMesh("3D domain")->addFaceVectorQuantity("euRT_X", dp.get_dp_euRT_X());
    // polyscope::getSurfaceMesh("3D domain")->addFaceVectorQuantity("euRT_Y", dp.get_dp_euRT_Y());
    // polyscope::getSurfaceMesh("3D domain")->addVertexScalarQuantity("lap_u", dp.get_dp_lap_u());

    // normalize
    // circle radius=1, center at (0,0)
    // mesh
    double u_bar, v_bar, max_v, min_v, bnd_bar;
    u_bar = 0;
    v_bar = 0;
    max_v = -999999;
    min_v = 999999;
    bnd_bar = 0;
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            bnd_bar++;
            u_bar += dp_uv[v].x;
            v_bar += dp_uv[v].y;
            if (dp_uv[v].y < min_v) {
                min_v = dp_uv[v].y;
            }
            if (dp_uv[v].y > max_v) {
                max_v = dp_uv[v].y;
            }
        }
    }
    double ratio = (max_v - min_v) * 0.5;
    u_bar /= bnd_bar;
    v_bar /= bnd_bar;
    for (Vertex v : mesh->vertices()) {
        dp_uv[v].x = (dp_uv[v].x - u_bar) / ratio;
        dp_uv[v].y = (dp_uv[v].y - v_bar) / ratio;
    }

    std::cout << "here" << std::endl;

    // mesh target
    double u_t_bar, v_t_bar, max_t_v, min_t_v, bnd_t_bar;
    u_t_bar = 0;
    v_t_bar = 0;
    max_t_v = -999999;
    min_t_v = 999999;
    bnd_t_bar = 0;
    for (Vertex v : tar_mesh->vertices()) {
        if (v.isBoundary()) {
            bnd_t_bar++;
            u_t_bar += dp_uv_target[v].x;
            v_t_bar += dp_uv_target[v].y;
            if (dp_uv_target[v].y < min_t_v) {
                min_t_v = dp_uv_target[v].y;
            }
            if (dp_uv_target[v].y > max_t_v) {
                max_t_v = dp_uv_target[v].y;
            }
        }
    }
    double ratio_t = (max_t_v - min_t_v) * 0.5;
    u_t_bar /= bnd_t_bar;
    v_t_bar /= bnd_t_bar;
    for (Vertex v : tar_mesh->vertices()) {
        dp_uv_target[v].x = (dp_uv_target[v].x - u_t_bar) / ratio_t + 0.01;
        dp_uv_target[v].y = (dp_uv_target[v].y - v_t_bar) / ratio_t - 0.015;
    }
    std::cout << "here" << std::endl;

    // mapping part
    // use igl::aabb to compute the mapping
    // ==================== AABB Mapping Part ====================

    // Convert parameterizations to Eigen matrices for libigl
    Eigen::MatrixXd V1(mesh->nVertices(), 2);
    Eigen::MatrixXi F1(mesh->nFaces(), 3);
    Eigen::MatrixXd V2(tar_mesh->nVertices(), 2);
    Eigen::MatrixXi F2(tar_mesh->nFaces(), 3);

    // Fill source mesh (mesh1)
    for (Vertex v : mesh->vertices()) {
        V1(v.getIndex(), 0) = dp_uv[v].x;
        V1(v.getIndex(), 1) = dp_uv[v].y;
    }
    for (Face f : mesh->faces()) {
        int i = 0;
        for (Vertex v : f.adjacentVertices()) {
            F1(f.getIndex(), i++) = v.getIndex();
        }
    }

    // Fill target mesh (mesh2)
    for (Vertex v : tar_mesh->vertices()) {
        V2(v.getIndex(), 0) = dp_uv_target[v].x;
        V2(v.getIndex(), 1) = dp_uv_target[v].y;
    }
    for (Face f : tar_mesh->faces()) {
        int i = 0;
        for (Vertex v : f.adjacentVertices()) {
            F2(f.getIndex(), i++) = v.getIndex();
        }
    }

    // Build AABB tree for target mesh
    igl::AABB<Eigen::MatrixXd, 2> tree;
    tree.init(V2, F2);

    // For each vertex in source mesh, find closest point on target mesh
    Eigen::VectorXi I; // Face indices in target mesh
    Eigen::MatrixXd C; // Closest points in parameter space
    Eigen::VectorXd D; // Squared distances

    // Find closest points - correct function signature
    tree.squared_distance(V2, F2, V1, D, I, C);

    // Now compute barycentric coordinates manually
    std::vector<Vector3> barycentric_coords(mesh->nVertices());

    for (int i = 0; i < mesh->nVertices(); i++) {
        int target_face_index = I(i);
        Face target_face = tar_mesh->face(target_face_index);

        // Get the three vertices of the target face
        std::vector<Vertex> face_vertices;
        for (Vertex v : target_face.adjacentVertices()) {
            face_vertices.push_back(v);
        }

        // Get the 2D parameter positions of the target face vertices
        Vector2 p0 = dp_uv_target[face_vertices[0]];
        Vector2 p1 = dp_uv_target[face_vertices[1]];
        Vector2 p2 = dp_uv_target[face_vertices[2]];
        Vector2 query_point = Vector2{C(i, 0), C(i, 1)};

        // Compute barycentric coordinates manually
        Vector2 v0 = p1 - p0, v1 = p2 - p0, v2 = query_point - p0;
        double d00 = dot(v0, v0);
        double d01 = dot(v0, v1);
        double d11 = dot(v1, v1);
        double d20 = dot(v2, v0);
        double d21 = dot(v2, v1);
        double denom = d00 * d11 - d01 * d01;

        double beta = (d11 * d20 - d01 * d21) / denom;
        double gamma = (d00 * d21 - d01 * d20) / denom;
        double alpha = 1.0 - beta - gamma;

        barycentric_coords[i] = Vector3{alpha, beta, gamma};
    }

    // Compute mapping results - map back to original 3D target mesh
    std::vector<Vector3> mapped_3d_points(mesh->nVertices());
    std::vector<double> distances(mesh->nVertices());
    VertexData<Vector2> mapped_3d(*mesh);

    for (int i = 0; i < mesh->nVertices(); i++) {
        int target_face_index = I(i);
        Face target_face = tar_mesh->face(target_face_index);

        // Get the three vertices of the target face
        std::vector<Vertex> face_vertices;
        for (Vertex v : target_face.adjacentVertices()) {
            face_vertices.push_back(v);
        }

        // Get original 3D positions of the target face vertices
        Vector3 v0 = tar_geom->vertexPositions[face_vertices[0]];
        Vector3 v1 = tar_geom->vertexPositions[face_vertices[1]];
        Vector3 v2 = tar_geom->vertexPositions[face_vertices[2]];

        // Get barycentric coordinates
        double alpha = barycentric_coords[i].x;
        double beta = barycentric_coords[i].y;
        double gamma = barycentric_coords[i].z;

        // Interpolate to get 3D position on target mesh
        Vector3 interpolated_pos = alpha * v0 + beta * v1 + gamma * v2;

        mapped_3d_points[i] = interpolated_pos;
        distances[i] = std::sqrt(D(i));
        mapped_3d[mesh->vertex(i)] = Vector2({interpolated_pos.x, interpolated_pos.y});
    }

    std::cout << "AABB mapping with barycentric interpolation completed" << std::endl;

    // ==================== Visualization ====================

    // Register original 3D meshes
    polyscope::registerSurfaceMesh(
        "Source 3D mesh",
        geom->vertexPositions,
        mesh->getFaceVertexList());
    polyscope::getSurfaceMesh("Source 3D mesh")
        ->addVertexScalarQuantity("Parameter space distance", distances);
    polyscope::getSurfaceMesh("Source 3D mesh")
        ->addVertexParameterizationQuantity("uv_to_target", mapped_3d);

    polyscope::registerSurfaceMesh(
        "Target 3D mesh",
        tar_geom->vertexPositions,
        tar_mesh->getFaceVertexList());

    // Register 2D parameterizations
    polyscope::registerSurfaceMesh2D("Source 2D domain", dp_uv, mesh->getFaceVertexList());
    polyscope::getSurfaceMesh("Source 2D domain")->addVertexParameterizationQuantity("uv", dp_uv);

    polyscope::registerSurfaceMesh2D(
        "Target 2D domain",
        dp_uv_target,
        tar_mesh->getFaceVertexList());
    polyscope::getSurfaceMesh("Target 2D domain")
        ->addVertexParameterizationQuantity("uv", dp_uv_target);

    // Visualize mapping results in 3D
    // polyscope::getSurfaceMesh("Source 3D mesh")
    //     ->addVertexVectorQuantity("Mapping vectors to target", mapped_3d_points);
    polyscope::registerSurfaceMesh(
        "Source 3D mesh_mapping to target",
        mapped_3d_points,
        mesh->getFaceVertexList());


    // polyscope::registerSurfaceMesh2D("2D domain", dp_uv, mesh->getFaceVertexList());
    // polyscope::getSurfaceMesh("2D domain")->addVertexParameterizationQuantity("uv", dp_uv);
    // // polyscope::getSurfaceMesh("2D domain")->addVertexScalarQuantity("u", dp.get_dp_u());
    // // polyscope::getSurfaceMesh("2D domain")->addVertexScalarQuantity("lap_u", dp.get_dp_lap_u());

    // polyscope::registerSurfaceMesh2D(
    //     "target 2D domain",
    //     dp_uv_target,
    //     tar_mesh->getFaceVertexList());
    // polyscope::getSurfaceMesh("target 2D domain")
    //     ->addVertexParameterizationQuantity("uv", dp_uv_target);


    polyscope::show();

    return 0;
}