#include "DP.hpp"
#include "unsupported/Eigen/MatrixFunctions"

namespace parameterization_playground {
DP::DP() {}
DP::~DP() {}

VertexData<Vector2> DP::solve(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    std::vector<Vector2>& deformed_shape)
{
    BFF bff;
    VertexData<Vector2> bff_uv = bff.solve(mesh, geom, deformed_shape);
    this->bff_uv = bff_uv;

    VertexData<double> bff_u = bff.get_u();
    this->bff_u = bff_u;

    std::cout << "here" << std::endl;

    VertexData<Vector2> uv(mesh);
    VertexData<int> bvid(mesh);
    VertexData<int> ivid(mesh);
    auto isInterior = Vector<bool>(mesh.nVertices());
    int bid = 0, iid = 0;
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            isInterior[bid + iid] = false;
            bvid[v] = bid++;
        } else {
            isInterior[bid + iid] = true;
            ivid[v] = iid++;
        }
    }
    geom.requireEdgeLengths();
    geom.requireCotanLaplacian();
    SparseMatrix<double> A(mesh.nVertices(), mesh.nVertices());
    SparseMatrix<double> L = geom.cotanLaplacian;
    shiftDiagonal(L, 1e-10);
    auto Adecomp = blockDecomposeSquare(geom.cotanLaplacian, isInterior);
    auto AII = Adecomp.AA;
    auto AIB = Adecomp.AB;
    auto ABI = Adecomp.BA;
    auto ABB = Adecomp.BB;
    Eigen::SimplicialLLT<SparseMatrix<double>> cholesky_AII;
    cholesky_AII.compute(AII);

    Eigen::VectorXd boundary_u_vec(bid);
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            int id = bvid[v];
            boundary_u_vec[id] = bff_u[v];
        }
    }

    Eigen::VectorXd interior_u_vec = cholesky_AII.solve(-AIB * boundary_u_vec);

    VertexData<double> PD_u(mesh);
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            int id = bvid[v];
            PD_u[v] = boundary_u_vec[id];
        } else {
            int id = ivid[v];
            PD_u[v] = interior_u_vec[id];
        }
    }

    this->dp_u = PD_u;

    geom.requireFaceNormals();
    geom.requireEdgeCotanWeights();
    EdgeData<double> w12_list(mesh);
    for (Edge edge : mesh.edges()) {
        if (!edge.isBoundary()) {
            Halfedge he = edge.halfedge();
            double ui = PD_u[he.tailVertex()];
            double uj = PD_u[he.tipVertex()];
            double cot_w = geom.edgeCotanWeight(edge);
            w12_list[edge] = cot_w * (uj - ui);
        }
    }

    std::cout << "here" << std::endl;

    EdgeData<Eigen::Matrix2<double>> w(mesh);
    for (Edge edge : mesh.edges()) {
        if (!edge.isBoundary()) {
            double w12_int = w12_list[edge];
            w[edge](0, 0) = 0;
            w[edge](0, 1) = w12_int;
            w[edge](1, 0) = -w12_int;
            w[edge](1, 1) = 0;
        }
    }
    std::cout << "here" << std::endl;

    // only work for 2D for now
    FaceData<Eigen::Matrix2<double>> R(mesh);
    R[mesh.face(229)] = Eigen::Matrix2<double>::Identity();
    FaceData<bool> visited_tag(mesh, false);
    visited_tag[mesh.face(229)] = true;

    while (true) {
        bool all_assigned = true;
        for (Edge edge : mesh.edges()) {
            if (!edge.isBoundary()) {
                Halfedge he = edge.halfedge();
                if (!visited_tag[he.face()] && visited_tag[he.twin().face()]) {
                    all_assigned = false;
                    Eigen::Matrix2d Rt = R[he.twin().face()];
                    Eigen::Matrix2d w_int = w[he.edge()];
                    R[he.face()] = (-w_int).exp() * Rt;
                    visited_tag[he.face()] = true;
                } else if (visited_tag[he.face()] && !visited_tag[he.twin().face()]) {
                    all_assigned = false;
                    Eigen::Matrix2<double> R0 = R[he.face()];
                    Eigen::Matrix2<double> w_int = w[he.edge()];
                    R[he.twin().face()] = w_int.exp() * R0;
                    visited_tag[he.twin().face()] = true;
                }
            }
        }
        if (all_assigned) {
            break;
        }
    }

    FaceData<Vector3> dir(mesh);

    for (Face f : mesh.faces()) {
        double w = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
        }
        Eigen::Matrix2d Rf = R[f].transpose();
        Eigen::Vector2d v = Rf.col(0);
        // v.normalize();
        dir[f].x = v.x();
        dir[f].y = v.y();
        dir[f].z = 0;
    }

    this->dp_wR = dir;

    // Eigen::SparseMatrix<double> grad_m = geom.polygonGradientMatrix;

    Eigen::VectorXd frame(3 * mesh.nFaces());

    for (Face f : mesh.faces()) {
        int id = f.getIndex();
        double eu = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            eu += PD_u[he.tipVertex()];
        }
        eu /= 3.0;

        frame[id * 3] = dir[f].x;
        frame[id * 3 + 1] = dir[f].y;
        frame[id * 3 + 2] = 0;
    }

    geom.requirePolygonDivergenceMatrix();
    geom.requireCotanLaplacian();

    auto lap_m = geom.cotanLaplacian;
    auto div_m = geom.polygonDivergenceMatrix;

    Eigen::VectorXd frame_x(mesh.nFaces() * 3);
    Eigen::VectorXd frame_y(mesh.nFaces() * 3);
    for (Face f : mesh.faces()) {
        int id = f.getIndex();
        frame_x(id * 3) = frame[id * 3];
        frame_x(id * 3 + 1) = 0;
        frame_x(id * 3 + 2) = 0;
        frame_y(id * 3) = frame[id * 3 + 1];
        frame_y(id * 3 + 1) = 0;
        frame_y(id * 3 + 2) = 0;
    }

    std::cout << frame_x.size() << std::endl;
    std::cout << div_m.rows() << " " << div_m.cols() << std::endl;


    Eigen::VectorXd div_x = div_m * frame_x;
    Eigen::VectorXd div_y = div_m * frame_y;


    VertexData<double> div_wR_x(mesh);
    VertexData<double> div_wR_y(mesh);
    for (Vertex v : mesh.vertices()) {
        int id = v.getIndex();
        div_wR_x[v] = div_x(id);
        div_wR_y[v] = div_y(id);
    }

    this->div_wR_x = div_wR_x;
    this->div_wR_y = div_wR_y;

    int fixed_id = 0;
    for (int i = 0; i < lap_m.cols(); i++) {
        lap_m.coeffRef(fixed_id, i) = 0;
    }
    lap_m.coeffRef(fixed_id, fixed_id) = 1.0;
    div_x(fixed_id) = 0.0;
    div_y(fixed_id) = 0.0;

    Eigen::SimplicialLLT<SparseMatrix<double>> cholesky_uv;
    cholesky_uv.compute(lap_m);


    std::cout << lap_m.rows() << " " << lap_m.cols() << std::endl;

    Eigen::VectorXd u = cholesky_uv.solve(div_x);
    Eigen::VectorXd v = cholesky_uv.solve(div_y);
    std::cout << "finished" << std::endl;

    VertexData<Vector2> ret_uv(mesh);
    for (Vertex vtx : mesh.vertices()) {
        int id = vtx.getIndex();
        ret_uv[vtx] = Vector2({u(id), v(id)});
    }

    this->dp_uv = ret_uv;
    return ret_uv;
}

} // namespace parameterization_playground