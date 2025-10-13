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

    FaceData<Vector3> dp_euR_X(mesh), dp_euR_Y(mesh);
    Eigen::VectorXd dp_euR_X_vec(mesh.nFaces() * 3), dp_euR_Y_vec(mesh.nFaces() * 3);

    for (Face f : mesh.faces()) {
        double u_bar = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            u_bar += PD_u[he.tipVertex()];
        }
        u_bar /= 3.0;
        double eu = exp(u_bar);

        Eigen::Matrix2d euRT = eu * R[f].transpose();

        Eigen::Vector2<double> euRT_X = euRT.col(0);
        Eigen::Vector2<double> euRT_Y = euRT.col(1);

        dp_euR_X[f].x = euRT_X.x();
        dp_euR_X[f].y = euRT_X.y();
        dp_euR_X[f].z = 0;

        dp_euR_Y[f].x = euRT_Y.x();
        dp_euR_Y[f].y = euRT_Y.y();
        dp_euR_Y[f].z = 0;

        int fid = f.getIndex();
        dp_euR_X_vec[fid * 3] = euRT_X.x();
        dp_euR_X_vec[fid * 3 + 1] = euRT_X.y();
        dp_euR_X_vec[fid * 3 + 2] = 0;

        dp_euR_Y_vec[fid * 3] = euRT_Y.x();
        dp_euR_Y_vec[fid * 3 + 1] = euRT_Y.y();
        dp_euR_Y_vec[fid * 3 + 2] = 0;
    }


    this->dp_euRT_X = dp_euR_X;
    this->dp_euRT_Y = dp_euR_Y;

    geom.requirePolygonDivergenceMatrix();
    geom.requireCotanLaplacian();
    SparseMatrix<double> div_m = geom.polygonDivergenceMatrix;
    SparseMatrix<double> lap_m = geom.cotanLaplacian;

    int fixed_id = 0;
    for (int i = 0; i < lap_m.cols(); i++) {
        lap_m.coeffRef(fixed_id, i) = 0;
    }
    lap_m.coeffRef(fixed_id, fixed_id) = 1.0;

    Eigen::SimplicialLLT<SparseMatrix<double>> cholesky_lap;
    cholesky_lap.compute(lap_m);

    Eigen::VectorXd rhs_u = div_m * dp_euR_X_vec;
    Eigen::VectorXd rhs_v = div_m * dp_euR_Y_vec;

    VertexData<double> div_ueRT_X(mesh), div_ueRT_Y(mesh);

    for (Vertex v : mesh.vertices()) {
        int vid = v.getIndex();
        div_ueRT_X[v] = rhs_u[vid];
        div_ueRT_Y[v] = rhs_v[vid];
    }

    this->div_ueRT_X = div_ueRT_X;
    this->div_ueRT_Y = div_ueRT_Y;

    rhs_u(fixed_id) = 0;
    rhs_v(fixed_id) = 0;

    Eigen::VectorXd PD_uv_u = cholesky_lap.solve(rhs_u);
    Eigen::VectorXd PD_uv_v = cholesky_lap.solve(rhs_v);

    VertexData<Vector2> ret_uv(mesh);

    for (Vertex v : mesh.vertices()) {
        int vid = v.getIndex();
        ret_uv[v].x = PD_uv_u[vid];
        ret_uv[v].y = PD_uv_v[vid];
    }

    this->dp_uv = ret_uv;
    return ret_uv;
}

} // namespace parameterization_playground