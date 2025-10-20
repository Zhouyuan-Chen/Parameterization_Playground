#include "DP.hpp"

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
        int vid = v.getIndex();
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
            // boundary_u_vec[id] = 0;
        }
    }

    VertexData<double> vertex_k = bff.get_k();
    Eigen::VectorXd int_k(iid);
    for (Vertex v : mesh.vertices()) {
        if (!v.isBoundary()) {
            int id = ivid[v];
            int_k[id] = vertex_k[v];
        }
    }

    Eigen::VectorXd interior_u_vec = cholesky_AII.solve(-AIB * boundary_u_vec - int_k);

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

    // PD_u = bff.get_u();

    this->dp_u = PD_u;

    Eigen::VectorXd u_vec(mesh.nVertices());
    for (Vertex v : mesh.vertices()) {
        u_vec[v.getIndex()] = PD_u[v];
    }
    Eigen::VectorXd lap_vec = geom.cotanLaplacian * u_vec;
    VertexData<double> lap_u(mesh);
    for (Vertex v : mesh.vertices()) {
        lap_u[v] = lap_vec[v.getIndex()];
    }
    this->dp_lap_u = lap_u;

    geom.requireFaceNormals();
    geom.requireEdgeCotanWeights();
    HalfedgeData<double> w12_list(mesh);
    for (Halfedge he : mesh.halfedges()) {
        if (!he.edge().isBoundary()) {
            double ui = PD_u[he.tailVertex()];
            double uj = PD_u[he.tipVertex()];
            double cot_w = geom.edgeCotanWeight(he.edge());
            w12_list[he] = cot_w * (uj - ui);
        }
    }

    HalfedgeData<Eigen::Matrix2<double>> w(mesh);
    for (Halfedge he : mesh.halfedges()) {
        if (!he.edge().isBoundary()) {
            double w12_int = w12_list[he];
            w[he](0, 0) = 0;
            w[he](0, 1) = w12_int;
            w[he](1, 0) = -w12_int;
            w[he](1, 1) = 0;
        }
    }

    geom.requireTransportVectorsAcrossHalfedge();

    geom.requireFaceTangentBasis();
    Eigen::SparseMatrix<double> R_face_basis(mesh.nEdges() * 2 + 2, mesh.nFaces() * 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(mesh.nEdges() * 8 + 4);
    for (Edge edge : mesh.edges()) {
        if (!edge.isBoundary()) {
            Halfedge he = edge.halfedge();
            Face f0 = he.face();
            Face f1 = he.twin().face();
            int fid0 = f0.getIndex();
            int fid1 = f1.getIndex();
            int eid = edge.getIndex();

            // Vector3 f0_bx = geom.faceTangentBasis[f0][0];
            // Vector3 f0_by = geom.faceTangentBasis[f0][1];
            // Vector3 f1_bx = geom.faceTangentBasis[f1][0];
            // Vector3 f1_by = geom.faceTangentBasis[f1][1];


            // triplets.emplace_back(eid * 3, fid0 * 2, f0_bx.x);
            // triplets.emplace_back(eid * 3 + 1, fid0 * 2, f0_bx.y);
            // triplets.emplace_back(eid * 3 + 2, fid0 * 2, f0_bx.z);
            // triplets.emplace_back(eid * 3, fid0 * 2 + 1, f0_by.x);
            // triplets.emplace_back(eid * 3 + 1, fid0 * 2 + 1, f0_by.y);
            // triplets.emplace_back(eid * 3 + 2, fid0 * 2 + 1, f0_by.z);


            // triplets.emplace_back(eid * 3, fid1 * 2, -f1_bx.x);
            // triplets.emplace_back(eid * 3 + 1, fid1 * 2, -f1_bx.y);
            // triplets.emplace_back(eid * 3 + 2, fid1 * 2, -f1_bx.z);
            // triplets.emplace_back(eid * 3, fid1 * 2 + 1, -f1_by.x);
            // triplets.emplace_back(eid * 3 + 1, fid1 * 2 + 1, -f1_by.y);
            // triplets.emplace_back(eid * 3 + 2, fid1 * 2 + 1, -f1_by.z);

            Vector2 rot = geom.transportVectorsAcrossHalfedge[he];
            triplets.emplace_back(eid * 2, fid0 * 2, rot.x);
            triplets.emplace_back(eid * 2, fid0 * 2 + 1, -rot.y);
            triplets.emplace_back(eid * 2 + 1, fid0 * 2, rot.y);
            triplets.emplace_back(eid * 2 + 1, fid0 * 2 + 1, rot.x);


            triplets.emplace_back(eid * 2, fid1 * 2, -1);
            triplets.emplace_back(eid * 2, fid1 * 2 + 1, 0);
            triplets.emplace_back(eid * 2 + 1, fid1 * 2, 0);
            triplets.emplace_back(eid * 2 + 1, fid1 * 2 + 1, -1);
        }
    }

    triplets.emplace_back(mesh.nEdges() * 2, 0, 1);
    triplets.emplace_back(mesh.nEdges() * 2 + 1, 1, 1);

    R_face_basis.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(mesh.nEdges() * 2 + 2);
    rhs[mesh.nEdges() * 2 + 0] = 1.0;
    rhs[mesh.nEdges() * 2 + 1] = 0.0;

    Eigen::SparseMatrix<double> AtA = R_face_basis.transpose() * R_face_basis;
    Eigen::VectorXd Atb = R_face_basis.transpose() * rhs;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(AtA);
    Eigen::VectorXd local_basis = solver.solve(Atb);


    if (solver.info() != Eigen::Success) {
        std::cerr << "[Error] LDLT decomposition also failed!" << std::endl;
    }

    FaceData<Eigen::Matrix2d> R(mesh);

    Eigen::SparseMatrix<double> A_R(mesh.nEdges() * 2 + 2, mesh.nFaces() * 2);
    Eigen::VectorXd b_x = Eigen::VectorXd::Zero(A_R.rows());
    Eigen::VectorXd b_y = Eigen::VectorXd::Zero(A_R.rows());
    triplets.clear();
    triplets.reserve(mesh.nEdges() * 8 + 4);

    int fixed_id = 0;
    for (Edge edge : mesh.edges()) {
        if (edge.isBoundary()) continue;

        Halfedge he = edge.halfedge();
        int i = he.twin().face().getIndex() * 2;
        int j = he.face().getIndex() * 2;

        int eq_row = edge.getIndex() * 2;

        Eigen::Matrix2d w_int = w[he].exp();

        triplets.emplace_back(eq_row, i, 1);
        triplets.emplace_back(eq_row + 1, i + 1, 1);
        triplets.emplace_back(eq_row, j, -w_int(0, 0));
        triplets.emplace_back(eq_row, j + 1, -w_int(0, 1));
        triplets.emplace_back(eq_row + 1, j, -w_int(1, 0));
        triplets.emplace_back(eq_row + 1, j + 1, -w_int(1, 1));
    }

    triplets.emplace_back(mesh.nEdges() * 2, 2 * fixed_id, 1);
    triplets.emplace_back(mesh.nEdges() * 2 + 1, 2 * fixed_id + 1, 1);
    A_R.setFromTriplets(triplets.begin(), triplets.end());
    A_R.makeCompressed();

    b_x(mesh.nEdges() * 2) = 1.0;
    b_x(mesh.nEdges() * 2 + 1) = 0.0;
    b_y(mesh.nEdges() * 2) = 0.0;
    b_y(mesh.nEdges() * 2 + 1) = 1.0;

    Eigen::SparseMatrix<double> ATA = A_R.transpose() * A_R;
    solver.compute(ATA);
    Eigen::VectorXd ATb_x = A_R.transpose() * b_x;
    Eigen::VectorXd r_x = solver.solve(ATb_x);
    Eigen::VectorXd ATb_y = A_R.transpose() * b_y;
    Eigen::VectorXd r_y = solver.solve(ATb_y);


    for (Face f : mesh.faces()) {
        Eigen::Vector2d rx(r_x(2 * f.getIndex()), r_x(2 * f.getIndex() + 1));
        Eigen::Vector2d ry(r_y(2 * f.getIndex()), r_y(2 * f.getIndex() + 1));

        // re-orthonormalize (to avoid drift from least squares)
        rx.normalize();
        ry = Eigen::Vector2d(-rx.y(), rx.x()); // enforce 90° rotation (counter-clockwise)

        Eigen::Matrix2d Rf;
        Rf.col(0) = rx;
        Rf.col(1) = ry;
        R[f] = Rf;
    }

    // try to propogate the field
    FaceData<bool> is_visted(mesh, false);
    is_visted[mesh.face(0)] = true;
    FaceData<Eigen::Vector2d> R_basis_x(mesh);
    R_basis_x[mesh.face(0)] = Eigen::Vector2d(1, 0);
    while (true) {
        bool is_changed = false;

        for (Face f : mesh.faces()) {
            for (Halfedge he : f.adjacentHalfedges()) {
                if (!he.edge().isBoundary()) {
                    Face fa = f;
                    Face fb = he.twin().face();
                    if (is_visted[fa] && !is_visted[fb]) {
                        is_changed = true;
                        is_visted[fb] = true;

                        Eigen::Vector2d vec_fa = R_basis_x[fa];
                        Eigen::Matrix2d rot = w[he];
                        rot = rot.exp();

                        Eigen::Vector2d vec_fb = rot * vec_fa;
                        Vector2 vec_fb_complex({vec_fa.x(), vec_fb.y()});
                        vec_fb_complex = geom.transportVectorsAcrossHalfedge[he] * vec_fb_complex;
                        vec_fb.x() = vec_fb_complex.x;
                        vec_fb.y() = vec_fb_complex.y;
                        R_basis_x[fb] = vec_fb;
                    }
                }
            }
        }

        if (!is_changed) {
            break;
        }
    }

    FaceData<Vector3> dp_euR_X(mesh), dp_euR_Y(mesh);
    Eigen::VectorXd dp_euR_X_vec(mesh.nFaces() * 3), dp_euR_Y_vec(mesh.nFaces() * 3);

    FaceData<Vector3> tangent_field(mesh);


    for (Face f : mesh.faces()) {
        double u_bar = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            u_bar += PD_u[he.tipVertex()];
        }
        u_bar /= 3.0;
        double eu = exp(u_bar);

        Eigen::Matrix2d RT = R[f].transpose();

        int fid = f.getIndex();

        // least square approach
        Eigen::Vector2d vec({local_basis[2 * fid + 0], local_basis[2 * fid + 1]});
        vec.normalized();
        vec = RT * vec;
        auto vecr = Eigen::Vector2d(-vec.y(), vec.x());


        // // propogation approach
        // Eigen::Vector2d vec = R_basis_x[f];
        // auto vecr = Eigen::Vector2d(-vec.y(), vec.x());


        dp_euR_X[f] =
            (geom.faceTangentBasis[f][0] * vec.x() + geom.faceTangentBasis[f][1] * vec.y())
                .normalize() *
            eu;
        dp_euR_Y[f] =
            (geom.faceTangentBasis[f][0] * vecr.x() + geom.faceTangentBasis[f][1] * vecr.y())
                .normalize() *
            eu;


        dp_euR_X_vec[fid * 3] = dp_euR_X[f].x;
        dp_euR_X_vec[fid * 3 + 1] = dp_euR_X[f].y;
        dp_euR_X_vec[fid * 3 + 2] = dp_euR_X[f].z;

        dp_euR_Y_vec[fid * 3] = dp_euR_Y[f].x;
        dp_euR_Y_vec[fid * 3 + 1] = dp_euR_Y[f].y;
        dp_euR_Y_vec[fid * 3 + 2] = dp_euR_Y[f].z;
    }


    this->dp_euRT_X = dp_euR_X;
    this->dp_euRT_Y = dp_euR_Y;


    geom.requirePolygonDivergenceMatrix();
    geom.requireCotanLaplacian();
    SparseMatrix<double> div_m = geom.polygonDivergenceMatrix;
    SparseMatrix<double> lap_m = geom.cotanLaplacian;


    fixed_id = 0;
    for (int i = 0; i < lap_m.cols(); i++) {
        lap_m.coeffRef(fixed_id, i) = 0;
    }
    lap_m.coeffRef(fixed_id, fixed_id) = 1.0;


    solver.compute(lap_m);
    if (solver.info() != Eigen::Success) {
        std::cerr << "[Error] LDLT decomposition also failed!" << std::endl;
    }

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

    Eigen::VectorXd PD_uv_u = solver.solve(rhs_u);
    Eigen::VectorXd PD_uv_v = solver.solve(rhs_v);

    VertexData<Vector2> ret_uv(mesh);

    for (Vertex v : mesh.vertices()) {
        int vid = v.getIndex();
        ret_uv[v].x = PD_uv_u[vid];
        ret_uv[v].y = PD_uv_v[vid];
    }

    this->dp_uv = ret_uv;

    // geom.requireHalfedgeVectorsInFace();

    // auto half_edge_projection = [&](Face face, Edge e, Eigen::Vector2d vec) {
    //     Halfedge he;
    //     for (Halfedge itr_he : face.adjacentHalfedges()) {
    //         if (itr_he.edge().getIndex() == e.getIndex()) {
    //             he = itr_he;
    //         }
    //     }
    //     Face f = he.face();
    //     Vector2 a = geom.halfedgeVectorsInFace[he].normalize();
    //     Vector2 b = Vector2({vec.x(), vec.y()}).normalize();
    //     return dot(a, b);
    // };

    // Eigen::SparseMatrix<double> uv_m(mesh.nEdges() * 2 + 2, mesh.nVertices() * 2);
    // Eigen::VectorXd euRTdp(mesh.nEdges() * 2 + 2);
    // triplets.clear();
    // for (Edge edge : mesh.edges()) {
    //     Halfedge he = edge.halfedge();
    //     int i = he.tailVertex().getIndex();
    //     int j = he.tipVertex().getIndex();
    //     int e = edge.getIndex();

    //     triplets.emplace_back(e * 2, i * 2, 1);
    //     triplets.emplace_back(e * 2 + 1, i * 2 + 1, 1);
    //     triplets.emplace_back(e * 2, j * 2, -1);
    //     triplets.emplace_back(e * 2 + 1, j * 2 + 1, -1);

    //     if (edge.isBoundary()) {
    //         Face face = he.face().isDead() || he.face() == Face() ? he.twin().face() : he.face();

    //         Eigen::Matrix2d RT = R[face].transpose();
    //         int fid = face.getIndex();
    //         Eigen::Vector2d vec({local_basis[2 * fid + 0], local_basis[2 * fid + 1]});
    //         vec.normalized();
    //         vec = RT * vec;
    //         auto vecr = Eigen::Vector2d(-vec.y(), vec.x());

    //         double euij = exp(0.5 * (PD_u[he.tipVertex()] + PD_u[he.tailVertex()]));

    //         Eigen::Vector2d Rijdp;

    //         Rijdp(0) = half_edge_projection(face, edge, vec);
    //         Rijdp(1) = half_edge_projection(face, edge, vecr);

    //         Eigen::Vector2d euRTdpij = euij * Rijdp;

    //         euRTdp[e * 2] = euRTdpij.x();
    //         euRTdp[e * 2 + 1] = euRTdpij.y();
    //     } else {
    //         Eigen::Matrix2d R_i = R[he.face()];
    //         Eigen::Matrix2d R_j = R[he.twin().face()];

    //         double euij = exp(0.5 * (PD_u[he.tipVertex()] + PD_u[he.tailVertex()]));

    //         // Fi
    //         int fidi = he.face().getIndex();
    //         Eigen::Vector2d veci({local_basis[2 * fidi + 0], local_basis[2 * fidi + 1]});
    //         veci.normalized();
    //         veci = R_i * veci;
    //         auto vecri = Eigen::Vector2d(-veci.y(), veci.x());
    //         double dpix = half_edge_projection(he.face(), edge, veci);
    //         double dpiy = half_edge_projection(he.face(), edge, vecri);
    //         // Fj
    //         int fidj = he.face().getIndex();
    //         Eigen::Vector2d vecj({local_basis[2 * fidj + 0], local_basis[2 * fidj + 1]});
    //         vecj.normalized();
    //         vecj = R_j * vecj;
    //         auto vecrj = Eigen::Vector2d(-vecj.y(), vecj.x());
    //         double dpjx = half_edge_projection(he.face(), edge, vecj);
    //         double dpjy = half_edge_projection(he.face(), edge, vecrj);


    //         Eigen::Vector2d Rijdp(0.5 * (dpix + dpjx), 0.5 * (dpiy + dpjy));

    //         Eigen::Vector2d euRTdpij = euij * Rijdp;

    //         euRTdp[e * 2] = euRTdpij.x();
    //         euRTdp[e * 2 + 1] = euRTdpij.y();
    //     }
    // }
    // triplets.emplace_back(mesh.nEdges() * 2, fixed_id * 2, 1);
    // triplets.emplace_back(mesh.nEdges() * 2 + 1, fixed_id * 2 + 1, 1);
    // euRTdp(mesh.nEdges() * 2) = 0;
    // euRTdp(mesh.nEdges() * 2 + 1) = 0;

    // uv_m.setFromTriplets(triplets.begin(), triplets.end());
    // uv_m.makeCompressed();

    // solver.compute(uv_m.transpose() * uv_m);
    // Eigen::VectorXd new_uv = solver.solve(uv_m.transpose() * euRTdp);

    // for (Vertex v : mesh.vertices()) {
    //     int vid = v.getIndex();
    //     ret_uv[v].x = new_uv[vid * 2];
    //     ret_uv[v].y = new_uv[vid * 2 + 1];
    // }

    // this->dp_uv = ret_uv;

    return ret_uv;
}

} // namespace parameterization_playground