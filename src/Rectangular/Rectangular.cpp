#include "Rectangular.hpp"

namespace parameterization_playground {


void RectangularParameterization::FlattenMesh(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom)
{
    geom.requireEdgeLengths();
    geom.requireCornerAngles();
    geom.requireFaceTangentBasis();
    geom.requireTransportVectorsAcrossHalfedge();


    SolveOptimizationProblem(mesh, geom);
}

void RectangularParameterization::hessObjective(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& u,
    VertexData<double>& v,
    FaceData<double>& theta,
    Eigen::SparseMatrix<double>& H)
{
    int nv = mesh.nVertices();
    int nf = mesh.nFaces();

    H.resize(2 * nv + nf, 2 * nv + nf);
    H.setZero();
    std::vector<Eigen::Triplet<double>> triplets;

    // default energy
    // phi = sum_i A_i (u_i^2 + v_i^2)
    for (Vertex vert : mesh.vertices()) {
        int vid = vert.getIndex();

        // u: uu uv
        triplets.emplace_back(vid, vid, 2 * geom.vertexDualAreas[vert]);
        triplets.emplace_back(vid, vid + nv, 2 * geom.vertexDualAreas[vert] * (u[vert] + v[vert]));

        // v: vu vv
        triplets.emplace_back(
            vid + nv,
            vid + nv,
            2 * geom.vertexDualAreas[vert] * (u[vert] + v[vert]));
        triplets.emplace_back(vid + nv, vid, 2 * geom.vertexDualAreas[vert]);
    }

    H.setFromTriplets(triplets.begin(), triplets.end());
}

void RectangularParameterization::gradObjective(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& u,
    VertexData<double>& v,
    FaceData<double>& theta,
    Eigen::VectorXd& g)
{
    int nv = mesh.nVertices();
    int nf = mesh.nFaces();

    g.setZero(2 * nv + nf);

    for (Vertex vert : mesh.vertices()) {
        int vid = vert.getIndex();
        g(vid) = 2 * geom.vertexDualAreas[vert] * u[vert];
        g(vid + nv) = 2 * geom.vertexDualAreas[vert] * v[vert];
    }
}

void BuildHessian(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& u,
    VertexData<double>& v,
    FaceData<double>& theta,
    EdgeData<double>& lambda,
    Eigen::SparseMatrix<double>& H)
{
    int nv = mesh.nVertices();
    int nf = mesh.nFaces();

    Eigen::SparseMatrix<double> H_v, H_theta;

    H_v.resize(nv, nf);
    H_theta.resize(nf, nf);

    auto signed_angle_2d = [](const Vector2& u, const Vector2& v) -> double {
        double dot_val = dot(u, v);
        double det_val = u.x * v.y - u.y * v.x;
        return std::atan2(det_val, dot_val);
    };

    auto hessian = [&](int i,
                       int j,
                       int k,
                       int ijk,
                       double lambda,
                       double eta,
                       double vi,
                       double vj,
                       double vk,
                       std::vector<double> cot_weights,
                       std::vector<double> corners,
                       std::vector<Eigen::Triplet<double>>& triplets_v_,
                       std::vector<Eigen::Triplet<double>>& triplets_theta_) {
        double cot_a_k_ij = cot_weights[0];
        double cot_a_i_jk = cot_weights[1];
        double cot_a_j_ki = cot_weights[2];
        double w = 0.5 * cot_a_k_ij;
        double lambda_ij = lambda;

        triplets_v_.emplace_back(
            i,
            ijk,
            2 * lambda_ij * w * (sin(2 * eta) - cos(2 * eta) * cot_a_j_ki));
        triplets_v_.emplace_back(
            j,
            ijk,
            2 * lambda_ij * w * (-sin(2 * eta) - cos(2 * eta) * cot_a_i_jk));
        triplets_v_.emplace_back(
            k,
            ijk,
            2 * lambda_ij * w * cos(2 - eta) * (cot_a_j_ki + cot_a_i_jk));

        triplets_theta_.emplace_back(ijk, ijk, -4 * lambda_ij * w * cos(2 * eta) * (vj - vi));
        triplets_theta_.emplace_back(
            ijk,
            ijk,
            -4 * lambda_ij * w * sin(2 * eta) * cot_a_j_ki * (vk - vi));
        triplets_theta_.emplace_back(
            ijk,
            ijk,
            -4 * lambda_ij * w * sin(2 * eta) * cot_a_i_jk * (vk - vj));
    };

    std::vector<Eigen::Triplet<double>> triplets_v, triplets_theta;

    for (Face face : mesh.faces()) {
        std::vector<Halfedge> halfedges;
        std::vector<Corner> corners;
        for (Halfedge he : face.adjacentHalfedges()) {
            halfedges.push_back(he);
            corners.push_back(he.corner());
        }
        // Vector2 X_basis = geom.faceTangentBasis[face][0];
        Vector2 dummy_basis({1, 0});
        Vector2 edge_relate_to_basis = geom.halfedgeVectorsInFace[halfedges[0]];
        double eta_ij = signed_angle_2d(edge_relate_to_basis, dummy_basis) + theta[face];
        double a_j = geom.cornerAngles[corners[1]];
        double eta_jk = eta_ij - (PI - a_j);
        double a_k = geom.cornerAngles[corners[2]];
        double eta_ki = eta_jk - (PI - a_k);

        hessian(
            halfedges[0].tailVertex().getIndex(),
            halfedges[1].tailVertex().getIndex(),
            halfedges[2].tailVertex().getIndex(),
            face.getIndex(),
            lambda[halfedges[0].edge()],
            eta_ij,
            v[halfedges[0].tailVertex()],
            v[halfedges[1].tailVertex()],
            v[halfedges[2].tailVertex()],
            {geom.halfedgeCotanWeights[halfedges[0]],
             geom.halfedgeCotanWeights[halfedges[1]],
             geom.halfedgeCotanWeights[halfedges[2]]},
            {geom.cornerAngles[corners[0]],
             geom.cornerAngles[corners[1]],
             geom.cornerAngles[corners[2]]},
            triplets_v,
            triplets_theta);


        hessian(
            halfedges[1].tailVertex().getIndex(),
            halfedges[2].tailVertex().getIndex(),
            halfedges[0].tailVertex().getIndex(),
            face.getIndex(),
            lambda[halfedges[1].edge()],
            eta_jk,
            v[halfedges[1].tailVertex()],
            v[halfedges[2].tailVertex()],
            v[halfedges[0].tailVertex()],
            {geom.halfedgeCotanWeights[halfedges[1]],
             geom.halfedgeCotanWeights[halfedges[2]],
             geom.halfedgeCotanWeights[halfedges[0]]},
            {geom.cornerAngles[corners[1]],
             geom.cornerAngles[corners[2]],
             geom.cornerAngles[corners[0]]},
            triplets_v,
            triplets_theta);


        hessian(
            halfedges[2].tailVertex().getIndex(),
            halfedges[0].tailVertex().getIndex(),
            halfedges[1].tailVertex().getIndex(),
            face.getIndex(),
            lambda[halfedges[2].edge()],
            eta_ij,
            v[halfedges[2].tailVertex()],
            v[halfedges[0].tailVertex()],
            v[halfedges[1].tailVertex()],
            {geom.halfedgeCotanWeights[halfedges[2]],
             geom.halfedgeCotanWeights[halfedges[0]],
             geom.halfedgeCotanWeights[halfedges[1]]},
            {geom.cornerAngles[corners[2]],
             geom.cornerAngles[corners[0]],
             geom.cornerAngles[corners[1]]},
            triplets_v,
            triplets_theta);
    }


    H_v.setFromTriplets(triplets_v.begin(), triplets_v.end());
    H_theta.setFromTriplets(triplets_v.begin(), triplets_v.end());

    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * nV + nF, 2 * nV + nF);
    H.resize(2 * nv + nf, 2 * nv + nf);
    H.setZero();

    H.block(nv, 2 * nv, nv, nf) = H_v.transpose();
    H.block(2 * nv, nv, nf, nv) = H_v;
    H.block(2 * nv, 2 * nv, nf, nf) = H_theta;
}


void RectangularParameterization::BuildJacobian(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& u,
    VertexData<double>& v,
    FaceData<double>& theta,
    EdgeData<double>& lambda,
    Eigen::SparseMatrix<double>& J)
{
    int ne = mesh.nEdges();
    int nv = mesh.nVertices();
    int nf = mesh.nFaces();

    Eigen::SparseMatrix<double> H_v, H_theta;
    std::vector<Eigen::Triplet<double>> triplets_u, triplets_v, triplets_theta;

    auto signed_angle_2d = [](const Vector2& u, const Vector2& v) -> double {
        double dot_val = dot(u, v);
        double det_val = u.x * v.y - u.y * v.x;
        return std::atan2(det_val, dot_val);
    };

    auto jacobian = [&](int i,
                        int j,
                        int k,
                        int ij,
                        int ijk,
                        double lambda,
                        double eta,
                        double vi,
                        double vj,
                        double vk,
                        std::vector<double> cot_weights,
                        std::vector<double> corners,
                        std::vector<Eigen::Triplet<double>>& triplets_u_,
                        std::vector<Eigen::Triplet<double>>& triplets_v_,
                        std::vector<Eigen::Triplet<double>>& triplets_theta_) {
        double w = 0.5 * cot_weights[0];
        double cot_a_k_ij = cot_weights[0];
        double cot_a_i_jk = cot_weights[1];
        double cot_a_j_ki = cot_weights[2];

        triplets_u_.emplace_back(ij, i, -w);
        triplets_u_.emplace_back(ij, j, w);

        triplets_v_.emplace_back(ij, i, -w * (cos(2 * eta) + sin(2 * eta) * cot_a_j_ki));
        triplets_v_.emplace_back(ij, j, -w * (cos(2 * eta) - sin(2 * eta) * cot_a_i_jk));
        triplets_v_.emplace_back(ij, k, -w * sin(2 * eta) * (cot_a_i_jk + cot_a_j_ki));

        triplets_theta_.emplace_back(ij, ijk, -2 * w * sin(2 * eta) * (vj - vi));
        triplets_theta_.emplace_back(ij, ijk, 2 * w * cos(2 * eta) * cot_a_j_ki * (vk - vi));
        triplets_theta_.emplace_back(ij, ijk, 2 * w * cos(2 * eta) * cot_a_i_jk * (vk - vj));
    };

    for (Edge edge : mesh.edges()) {
        if (edge.isBoundary()) {
            Halfedge he = edge.halfedge();
            if (he.face().isDead() || he.face() == Face()) {
                he = he.twin();
            }
            int fid = he.face().getIndex();
            triplets_theta.emplace_back(edge.getIndex(), fid, 1);
        } else {
            int fid = edge.halfedge().face().getIndex();
            int fid_twin = edge.halfedge().twin().face().getIndex();
            triplets_theta.emplace_back(edge.getIndex(), fid, 1);
            triplets_theta.emplace_back(edge.halfedge().twin().getIndex(), fid_twin, -1);
        }
    }

    for (Face face : mesh.faces()) {
        std::vector<Halfedge> halfedges;
        std::vector<Corner> corners;
        for (Halfedge he : face.adjacentHalfedges()) {
            halfedges.push_back(he);
            corners.push_back(he.corner());
        }
        // Vector2 X_basis = geom.faceTangentBasis[face][0];
        Vector2 dummy_basis({1, 0});
        Vector2 edge_relate_to_basis = geom.halfedgeVectorsInFace[halfedges[0]];
        double eta_ij = signed_angle_2d(edge_relate_to_basis, dummy_basis) + theta[face];
        double a_j = geom.cornerAngles[corners[1]];
        double eta_jk = eta_ij - (PI - a_j);
        double a_k = geom.cornerAngles[corners[2]];
        double eta_ki = eta_jk - (PI - a_k);

        jacobian(
            halfedges[0].tailVertex().getIndex(),
            halfedges[1].tailVertex().getIndex(),
            halfedges[2].tailVertex().getIndex(),
            halfedges[0].edge().getIndex(),
            face.getIndex(),
            lambda[halfedges[0].edge()],
            eta_ij,
            v[halfedges[0].tailVertex()],
            v[halfedges[1].tailVertex()],
            v[halfedges[2].tailVertex()],
            {geom.halfedgeCotanWeights[halfedges[0]],
             geom.halfedgeCotanWeights[halfedges[1]],
             geom.halfedgeCotanWeights[halfedges[2]]},
            {geom.cornerAngles[corners[0]],
             geom.cornerAngles[corners[1]],
             geom.cornerAngles[corners[2]]},
            triplets_u,
            triplets_v,
            triplets_theta);

        jacobian(
            halfedges[1].tailVertex().getIndex(),
            halfedges[2].tailVertex().getIndex(),
            halfedges[0].tailVertex().getIndex(),
            halfedges[1].edge().getIndex(),
            face.getIndex(),
            lambda[halfedges[1].edge()],
            eta_ij,
            v[halfedges[1].tailVertex()],
            v[halfedges[2].tailVertex()],
            v[halfedges[0].tailVertex()],
            {geom.halfedgeCotanWeights[halfedges[1]],
             geom.halfedgeCotanWeights[halfedges[2]],
             geom.halfedgeCotanWeights[halfedges[0]]},
            {geom.cornerAngles[corners[1]],
             geom.cornerAngles[corners[2]],
             geom.cornerAngles[corners[0]]},
            triplets_u,
            triplets_v,
            triplets_theta);

        jacobian(
            halfedges[2].tailVertex().getIndex(),
            halfedges[0].tailVertex().getIndex(),
            halfedges[1].tailVertex().getIndex(),
            halfedges[2].edge().getIndex(),
            face.getIndex(),
            lambda[halfedges[2].edge()],
            eta_ij,
            v[halfedges[2].tailVertex()],
            v[halfedges[0].tailVertex()],
            v[halfedges[1].tailVertex()],
            {geom.halfedgeCotanWeights[halfedges[2]],
             geom.halfedgeCotanWeights[halfedges[0]],
             geom.halfedgeCotanWeights[halfedges[1]]},
            {geom.cornerAngles[corners[2]],
             geom.cornerAngles[corners[0]],
             geom.cornerAngles[corners[1]]},
            triplets_u,
            triplets_v,
            triplets_theta);
    }

    J.resize(ne, nv * 2 + nf);
    J.setZero();

    Eigen::SparseMatrix<double> J_u, J_v, J_theta;
    J_u.setFromTriplets(triplets_u.begin(), triplets_u.end());
    J_v.setFromTriplets(triplets_v.begin(), triplets_v.end());
    J_theta.setFromTriplets(triplets_theta.begin(), triplets_theta.end());

    J.block(0, 0, ne, nv) = J_u;
    J.block(0, nv, ne, nv) = J_v;
    J.block(0, 2 * nv, ne, nf) = J_theta;
}


void RectangularParameterization::BuildSystem(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& u,
    VertexData<double>& v,
    FaceData<double>& theta,
    EdgeData<double>& lambda,
    Eigen::VectorXd& F)
{
    F.setZero(mesh.nEdges());

    auto residual = [&](int i,
                        int j,
                        int k,
                        double eta,
                        double ui,
                        double uj,
                        double uk,
                        double vi,
                        double vj,
                        double vk,
                        std::vector<double> cot_weights,
                        std::vector<double> corners) -> double {
        double rho = uj - ui;
        double cot_a_k_ij = cot_weights[0];
        double cot_a_i_jk = cot_weights[1];
        double cot_a_j_ki = cot_weights[2];

        rho += cos(2 * eta) * (vj - vi);
        rho += sin(2 * eta) * cot_a_j_ki * (vk - vi);
        rho += sin(2 * eta) * cot_a_i_jk * (vk - vj);
        rho += 0.5 * cot_a_k_ij;

        return rho;
    };

    auto signed_angle_2d = [](const Vector2& u, const Vector2& v) -> double {
        double dot_val = dot(u, v);
        double det_val = u.x * v.y - u.y * v.x;
        return std::atan2(det_val, dot_val);
    };

    for (Edge edge : mesh.edges()) {
        Halfedge he = edge.halfedge();
        int edge_id = edge.getIndex();
        if (edge.isBoundary()) {
            if (he.face() == Face() || he.face().isDead()) {
                he = he.twin();
            }
            // not sure if I need to do something in this block...
        } else {
            std::vector<Halfedge> halfedges({he, he.next(), he.next().next()});
            std::vector<Halfedge> halfedges_twin(
                {he.twin(), he.twin().next(), he.twin().next().next()});
            std::vector<Corner> corners(
                {halfedges[0].corner(), halfedges[1].corner(), halfedges[2].corner()});
            std::vector<Corner> corners_twin(
                {halfedges_twin[0].corner(),
                 halfedges_twin[1].corner(),
                 halfedges_twin[2].corner()});

            double angle_f0 = theta[he.face()];
            double angle_f1 = theta[he.twin().face()];
            Vector2 vec_f0({cos(angle_f0), sin(angle_f0)});
            Vector2 vec_f1({cos(angle_f1), sin(angle_f1)});

            Vector2 rot = geom.transportVectorsAcrossHalfedge[he];
            double rho_ij = angle(rot * vec_f0, vec_f1);

            Vertex vi = he.tailVertex();
            Vertex vj = he.tipVertex();
            Vertex vk = he.next().tipVertex();
            Vertex vl = he.twin().next().tipVertex();

            int i = vi.getIndex();
            int j = vj.getIndex();
            int k = vk.getIndex();
            int l = vl.getIndex();

            double eta_f0, eta_f1;
            eta_f0 = signed_angle_2d(geom.halfedgeVectorsInFace[halfedges[0]], vec_f0);
            eta_f1 = signed_angle_2d(geom.halfedgeVectorsInFace[halfedges_twin[0]], vec_f1);

            double rho_ij_k = residual(
                i,
                j,
                k,
                eta_f0,
                u[vi],
                u[vj],
                u[vk],
                v[vi],
                v[vj],
                v[vk],
                {geom.halfedgeCotanWeights[halfedges[0]],
                 geom.halfedgeCotanWeights[halfedges[1]],
                 geom.halfedgeCotanWeights[halfedges[2]]},
                {geom.cornerAngles[corners[0]],
                 geom.cornerAngles[corners[1]],
                 geom.cornerAngles[corners[2]]});

            double rho_ji_l = residual(
                j,
                i,
                l,
                eta_f1,
                u[vj],
                u[vi],
                u[vl],
                v[vj],
                v[vi],
                v[vl],
                {geom.halfedgeCotanWeights[halfedges_twin[0]],
                 geom.halfedgeCotanWeights[halfedges_twin[1]],
                 geom.halfedgeCotanWeights[halfedges_twin[2]]},
                {geom.cornerAngles[halfedges_twin[0]],
                 geom.cornerAngles[halfedges_twin[1]],
                 geom.cornerAngles[halfedges_twin[2]]});


            F(edge_id) = rho_ij + rho_ij_k + rho_ji_l;
        }
    }
}

void RectangularParameterization::SolveOptimizationProblem(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom)
{
    VertexData<double> u(mesh);
    VertexData<double> v(mesh);
    FaceData<double> theta(mesh);
    EdgeData<double> lambda(mesh);

    geom.requireVertexDualAreas();
    geom.requireFaceTangentBasis();
    geom.requireCornerAngles();
    geom.requireHalfedgeVectorsInFace();
    geom.requireHalfedgeCotanWeights();

    int nv = mesh.nVertices();
    int ne = mesh.nEdges();
    int nf = mesh.nFaces();

    while (true) {
        Eigen::SparseMatrix<double> H;
        hessObjective(mesh, geom, u, v, theta, H);
        Eigen::VectorXd g;
        gradObjective(mesh, geom, u, v, theta, g);
        Eigen::SparseMatrix<double> D;
        BuildHessian(mesh, geom, u, v, theta, lambda, D);

        Eigen::SparseMatrix<double> J;
        BuildJacobian(mesh, geom, u, v, theta, lambda, J);
        Eigen::VectorXd F;
        BuildSystem(mesh, geom, u, v, theta, lambda, F);

        // [u,v,theta,lambda]
        Eigen::VectorXd x(mesh.nVertices() * 2 + mesh.nEdges() + mesh.nFaces());
        Eigen::VectorXd l(mesh.nEdges());
        {
            for (Vertex vert : mesh.vertices()) {
                int vid = vert.getIndex();
                x(vid) = u[vert];
                x(vid + nv) = v[vert];
            }
            for (Face face : mesh.faces()) {
                int fid = face.getIndex();
                x(fid + 2 * nv) = theta[face];
            }
            for (Edge edge : mesh.edges()) {
                int eid = edge.getIndex();
                x(eid + 2 * nv + nf) = lambda[edge];
                l(eid) = lambda[edge];
            }
        }

        Eigen::SparseMatrix<double> A(2 * nv + ne + nf, 2 * nv + ne + nf);
        A.setZero();
        A.block(0, 0, 2 * nv + nf, 2 * nv + nf) = H + D;
        A.block(0, 2 * nv + nf, 2 * nv + nf, ne) = J.transpose();
        A.block(2 * nv + nf, 0, ne, 2 * nv + nf) = J;

        Eigen::VectorXd b(2 * mesh.nEdges());
        b << g + J.transpose() * l, F;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd y = solver.solve(b);

        double tau = 1;

        // line-search
        {
            double energy = (g + J.transpose() * l).norm() + F.norm();

            VertexData<double> u_tild(mesh);
            VertexData<double> v_tild(mesh);
            FaceData<double> theta_tild(mesh);
            EdgeData<double> lambda_tild(mesh);

            while (true) {
                Eigen::VectorXd x_tild = x + tau * y;
                {
                    for (Vertex vert : mesh.vertices()) {
                        int vid = vert.getIndex();
                        u_tild[vert] = x_tild(vid);
                        v_tild[vert] = x_tild(vid + nv);
                    }
                    for (Face face : mesh.faces()) {
                        int fid = face.getIndex();
                        theta_tild[face] = x_tild(fid + 2 * nv);
                    }
                    for (Edge edge : mesh.edges()) {
                        int eid = edge.getIndex();
                        lambda_tild[edge] = x_tild(eid + 2 * nv + nf);
                    }
                }

                Eigen::VectorXd g;
                gradObjective(mesh, geom, u, v, theta, g);
                Eigen::SparseMatrix<double> J;
                BuildJacobian(mesh, geom, u, v, theta, lambda, J);
                Eigen::VectorXd F;
                BuildSystem(mesh, geom, u, v, theta, lambda, F);

                double new_energy = (g + J.transpose() * l).norm() + F.norm();
                if (new_energy < (1 - 0.5 * tau) * energy) {
                    break;
                } else {
                    tau *= 0.9;
                }
            }
        }

        x = x + tau * y;

        {
            for (Vertex vert : mesh.vertices()) {
                int vid = vert.getIndex();
                u[vert] = x(vid);
                v[vert] = x(vid + nv);
            }
            for (Face face : mesh.faces()) {
                int fid = face.getIndex();
                theta[face] = x(fid + 2 * nv);
            }
            for (Edge edge : mesh.edges()) {
                int eid = edge.getIndex();
                lambda[edge] = x(eid + 2 * nv + nf);
            }
        }

        if ((g + J.transpose() * l).norm() + F.norm() < 10) {
            break;
        }
    }

    ret_u = u;
    ret_v = v;
    ret_theta = theta;
    ret_lambda = lambda;
}

} // namespace parameterization_playground