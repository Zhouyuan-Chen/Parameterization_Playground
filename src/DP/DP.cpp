#include "DP.hpp"

namespace parameterization_playground {
DP::DP() {}
DP::~DP() {}

VertexData<Vector2> DP::solve(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& constraint,
    ConstraintType type)
{
    VertexData<int> bvid(mesh);
    VertexData<int> ivid(mesh);

    auto isInterior = Vector<bool>(mesh.nVertices());

    // compute curvatures
    geom.requireVertexGaussianCurvatures();
    geom.requireVertexAngleSums();
    geom.requireEdgeLengths();
    VertexData<double> vertex_gaussian_curvatures(mesh, 0);
    int bid = 0, iid = 0;
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            vertex_gaussian_curvatures[v] = 1. * PI - geom.vertexAngleSums[v];
            isInterior[bid + iid] = false;
            bvid[v] = bid++;
        } else {
            vertex_gaussian_curvatures[v] = 2. * PI - geom.vertexAngleSums[v];
            isInterior[bid + iid] = true;
            ivid[v] = iid++;
        }
    }

    // aproximate Yamabe
    // Au = Omega - [0 , h]^T
    // h = k - k^*

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
    SparseMatrix<double> schur = ABB - ABI * cholesky_AII.solve(AIB);
    Eigen::SimplicialLDLT<SparseMatrix<double>> cholesky_schur;
    cholesky_schur.compute(schur);

    Vector<double> OmegaI(iid);
    for (Vertex v : mesh.vertices()) {
        if (!v.isBoundary()) {
            int id = ivid[v];
            OmegaI[id] = vertex_gaussian_curvatures[v];
        }
    }
    Vector<double> OmegaB(bid);
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            int id = bvid[v];
            OmegaB[id] = vertex_gaussian_curvatures[v];
        }
    }


    Eigen::VectorXd u_boundary(bid);
    Eigen::VectorXd u_interior(iid);
    Eigen::VectorXd h(bid);
    VertexData<double> target_curvature(mesh, 0);

    if (type == BND_SCALE) {
        for (Vertex v : mesh.vertices()) {
            if (v.isBoundary()) {
                int id = bvid[v];
                u_boundary[id] = constraint[v];
            }
        }
        u_interior = cholesky_AII.solve((OmegaI - AIB * u_boundary));
        h = OmegaB - ABI * u_interior - ABB * u_boundary;


        int cnt = 0;
        for (Vertex v : mesh.vertices()) {
            if (v.isBoundary()) {
                // target_curvature[v] = OmegaB[cnt] - h[cnt];
                target_curvature[v] = h[cnt++];
                // target_curvature[v] = 6.28 / (boundary_num * 0.1);
            }
        }

    } else if (type == BND_CURVATURE) {
        for (Vertex v : mesh.vertices()) {
            if (v.isBoundary()) {
                int id = bvid[v];
                h[id] = constraint[v];
            }
        }

        int interior_num = iid;
        int boundary_num = bid;

        Vector<double> schur_rhs = OmegaB - h - ABI * cholesky_AII.solve(OmegaI);
        u_boundary = -cholesky_schur.solve(schur_rhs);

        double mean_u_boundary = u_boundary.mean();
        for (int i = 0; i < u_boundary.size(); i++) {
            u_boundary[i] -= mean_u_boundary;
        }

        // geometrycentral::surface::BFF ref_bff(mesh, geom);
        // u_boundary = ref_bff.neumannToDirichlet(h);

        int cnt = 0;
        for (Vertex v : mesh.vertices()) {
            if (v.isBoundary()) {
                target_curvature[v] = h[bvid[v]];
            }
        }
    } else {
        std::cout << "error: wrong constraint type for BFF" << std::endl;
        exit(1);
    }

    // DenseMatrix<double> T(2, boundary_num);
    // SparseMatrix<double> N(boundary_num, boundary_num);
    // Vector<double> targetLength(mesh);

    Eigen::MatrixXd T(2, bid);
    Eigen::SparseMatrix<double> N(bid, bid);
    Eigen::VectorXd targetLength(bid);

    double phi = 0;
    for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
        int id = bvid[he.vertex()];
        N.coeffRef(id, id) = geom.edgeLengths[he.edge()];
        T(0, id) = cos(phi);
        T(1, id) = sin(phi);
        targetLength(id) =
            exp((u_boundary[bvid[he.tipVertex()]] + u_boundary(bvid[he.tailVertex()])) * 0.5) *
            geom.edgeLength(he.edge());
        phi += target_curvature[he.tipVertex()];
    }

    Vector<double> roundedLength =
        targetLength - N * T.transpose() * (T * N * T.transpose()).inverse() * (T * targetLength);


    VertexData<Vector2> uv(mesh);
    double u = 0;
    double v = 0;
    for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
        int id = bvid[he.vertex()];
        uv[he.tailVertex()].x = u;
        uv[he.tailVertex()].y = v;

        u += roundedLength(id) * T(0, id);
        v += roundedLength(id) * T(1, id);
    }

    Eigen::VectorXd U_boundary_u(bid);
    Eigen::VectorXd U_boundary_v(bid);
    Vector<double> position_u(bid);
    Vector<double> position_v(bid);
    {
        for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
            U_boundary_u(bvid[he.tailVertex()]) = -uv[he.tailVertex()].x;
            U_boundary_v(bvid[he.tailVertex()]) = uv[he.tailVertex()].y;
        }
    }

    Eigen::VectorXd U_interior_u = cholesky_AII.solve(-AIB * U_boundary_u);
    Eigen::VectorXd U_interior_v = cholesky_AII.solve(-AIB * U_boundary_v);

    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            int id = bvid[v];
            uv[v] = Vector2{U_boundary_u[id], U_boundary_v[id]};
        } else {
            int id = ivid[v];
            uv[v] = Vector2{U_interior_u[id], U_interior_v[id]};
        }
    }

    return uv;
}

} // namespace parameterization_playground