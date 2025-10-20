#include "BFF.hpp"


namespace parameterization_playground {
BFF::BFF() {}
BFF::~BFF() {}

inline double angle2pi(const Vector2& u, const Vector2& v)
{
    Vector2 uu = unit(u);
    Vector2 vv = unit(v);

    double dot_val = std::fmax(-1., std::fmin(1., dot(uu, vv)));
    double ang = std::acos(dot_val);

    double cross_val = uu.x * vv.y - uu.y * vv.x;
    if (cross_val < 0) {
        ang = 2.0 * M_PI - ang;
    }
    return ang;
}

VertexData<Vector2> BFF::solve(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& constraint,
    ConstraintType type)
{
    std::vector<Vector2> shape;
    return solve(mesh, geom, constraint, shape, type);
}

VertexData<Vector2>
BFF::solve(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, std::vector<Vector2>& shape)
{
    VertexData<double> constraint(mesh, 0);
    return solve(mesh, geom, constraint, shape, BND_FIT_SHAPE);
}

VertexData<Vector2> BFF::solve(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geom,
    VertexData<double>& constraint,
    std::vector<Vector2>& shape,
    ConstraintType type)
{
    VertexData<Vector2> uv(mesh);

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

    this->curvature = vertex_gaussian_curvatures;

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
                target_curvature[v] = h[cnt++];
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
    } else if (type == BND_FIT_SHAPE) {
        Eigen::VectorXd target_h(bid);
        Eigen::VectorXd target_h_pre(bid);
        target_h.setZero();
        target_h_pre.setZero();

        std::vector<double> shape_accumulate_len;
        double target_len_sum = 0.0;
        for (int i = 0; i < (int)shape.size(); i++) {
            shape_accumulate_len.push_back(target_len_sum);
            int j = i == (int)shape.size() - 1 ? 0 : i + 1;
            target_len_sum += (shape[i] - shape[j]).norm();
        }
        shape_accumulate_len.push_back(target_len_sum);

        for (int itr = 0; itr < 10; itr++) {
            std::cout << "itr=" << itr << std::endl;
            target_h_pre = target_h;
            // compute lengths
            double cur_len_sum = 0.0;
            if (itr == 0) {
                for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                    cur_len_sum += geom.edgeLength(he.edge());
                }
            } else {
                for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                    cur_len_sum += (uv[he.tipVertex()] - uv[he.tailVertex()]).norm();
                }
            }

            VertexData<double> cur_cumlen(mesh, 0);
            {
                double cur_accumulate_len = 0;
                for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                    cur_cumlen[he.tailVertex()] = cur_accumulate_len;
                    if (itr == 0) {
                        cur_accumulate_len += geom.edgeLength(he.edge());
                    } else {
                        cur_accumulate_len += (uv[he.tipVertex()] - uv[he.tailVertex()]).norm();
                    }
                }
            }

            // sample
            std::vector<Vector2> z_bnd(bid);

            for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                double t = cur_cumlen[he.tailVertex()] / cur_len_sum * target_len_sum;

                for (int i = 0; i < shape_accumulate_len.size(); i++) {
                    double low_bound = shape_accumulate_len[i];
                    double high_bound = (i == shape_accumulate_len.size() - 1)
                                            ? target_len_sum
                                            : shape_accumulate_len[i + 1];

                    if (t >= low_bound && t <= high_bound) {
                        double alpha = (t - low_bound) / (high_bound - low_bound);
                        int id = bvid[he.tailVertex()];

                        int next_i = (i + 1) % shape.size();
                        z_bnd[id] = (1 - alpha) * shape[i] + alpha * shape[next_i];
                        break;
                    }
                }
            }

            // recompute target curvatures
            double totalAngle = 0.0;
            VertexData<Halfedge> prev_edge(mesh);
            for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                prev_edge[he.tipVertex()] = he;
            }

            // compute angles
            for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                Vertex vi = prev_edge[he.tailVertex()].tailVertex();
                Vertex vj = he.tailVertex();
                Vertex vk = he.tipVertex();
                int id_i = bvid[vi];
                int id_j = bvid[vj];
                int id_k = bvid[vk];

                // safe: if any z_bnd entry is uninitialized (shouldn't be), guard
                Vector2 vec1 = z_bnd[id_k] - z_bnd[id_j];
                Vector2 vec2 = z_bnd[id_i] - z_bnd[id_j];
                double theta = angle2pi(vec1, vec2);
                // target_h_pre[id_j] = target_h[id_j]; // store previous
                target_h[id_j] = PI - theta;
                totalAngle += target_h[id_j];
                // std::cout << theta / M_PI << std::endl;
                // std::cout << z_bnd[id_j].x << " " << z_bnd[id_j].y << std::endl;
            }

            // orientation fix
            if (totalAngle < 0.0) {
                target_h *= -1.0;
                totalAngle *= -1.0;
            }

            // average with previous (0.5*(cur + prev))
            double interpolate_rate = 1.0;
            if (itr != 0)
                for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                    int id_j = bvid[he.tailVertex()];
                    target_h[id_j] = target_h[id_j] * (interpolate_rate) +
                                     (1 - interpolate_rate) * target_h_pre[id_j];
                }

            // closure correction: distribute (2π - totalAngle) weighted by dual lengths of sampled
            // z
            if (true) {
                double check_sum = 0;
                for (Halfedge he : mesh.boundaryLoop(0).adjacentHalfedges()) {
                    Vertex vi = prev_edge[he.tailVertex()].tailVertex();
                    Vertex vj = he.tailVertex();
                    Vertex vk = he.tipVertex();
                    int id_i = bvid[vi];
                    int id_j = bvid[vj];
                    int id_k = bvid[vk];

                    double len1 = (z_bnd[id_k] - z_bnd[id_j]).norm();
                    double len2 = (z_bnd[id_i] - z_bnd[id_j]).norm();
                    double ldual = 0.5 * (len1 + len2);

                    target_h[id_j] += (ldual / target_len_sum) * (2.0 * M_PI - totalAngle);

                    check_sum += target_h[id_j];
                }

                // std::cout << check_sum << std::endl;
                // if (target_h.sum() - 2 * M_PI > 1e-3) {
                //     double scale = (2.0 * M_PI) / check_sum;
                //     check_sum = 0;
                //     for (int j = 0; j < target_h.size(); j++) {
                //         target_h[j] *= scale;

                //         std::cout << target_h[j] / PI << std::endl;
                //     }
                //     // if (itr == 0) {
                //     //     double scale = (2.0 * M_PI) / check_sum;
                //     //     check_sum = 0;
                //     //     for (int j = 0; j < target_h.size(); j++) {
                //     //         target_h[j] *= scale;
                //     //         check_sum += target_h[j];
                //     //     }
                //     // } else {
                //     //     std::cout << "end at itr==" << itr << std::endl;
                //     //     return uv;
                //     // }
                // }
            }


            // build constraint and solve
            VertexData<double> constraint_h(mesh, 0.0);
            for (Vertex v : mesh.vertices()) {
                if (v.isBoundary()) {
                    int id = bvid[v];
                    constraint_h[v] = target_h[id];
                } else {
                    constraint_h[v] = 0;
                }
            }
            // std::cout << target_h.sum() << std::endl;
            uv = solve(mesh, geom, constraint_h, BND_CURVATURE);
        }

        for (Vertex v : mesh.vertices()) {
            uv[v].x *= -1;
        }

        return uv;

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

    result_u = VertexData<double>(mesh);
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) {
            int id = bvid[v];
            result_u[v] = u_boundary(id);
        } else {
            int id = ivid[v];
            result_u[v] = u_interior(id);
        }
    }

    return uv;
}


} // namespace parameterization_playground