//
// Created by richard on 12/09/23.
//

#include "../include/qcprot512.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/qcprot512.cpp"

#include "hwy/foreach_target.h"

#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();

namespace distopia {
    namespace HWY_NAMESPACE {
        namespace hn = hwy::HWY_NAMESPACE;

        static double InnerProduct(double *A, double **coords1, double **coords2, int len, double *weight) {
            hn::ScalableTag<double> d;
            using VecT = decltype(hn::Zero(d));

            size_t nlanes = hn::Lanes(d);
            size_t remainder = len % nlanes;
            const double *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
            const double *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
            double G1 = 0.0, G2 = 0.0;

            A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

            auto vG1 = hn::Zero(d);
            auto vG2 = hn::Zero(d);
            VecT vA[9];

            if (weight != nullptr) {
                for (size_t i=0; i<len-remainder; i+=nlanes) {
                    auto vw = hn::LoadU(d, weight + i);

                    auto vx1 = hn::LoadU(d, fx1 + i);
                    vx1 = vx1 * vx1;
                    auto vy1 = hn::LoadU(d, fy1 + i);
                    vy1 = vy1 * vy1;
                    auto vz1 = hn::LoadU(d, fz1 + i);
                    vz1 = vz1 * vz1;

                    auto vx2 = hn::LoadU(d, fx2 + i);
                    vx2 = vx2 * vx2;
                    auto vy2 = hn::LoadU(d, fy2 + i);
                    vy2 = vy2 * vy2;
                    auto vz2 = hn::LoadU(d, fz2 + i);
                    vz2 = vz2 * vz2;

                    vG1 += vw * (vx1 + vy1 + vz1);
                    vG2 += vw * (vx2 + vy2 + vz2);

                    vA[0] += vx1 * vx2;
                    vA[1] += vx1 * vy2;
                    vA[2] += vx1 * vz2;

                    vA[3] += vy1 * vx2;
                    vA[4] += vy1 * vy2;
                    vA[5] += vy1 * vz2;

                    vA[6] += vz1 * vx2;
                    vA[7] += vz1 * vy2;
                    vA[8] += vz1 * vz2;
                }
                for (size_t i=len-remainder; i<len; ++i) {

                }
            }
            else {
                for (size_t i=0; i<len-remainder; i+=nlanes) {

                }
                for (size_t i=len-remainder; i<len; ++i) {

                }
            }

            return 0.5 * (G1 + G2);
        }

        static HWY_INLINE double Sum(const double *src, int len) {
            double s = 0.0;
            hn::ScalableTag<double> d;
            int nlanes = hn::Lanes(d);
            int remainder = len % nlanes;

            auto vs = hn::Zero(d);
            for (auto i=0; i<len-nlanes; i+=nlanes) {
                vs += hn::LoadU(d, src + i);
            }
            s += hn::ReduceSum(d, vs);
            for (auto i=len-remainder; i<len; ++i) {
                s += src[i];
            }

            return s;
        }

        static void CenterCoords(double **coords, int len, const double *weight, double wsum) {
            hn::ScalableTag<double> d;
            auto nlanes = hn::Lanes(d);
            size_t remainder = len % nlanes;
            double *x = coords[0];
            double *y = coords[1];
            double *z = coords[2];

            double xsum, ysum, zsum;
            auto vxsum = hn::Zero(d);
            auto vysum = hn::Zero(d);
            auto vzsum = hn::Zero(d);

            if (weight != nullptr) {
                for (size_t i=0; i<len - remainder; i += nlanes) {
                    auto lx = hn::LoadU(d, x + i);
                    auto ly = hn::LoadU(d, y + i);
                    auto lz = hn::LoadU(d, z + i);
                    auto w = hn::LoadU(d, weight + i);
                    lx *= w;
                    ly *= w;
                    lz *= w;
                    vxsum += lx;
                    vysum += ly;
                    vzsum += lz;
                }
                xsum = hn::ReduceSum(d, vxsum);
                ysum = hn::ReduceSum(d, vysum);
                zsum = hn::ReduceSum(d, vzsum);

                for (auto i=len-remainder; i<len; ++i) {
                    xsum += x[i] * weight[i];
                    ysum += y[i] * weight[i];
                    zsum += z[i] * weight[i];
                }
            }
            else {
                for (size_t i=0; i<len - remainder; i += nlanes) {
                    auto lx = hn::LoadU(d, x + i);
                    auto ly = hn::LoadU(d, y + i);
                    auto lz = hn::LoadU(d, z + i);
                    vxsum += lx;
                    vysum += ly;
                    vzsum += lz;
                }
                xsum = hn::ReduceSum(d, vxsum);
                ysum = hn::ReduceSum(d, vysum);
                zsum = hn::ReduceSum(d, vzsum);

                for (auto i=len-remainder; i<len; ++i) {
                    xsum += x[i];
                    ysum += y[i];
                    zsum += z[i];
                }
            }
            xsum /= wsum;
            ysum /= wsum;
            zsum /= wsum;

            vxsum = hn::Set(d, xsum);
            vysum = hn::Set(d, ysum);
            vzsum = hn::Set(d, zsum);
            for (size_t i=0; i<len-remainder; i += nlanes) {
                auto lx = hn::LoadU(d, x + i);
                auto ly = hn::LoadU(d, y + i);
                auto lz = hn::LoadU(d, z + i);

                lx -= vxsum;
                ly -= vysum;
                lz -= vzsum;

                hn::StoreU(lx, d, x + i);
                hn::StoreU(ly, d, y + i);
                hn::StoreU(lz, d, z + i);
            }
            for (auto i=len-remainder; i<len; ++i) {
                x[i] -= xsum;
                y[i] -= ysum;
                z[i] -= zsum;
            }
        }

        int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, double len, double minScore) {
            return 1;
        }


        double CalcRMSDRotationalMatrix(double **coords1, double **coords2, int len, double *rot, double *weight) {
            hn::ScalableTag<double> d;
            double A[9], rmsd, wsum, E0;

            wsum = (weight == nullptr) ? len : Sum(weight, len);

            CenterCoords(coords1, len, weight, wsum);
            CenterCoords(coords2, len, weight, wsum);

            E0 = InnerProduct(A, coords1, coords2, len, weight);

            FastCalcRMSDAndRotation(rot, A, &rmsd, E0, wsum, -1);

            return rmsd;
        }
    }
}

HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace distopia {
    HWY_EXPORT(CalcRMSDRotationalMatrix);

    HWY_DLLEXPORT double CalcRMSDRotationalMatrix(double **coords1, double **coords2, int len, double *rot, double *weight) {
        return HWY_DYNAMIC_DISPATCH(CalcRMSDRotationalMatrix)(coords1, coords2, len, rot, weight);
    }
}

#endif