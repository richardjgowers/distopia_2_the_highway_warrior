//
// Created by richard on 13/08/23.
//
#include <distopia.h>

#include <iostream>
#include <cstring>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/distopia.cpp"

#include "hwy/foreach_target.h"

#include "hwy/highway.h"
#include "hwy/print-inl.h"

#define DEBUG_DIST 0

HWY_BEFORE_NAMESPACE();

namespace roadwarrior {
    namespace HWY_NAMESPACE {
        namespace hn = hwy::HWY_NAMESPACE;

        template <class D, typename T = hn::TFromD<D>>
        struct NoBox {
            explicit NoBox(D) {};

            void MinimiseVectors(hn::VFromD<D> &vx,
                                 hn::VFromD<D> &vy,
                                 hn::VFromD<D> &vz) const {};
        };

        template <class D, typename T>
        struct OrthogonalBox {
            hn::VFromD<D> lx, ly, lz, ix, iy, iz;
            explicit OrthogonalBox(D d, const T *sbox) {
                this->lx = hn::Set(d, sbox[0]);
                this->ly = hn::Set(d, sbox[1]);
                this->lz = hn::Set(d, sbox[2]);
                this->ix = hn::Set(d, 1 / sbox[0]);
                this->iy = hn::Set(d, 1 / sbox[1]);
                this->iz = hn::Set(d, 1 / sbox[2]);
            };

            void MinimiseVectors(hn::VFromD<D> &vx,
                                 hn::VFromD<D> &vy,
                                 hn::VFromD<D> &vz) const {
                auto sx = ix * vx;
                auto dsx = sx - hn::Round(sx);
                auto sy = iy * vy;
                auto dsy = sy - hn::Round(sy);
                auto sz = iz * vz;
                auto dsz = sz - hn::Round(sz);
                vx = lx * dsx;
                vy = ly * dsy;
                vz = lz * dsz;
            };
        };
        template <class D, typename T = hn::TFromD<D>>
        struct TriclinicBox {
            hn::VFromD<D> xx, xy, yy, xz, yz, zz;

            explicit TriclinicBox(D d, const T *sbox) {
                this->xx = hn::Set(d, sbox[0]);
                this->xy = hn::Set(d, sbox[1]); this->yy = hn::Set(d, sbox[2]);
                this->xz = hn::Set(d, sbox[3]); this->yz = hn::Set(d, sbox[4]); this->zz = hn::Set(d, sbox[5]);
            };

            void MinimiseVectors(hn::VFromD<D> &vx,
                                 hn::VFromD<D> &vy,
                                 hn::VFromD<D> &vz) const {

            };
        };

        template <class V, typename T = hn::TFromV<V>, class B>
        HWY_INLINE V distance(const V &ax, const V &ay, const V &az,
                              const V &bx, const V &by, const V &bz,
                              const B &box) {
            auto dx = ax - bx;
            auto dy = ay - by;
            auto dz = az - bz;

            dx = dx * dx;
            dy = dy * dy;
            dz = dz * dz;

            box.MinimiseVectors(dx, dy, dz);

            auto acc = dx + dy;
            acc = acc + dz;

            auto out = hn::Sqrt(acc);

            return out;
        }

        template <typename T, typename B>
        void calc_bonds(const T* a, const T* b, int n, T* out, B &box) {
            const hn::ScalableTag<T> d;
            int nlanes = hn::Lanes(d);

            // temporary arrays used for problem sizes smaller than nlanes
            T a_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T b_sub[3 * HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            T out_sub[HWY_MAX_LANES_D(hn::ScalableTag<T>)];
            const T *a_src, *b_src;
            T *dst;

            if (HWY_UNLIKELY(n < nlanes)) {
                // input problem was too small to bother with
                memcpy(a_sub, a, 3 * n * sizeof(T));
                memcpy(b_sub, b, 3 * n * sizeof(T));

                a_src = &a_sub[0];
                b_src = &b_sub[0];
                dst = &out_sub[0];
            } else {
                a_src = &a[0];
                b_src = &b[0];
                dst = out;
            }

            auto a_x = hn::Undefined(d);
            auto a_y = hn::Undefined(d);
            auto a_z = hn::Undefined(d);
            auto b_x = hn::Undefined(d);
            auto b_y = hn::Undefined(d);
            auto b_z = hn::Undefined(d);

            for (int i=0; i<n; i += nlanes) {
                // to deal with end of loop, can't load starting further than a[n - nlanes] so clip i to that
                size_t p = HWY_MIN(i, n - nlanes);

                hn::LoadInterleaved3(d, a_src + 3 * p, a_x, a_y, a_z);
                hn::LoadInterleaved3(d, b_src + 3 * p, b_x, b_y, b_z);

                auto result = distance(a_x, a_y, a_z, b_x, b_y, b_z, box);
                #ifndef DEBUG_DIST
                hn::Print(d, "ax is: ", a_x, 0, nlanes);
                hn::Print(d, "ay is: ", a_y, 0, nlanes);
                hn::Print(d, "az is: ", a_z, 0, nlanes);
                std::cout << std::endl;
                hn::Print(d, "bx is: ", b_x, 0, nlanes);
                hn::Print(d, "by is: ", b_y, 0, nlanes);
                hn::Print(d, "bz is: ", b_z, 0, nlanes);
                hn::Print(d, "result is: ", result, 0, 16);
                #endif

                hn::StoreU(result, d, dst + p);
            }
            if (HWY_UNLIKELY(n < nlanes)) {
                memcpy(out, dst, n * sizeof(T));
            }
        }

        void calc_bonds_double(const double *a, const double *b, int n, double *out) {
            hn::ScalableTag<double> d;
            const NoBox vbox(d);
            calc_bonds(a, b, n, out, vbox);
        }
        void calc_bonds_single(const float *a, const float *b, int n, float *out) {
            hn::ScalableTag<float> d;
            const NoBox vbox(d);
            calc_bonds(a, b, n, out, vbox);
        }
        void calc_bonds_ortho_double(const double *a, const double *b, int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const OrthogonalBox vbox(d, box);
            calc_bonds(a, b, n, out, vbox);
        }
        void calc_bonds_ortho_single(const float *a, const float *b, int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const OrthogonalBox vbox(d, box);
            calc_bonds(a, b, n, out, vbox);
        }
        void calc_bonds_triclinic_double(const double *a, const double *b, int n, const double *box, double *out) {
            hn::ScalableTag<double> d;
            const TriclinicBox vbox(d, box);
            calc_bonds(a, b, n, out, vbox);
        }
        void calc_bonds_triclinic_single(const float *a, const float *b, int n, const float *box, float *out) {
            hn::ScalableTag<float> d;
            const TriclinicBox vbox(d, box);
            calc_bonds(a, b, n, out, vbox);
        }
    }
}

HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace roadwarrior {
    HWY_EXPORT(calc_bonds_double);
    HWY_EXPORT(calc_bonds_single);
    HWY_EXPORT(calc_bonds_ortho_double);
    HWY_EXPORT(calc_bonds_ortho_single);
    HWY_EXPORT(calc_bonds_triclinic_double);
    HWY_EXPORT(calc_bonds_triclinic_single);

    HWY_DLLEXPORT template <> void calc_bonds(const float* a, const float* b, int n, float* out) {
        // TODO: Could instead put small problem handling here, if n<16 manually dispatch to non-vector route
        //       Would benefit the codepath in all vector versions
        return HWY_DYNAMIC_DISPATCH(calc_bonds_single)(a, b, n, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds(const double* a, const double* b, int n, double* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_double)(a, b, n, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds_orthogonal(const float* a, const float* b, int n, const float *box, float* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_ortho_single)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds_orthogonal(const double* a, const double* b, int n, const double *box, double* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_ortho_double)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds_triclinic(const float* a, const float* b, int n, const float *box, float* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_triclinic_single)(a, b, n, box, out);
    }
    HWY_DLLEXPORT template <> void calc_bonds_triclinic(const double* a, const double* b, int n, const double *box, double* out) {
        return HWY_DYNAMIC_DISPATCH(calc_bonds_triclinic_double)(a, b, n, box, out);
    }
}

#endif