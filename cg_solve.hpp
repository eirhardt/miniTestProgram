#ifndef _cg_solve_hpp_
#define _cg_solve_hpp_

//@HEADER
// ************************************************************************
//
// MiniFE: Simple Finite Element Assembly and Solve
// Copyright (2006-2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// ************************************************************************
//@HEADER

#include <cmath>
#include <limits>

#include "Vector_functions.hpp"
#include "utils/mytimer.hpp"

#include "utils/outstream.hpp"
#include "posit/blas.hpp"
#include "jack_settings.hpp"
#include "SparseMatrix_functions.hpp"
#include "TypeTraits.hpp"

namespace miniFE {

    template<typename Scalar>
    void print_vec(const std::vector<Scalar>& vec, const std::string& name)
    {
        for(size_t i=0; i<vec.size(); ++i) {
            std::cout << name << "["<<i<<"]: " << vec[i] << std::endl;
        }
    }

    template<typename VectorType>
    bool breakdown(typename VectorType::ScalarType inner, const VectorType& v, const VectorType& w)
    {
        typedef typename VectorType::ScalarType Scalar;
        typedef typename TypeTraits<Scalar>::magnitude_type magnitude;

//This is code that was copied from Aztec, and originally written
//by my hero, Ray Tuminaro.
//
//Assuming that inner = <v,w> (inner product of v and w),
//v and w are considered orthogonal if
//  |inner| < 100 * ||v||_2 * ||w||_2 * epsilon

#if USE_POSITX
        magnitude vnorm = sw::unum::sqrt(dot(v,v));
        magnitude wnorm = sw::unum::sqrt(dot(w,w));
        return sw::unum::abs(inner) <= 100*vnorm*wnorm*std::numeric_limits<magnitude>::epsilon();
#else
        magnitude vnorm = std::sqrt(dot(v,v));
        magnitude wnorm = std::sqrt(dot(w,w));
        return std::abs(inner) <= 100*vnorm*wnorm*std::numeric_limits<magnitude>::epsilon();
#endif

    }

    template<typename OperatorType, typename VectorType, typename Matvec>
    void cg_solve(OperatorType& A, const VectorType& b, VectorType& x, Matvec matvec, typename OperatorType::LocalOrdinalType max_iter,
                  typename TypeTraits<typename OperatorType::ScalarType>::magnitude_type& tolerance, typename OperatorType::LocalOrdinalType& num_iters,
                  typename TypeTraits<typename OperatorType::ScalarType>::magnitude_type& normr, timer_type* my_cg_times)
    {
        typedef typename OperatorType::ScalarType ScalarType;
        typedef typename OperatorType::GlobalOrdinalType GlobalOrdinalType;
        typedef typename OperatorType::LocalOrdinalType LocalOrdinalType;
        typedef typename TypeTraits<ScalarType>::magnitude_type magnitude_type;

        timer_type t0 = 0, tWAXPY = 0, tDOT = 0, tMATVEC = 0, tMATVECDOT = 0;
        timer_type total_time = mytimer();

        int myproc = 0;
#ifdef HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#ifdef USE_MPI_PCONTROL
  MPI_Pcontrol(1);
#endif
#endif

        if (!A.has_local_indices) {
            std::cerr << "miniFE::cg_solve ERROR, A.has_local_indices is false, needs to be true. This probably means "
                      << "miniFE::make_local_matrix(A) was not called prior to calling miniFE::cg_solve."
                      << std::endl;
            return;
        }

        size_t nrows = A.rows.size();

        LocalOrdinalType ncols = A.num_cols;

#if MINIFE_DEBUG
        std::cout << "Number of rows = "<< nrows << std::endl;
        std::cout << "Number of columns = "<< ncols << std::endl;
#endif

        VectorType r(b.startIndex, nrows);
        VectorType p(0, ncols);
        VectorType Ap(b.startIndex, nrows);

        normr = 0;
        magnitude_type rtrans = 0;
        magnitude_type oldrtrans = 0;

        LocalOrdinalType print_freq = max_iter/10;
        if (print_freq>50) print_freq = 50;
        if (print_freq<1)  print_freq = 1;


        //More debugging:
#if MINIFE_DEBUG
//        std::cout << "Vector r: = "<< std::endl;
//        print_vec(r.coefs, "r");
//
//        std::cout << "Vector p: = "<< std::endl;
//        print_vec(p.coefs, "p");
//
//        std::cout << "Vector Ap: = "<< std::endl;
//        print_vec(Ap.coefs, "Ap");

#endif

        ScalarType one = 1.0;
        ScalarType zero = 0.0;
        //------------------------------------------------------------
//Compute the update of a vector with the sum of two scaled vectors where:
//
// w = alpha*x + beta*y
//
// x,y - input vectors
//
// alpha,beta - scalars applied to x and y respectively
//
// w - output vector
//
//    waxpby(typename VectorType::ScalarType alpha, const VectorType& x, typename VectorType::ScalarType beta, const VectorType& y, VectorType& w)
        TICK(); waxpby(one, x, zero, x, p); TOCK(tWAXPY);
#if MINIFE_DEBUG_VERBOSE_2
        std::cout << "vector p after waxpby(one, x, zero, x, p) = "<< std::endl;
        print_vec(p.coefs, "p");
#endif

        TICK(); matvec(A, p, Ap); TOCK(tMATVEC);

#if MINIFE_DEBUG_VERBOSE_2
        std::cout << "vector Ap after matvec(A, p, Ap) = "<< std::endl;
        print_vec(Ap.coefs, "Ap");
#endif

        TICK(); waxpby(one, b, -one, Ap, r); TOCK(tWAXPY);

#if MINIFE_DEBUG_VERBOSE_2
        std::cout << "vector r after waxpby(one, b, -one, Ap, r) = "<< std::endl;
        print_vec(r.coefs, "r");
#endif

        TICK(); rtrans = dot(r, r); TOCK(tDOT); //squaring the vector r's magnitude.

#if MINIFE_DEBUG
        std::cout << "rtrans value after dot(r, r) = "<< rtrans << std::endl;
#endif

#if USE_POSITX
        normr = sw::unum::sqrt(rtrans);
#else
        normr = std::sqrt(rtrans);
#endif

        if (myproc == 0) {
            std::cout << "Initial Residual = "<< normr << std::endl;
        }
#if USE_POSITX
        magnitude_type brkdown_tol = 0;
#else
        magnitude_type brkdown_tol = std::numeric_limits<magnitude_type>::epsilon();
#endif

#ifdef MINIFE_DEBUG
        std::ostream& os = outstream();
        os << "brkdown_tol = " << brkdown_tol << std::endl;
#endif

        //MASTER LOOP WHERE HEAVY LIFTING IS DONE
        for(LocalOrdinalType k=1; k <= max_iter && normr > tolerance; ++k) {
            if (k == 1) {
                TICK(); waxpby(one, r, zero, r, p); TOCK(tWAXPY);
#if MINIFE_DEBUG
                std::cout << "[First run inside for-loop]: vector p after waxpby(one, r, zero, r, p) = "<< std::endl;
                print_vec(p.coefs, "p");
#endif
            }
            else {
                oldrtrans = rtrans;
                TICK(); rtrans = dot(r, r); TOCK(tDOT);
#if MINIFE_DEBUG
                std::cout << "oldrtrans == "<< oldrtrans << std::endl;
                std::cout << "rtrans == dot(r, r) == "<< rtrans << std::endl;
#endif
//      TICK(); rtrans = sw::blas::positX_fused_dot(nrows, r.coefs, 1, r.coefs, 1); TOCK(tDOT); //*dot product alteration - using Stillwater's fused_dot instead of miniFE as test.
                magnitude_type beta = rtrans/oldrtrans;
#if MINIFE_DEBUG
                std::cout << "magnitude_type beta = rtrans/oldrtrans: = "<< beta << std::endl;
#endif

                TICK(); waxpby(one, r, beta, p, p); TOCK(tWAXPY);

            }

#if USE_POSITX
            normr = sw::unum::sqrt(rtrans);
#else
            normr = std::sqrt(rtrans);
#endif

#if MINIFE_DEBUG
            std::cout << "residual (sqrt of rtrans) = "<< normr << std::endl;
#endif

            if (myproc == 0 && (k%print_freq==0 || k==max_iter)) {
                std::cout << "Iteration = "<<k<<"   Residual = "<<normr<<std::endl;
#if MINIFE_DEBUG
                std::cout << "Max iterations? " << max_iter << ", residual? " << normr << ", tolerance? " << tolerance << std::endl; //debugging
#endif
            }

            magnitude_type alpha = 0;
            magnitude_type p_ap_dot = 0;

#ifdef MINIFE_FUSED
            TICK(); p_ap_dot = matvec_and_dot(A, p, Ap); TOCK(tMATVECDOT);
#if MINIFE_DEBUG
            std::cout << "p_ap_dot after matvec_and_dot = "<< p_ap_dot << std::endl;
#endif
#else
            TICK(); matvec(A, p, Ap); TOCK(tMATVEC);
#if MINIFE_DEBUG_VERBOSE_2
            std::cout << "vector Ap after matvec(A, p, Ap) = "<< std::endl;
            print_vec(Ap.coefs, "Ap");
#endif

            TICK(); p_ap_dot = dot(Ap, p); TOCK(tDOT);
#if MINIFE_DEBUG
            std::cout << "p_ap_dot after dot(Ap, p) = "<< p_ap_dot << std::endl;
#endif

#endif

#if MINIFE_DEBUG
            os << "iter " << k << ", p_ap_dot = " << p_ap_dot;
            os.flush();
#endif
            if (p_ap_dot < brkdown_tol) {
                if (p_ap_dot < 0 || breakdown(p_ap_dot, Ap, p)) {
                    std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"<<std::endl;
#ifdef MINIFE_DEBUG
                    os << "ERROR, numerical breakdown!"<<std::endl;
#endif
                    //update the timers before jumping out.
                    my_cg_times[WAXPY] = tWAXPY;
                    my_cg_times[DOT] = tDOT;
                    my_cg_times[MATVEC] = tMATVEC;
                    my_cg_times[TOTAL] = mytimer() - total_time;
                    return;

                } else {
                    brkdown_tol = 0.1 * p_ap_dot;
                }
            }
            alpha = rtrans/p_ap_dot;
#if MINIFE_DEBUG
            os << ", rtrans = " << rtrans << ", alpha = " << alpha << std::endl;
#endif

#ifdef MINIFE_FUSED
            TICK();
            fused_waxpby(one, x, alpha, p, x, one, r, -alpha, Ap, r);
            TOCK(tWAXPY);
#else
            TICK(); waxpby(one, x, alpha, p, x);
            waxpby(one, r, -alpha, Ap, r); TOCK(tWAXPY);
#endif

            num_iters = k;
        }

#ifdef HAVE_MPI
        #ifdef USE_MPI_PCONTROL
  MPI_Pcontrol(0);
#endif
#endif

        my_cg_times[WAXPY] = tWAXPY;
        my_cg_times[DOT] = tDOT;
        my_cg_times[MATVEC] = tMATVEC;
        my_cg_times[MATVECDOT] = tMATVECDOT;
        my_cg_times[TOTAL] = mytimer() - total_time;
    }

}//namespace miniFE

#endif

