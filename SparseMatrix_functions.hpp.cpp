//
// Created by Jack on 6/19/2018.
//

#include <TypeTraits.hpp>

//------------------------------------------------------------------------
//Compute matrix vector product y = A*x and return dot(x,y), where:
//
// A - input matrix
// x - input vector
// y - result vector
//

template<typename MatrixType, typename VectorType>
typename TypeTraits<typename VectorType::ScalarType>::magnitude_type matvec_and_dot(MatrixType& A, VectorType& x, VectorType& y)
{

    typedef typename TypeTraits<typename VectorType::ScalarType>::magnitude_type magnitude;
    typedef typename MatrixType::ScalarType ScalarType;
    typedef typename MatrixType::GlobalOrdinalType GlobalOrdinalType;
    typedef typename MatrixType::LocalOrdinalType LocalOrdinalType;

    int n = A.rows.size();
    const LocalOrdinalType* Arowoffsets = &A.row_offsets[0];
    const GlobalOrdinalType* Acols      = &A.packed_cols[0];
    const ScalarType* Acoefs            = &A.packed_coefs[0];
    const ScalarType* xcoefs = &x.coefs[0];
    ScalarType* ycoefs = &y.coefs[0];
    ScalarType beta = 0;
#if USES_QUIRE
    positX posit_result = 0;
        quireX result = 0;
        positX testCompare = -1.0;
        quireX tempQuire;
        quireX tempQuireMult;
        positX tempQuireMultholder;

        quireX sumHolder;
        positX sumHolderPosit;



        for(int row=0; row<n; ++row) {
            quireX sum = sw::unum::quire_mul(beta, ycoefs[row]);

            for(LocalOrdinalType i=Arowoffsets[row]; i<Arowoffsets[row+1]; ++i) {

                sum += sw::unum::quire_mul(Acoefs[i], xcoefs[Acols[i]]);

            }
            positX tempSum;
            tempSum.convert(sum.to_value());
            ycoefs[row] = tempSum;

            positX oldResult;
            positX newResult;

            oldResult.convert(result.to_value()); //set a temp Scalar = result

            result += sw::unum::quire_mul(xcoefs[row], tempSum); //update the value of result (a Quire)
            newResult.convert(result.to_value()); //snapshot the accumulated result in a value called tempResult (a posit)
            //to summarize, result has been snap-shotted before and after an accumulation has occured.

        }
        std::cout << "(matvec_and_dot) result (as quire):  = " << result << std::endl;
        posit_result.convert(result.to_value());
        std::cout << "(matvec_and_dot) result (as posit):  = " << posit_result << std::endl;


        return posit_result;
#else
    magnitude result = 0;

    for(int row=0; row<n; ++row) {
        ScalarType sum = beta*ycoefs[row];

        for(LocalOrdinalType i=Arowoffsets[row]; i<Arowoffsets[row+1]; ++i) {
            sum += Acoefs[i]*xcoefs[Acols[i]];
            if (row > 263 && row < 269) {
                std::cout <<"Row is " << row << ". Acoefs[" << row << "] = " << Acoefs[i] <<", xcoefs[Acols[i]] = " << xcoefs[Acols[i]] << std::endl;
            }
#if MINIFE_DEBUG
            //std::cout << "Acoefs[i] (" << Acoefs[i] << ") * xcoefs[Acols[i]] (" << xcoefs[Acols[i]] << ") = " << Acoefs[i]*xcoefs[Acols[i]]<< std::endl;
               //std::cout << "Same, but using quire_mul this time: " << tempQuire << std::endl;
#endif

        }
//            std::cout << "sum = " << sum << std::endl;

        ycoefs[row] = sum;
        ScalarType tempScalar = result;
        result += xcoefs[row]*sum;
        std::cout << result << " = (result + (xcoefs[row] * sum)) == (" << tempScalar << " + (" << xcoefs[row] << " * " << sum << ") " << std::endl;
    }

#ifdef HAVE_MPI
    magnitude local_dot = result, global_dot = 0;
  MPI_Datatype mpi_dtype = TypeTraits<magnitude>::mpi_type();
  MPI_Allreduce(&local_dot, &global_dot, 1, mpi_dtype, MPI_SUM, MPI_COMM_WORLD);
  return global_dot;
#else
    return result;
#endif
#endif
}