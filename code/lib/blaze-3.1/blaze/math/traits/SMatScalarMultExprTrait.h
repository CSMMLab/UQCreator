//=================================================================================================
/*!
//  \file blaze/math/traits/SMatScalarMultExprTrait.h
//  \brief Header file for the SMatScalarMultExprTrait class template
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_TRAITS_SMATSCALARMULTEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_SMATSCALARMULTEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/Decay.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the SMatScalarMultExprTrait trait.
// \ingroup math_traits
*/
template< typename MT
        , typename ST
        , bool Condition >
struct SMatScalarMultExprTraitHelper
{
 private:
   //**********************************************************************************************
   using ScalarType = If_< And< IsComplex< UnderlyingNumeric_<MT> >, IsBuiltin<ST> >
                         , MultTrait_< UnderlyingBuiltin_<MT>, ST >
                         , MultTrait_< UnderlyingNumeric_<MT>, ST > >;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   using Type = SMatScalarMultExpr<MT,ScalarType,false>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SMatScalarMultExprTraitHelper class template.
// \ingroup math_traits
*/
template< typename MT
        , typename ST >
struct SMatScalarMultExprTraitHelper<MT,ST,false>
{
 public:
   //**********************************************************************************************
   using Type = INVALID_TYPE;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluation of the expression type of a sparse matrix/scalar multiplication.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a sparse
// matrix/scalar multiplication. Given the row-major sparse matrix type \a MT and the scalar
// type \a ST, the nested type \a Type corresponds to the resulting expression type. In case
// either \a MT is not a row-major sparse matrix type or \a ST is not a scalar type, the
// resulting \a Type is set to \a INVALID_TYPE.
*/
template< typename MT        // Type of the left-hand side sparse matrix
        , typename ST        // Type of the right-hand side scalar
        , typename = void >  // Restricting condition
struct SMatScalarMultExprTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum : bool { condition = And< IsSparseMatrix<MT>, IsRowMajorMatrix<MT>, IsNumeric<ST> >::value };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<MT>, IsVolatile<MT>, IsReference<MT>
                                , IsConst<ST>, IsVolatile<ST>, IsReference<ST> >
                            , SMatScalarMultExprTrait< Decay_<MT>, Decay_<ST> >
                            , SMatScalarMultExprTraitHelper<MT,ST,condition> >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the SMatScalarMultExprTrait class template.
// \ingroup math_traits
//
// The SMatScalarMultExprTrait_ alias declaration provides a convenient shortcut to access the
// nested \a Type of the SMatScalarMultExprTrait class template. For instance, given the row-major
// sparse matrix type \a MT and the scalar type \a ST the following two type definitions are
// identical:

   \code
   using Type1 = typename SMatScalarMultExprTrait<MT,ST>::Type;
   using Type2 = SMatScalarMultExprTrait_<MT,ST>;
   \endcode
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , typename ST >  // Type of the right-hand side scalar
using SMatScalarMultExprTrait_ = typename SMatScalarMultExprTrait<MT,ST>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
