//=================================================================================================
/*!
//  \file blaze/math/traits/TSVecScalarDivExprTrait.h
//  \brief Header file for the TSVecScalarDivExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_TSVECSCALARDIVEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_TSVECSCALARDIVEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsSparseVector.h>
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
/*!\brief Auxiliary helper struct for the TSVecScalarDivExprTrait trait.
// \ingroup math_traits
*/
template< typename VT
        , typename ST
        , bool Condition >
struct TSVecScalarDivExprTraitHelper
{
 private:
   //**********************************************************************************************
   using ScalarType = If_< Or< IsFloatingPoint< UnderlyingBuiltin_<VT> >
                             , IsFloatingPoint< UnderlyingBuiltin_<ST> > >
                         , If_< And< IsComplex< UnderlyingNumeric_<VT> >
                                   , IsBuiltin<ST> >
                              , DivTrait_< UnderlyingBuiltin_<VT>, ST >
                              , DivTrait_< UnderlyingNumeric_<VT>, ST > >
                         , ST >;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   using Type = If_< IsInvertible<ScalarType>
                   , SVecScalarMultExpr<VT,ScalarType,true>
                   , SVecScalarDivExpr<VT,ScalarType,true> >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TSVecScalarDivExprTraitHelper class template.
// \ingroup math_traits
*/
template< typename VT
        , typename ST >
struct TSVecScalarDivExprTraitHelper<VT,ST,false>
{
 public:
   //**********************************************************************************************
   using Type = INVALID_TYPE;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluation of the expression type of a transpose sparse vector/scalar division.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a transpose
// sparse vector/scalar division. Given the transpose sparse vector type \a VT and the scalar
// type \a ST, the nested type \a Type corresponds to the resulting expression type. In case
// either \a VT is not a transpose sparse vector type or \a ST is not a scalar type, the
// resulting \a Type is set to \a INVALID_TYPE.
*/
template< typename VT        // Type of the left-hand side sparse vector
        , typename ST        // Type of the right-hand side scalar
        , typename = void >  // Restricting condition
struct TSVecScalarDivExprTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum : bool { condition = And< IsSparseVector<VT>, IsRowVector<VT>, IsNumeric<ST> >::value };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<VT>, IsVolatile<VT>, IsReference<VT>
                                , IsConst<ST>, IsVolatile<ST>, IsReference<ST> >
                            , TSVecScalarDivExprTrait< Decay_<VT>, Decay_<ST> >
                            , TSVecScalarDivExprTraitHelper<VT,ST,condition> >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the TSVecScalarDivExprTrait class template.
// \ingroup math_traits
//
// The TSVecScalarDivExprTrait_ alias declaration provides a convenient shortcut to access
// the nested \a Type of the TSVecScalarDivExprTrait class template. For instance, given
// the transpose sparse vector type \a VT and the scalar type \a ST the following two type
// definitions are identical:

   \code
   using Type1 = typename TSVecScalarDivExprTrait<VT,ST>::Type;
   using Type2 = TSVecScalarDivExprTrait_<VT,ST>;
   \endcode
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename ST >  // Type of the right-hand side scalar
using TSVecScalarDivExprTrait_ = typename TSVecScalarDivExprTrait<VT,ST>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
