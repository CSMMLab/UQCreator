//=================================================================================================
/*!
//  \file blaze/math/traits/TDVecTDVecMapExprTrait.h
//  \brief Header file for the TDVecTDVecMapExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_TDVECTDVECMAPEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_TDVECTDVECMAPEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/typetraits/Decay.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>
#include <blaze/util/typetraits/IsVolatile.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the expression type of a transpose dense vector/transpose dense vector
//        map operation.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a transpose
// dense vector/transpose dense vector map operation. Given the two transpose dense vector types
// \a VT1 and \a VT2 and the custom operation type \a OP, the nested type \a Type corresponds to
// the resulting expression type. In case either \a VT1 or \a VT2 is not a transpose dense vector
// type, the resulting \a Type is set to \a INVALID_TYPE.
*/
template< typename VT1       // Type of the left-hand side transpose dense vector
        , typename VT2       // Type of the right-hand side transpose dense vector
        , typename OP        // Type of the custom operation
        , typename = void >  // Restricting condition
struct TDVecTDVecMapExprTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Tmp = If< And< IsDenseVector<VT1>, IsRowVector<VT1>
                      , IsDenseVector<VT2>, IsRowVector<VT2> >
                 , DVecDVecMapExpr<VT1,VT2,OP,true>
                 , INVALID_TYPE >;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<VT1>, IsVolatile<VT1>, IsReference<VT1>
                                , IsConst<VT2>, IsVolatile<VT2>, IsReference<VT2> >
                            , TDVecTDVecMapExprTrait< Decay_<VT1>, Decay_<VT2>, OP >
                            , Tmp >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the TDVecTDVecMapExprTrait class template.
// \ingroup math_traits
//
// The TDVecTDVecMapExprTrait_ alias declaration provides a convenient shortcut to access
// the nested \a Type of the TDVecTDVecMapExprTrait class template. For instance, given the
// transpose dense vector types \a VT1 and \a VT2 and the custom operation type \a OP the
// following two type definitions are identical:

   \code
   using Type1 = typename TDVecTDVecMapExprTrait<VT1,VT2,OP>::Type;
   using Type2 = TDVecTDVecMapExprTrait_<VT1,VT2,OP>;
   \endcode
*/
template< typename VT1   // Type of the left-hand side transpose dense vector
        , typename VT2   // Type of the right-hand side transpose dense vector
        , typename OP >  // Type of the custom operation
using TDVecTDVecMapExprTrait_ = typename TDVecTDVecMapExprTrait<VT1,VT2,OP>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
