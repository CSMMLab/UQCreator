//=================================================================================================
/*!
//  \file blaze/math/traits/TSMatDeclUppExprTrait.h
//  \brief Header file for the TSMatDeclUppExprTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_TSMATDECLUPPEXPRTRAIT_H_
#define _BLAZE_MATH_TRAITS_TSMATDECLUPPEXPRTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUpper.h>
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
/*!\brief Evaluation of the expression type of a sparse matrix declupp operation.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting expression type of a sparse
// matrix declupp operation. Given the column-major sparse matrix type \a MT, the nested type
// \a Type corresponds to the resulting expression type. In case \a MT is not a column-major
// sparse matrix type, the resulting \a Type is set to \a INVALID_TYPE.
*/
template< typename MT        // Type of the sparse matrix
        , typename = void >  // Restricting condition
struct TSMatDeclUppExprTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Tmp = If< And< IsSparseMatrix<MT>, IsColumnMajorMatrix<MT> >
                 , If_< IsUniLower<MT>
                      , IdentityMatrix< ElementType_<MT>, true >
                      , If_< IsUpper<MT>
                           , const MT&
                           , SMatDeclUppExpr<MT,true> > >
                 , INVALID_TYPE >;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<MT>, IsVolatile<MT>, IsReference<MT> >
                            , TSMatDeclUppExprTrait< Decay_<MT> >
                            , Tmp >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the TSMatDeclUppExprTrait class template.
// \ingroup math_traits
//
// The TSMatDeclUppExprTrait_ alias declaration provides a convenient shortcut to access
// the nested \a Type of the TSMatDeclUppExprTrait class template. For instance, given the
// column-major sparse matrix type \a MT the following two type definitions are identical:

   \code
   using Type1 = typename TSMatDeclUppExprTrait<MT>::Type;
   using Type2 = TSMatDeclUppExprTrait_<MT>;
   \endcode
*/
template< typename MT >  // Type of the sparse matrix
using TSMatDeclUppExprTrait_ = typename TSMatDeclUppExprTrait<MT>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
