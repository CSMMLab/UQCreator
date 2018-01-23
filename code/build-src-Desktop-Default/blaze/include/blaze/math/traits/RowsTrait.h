//=================================================================================================
/*!
//  \file blaze/math/traits/RowsTrait.h
//  \brief Header file for the rows trait
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_TRAITS_ROWSTRAIT_H_
#define _BLAZE_MATH_TRAITS_ROWSTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Types.h>
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
/*!\brief Base template for the RowsTrait class.
// \ingroup math_traits
//
// \section rowstrait_general General
//
// The RowsTrait class template offers the possibility to select the resulting data type when
// creating a view on a set of rows of a dense or sparse matrix. RowsTrait defines the nested
// type \a Type, which represents the resulting data type of the rows operation. In case the
// given data type is not a dense or sparse matrix type, the resulting data type \a Type is
// set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and reference modifiers
// are generally ignored.
//
//
// \section rowstrait_specializations Creating custom specializations
//
// Per default, RowsTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the RowsTrait template. The
// following example shows the according specialization for the DynamicMatrix class template:

   \code
   template< typename T1, bool SO, size_t... CRAs >
   struct RowsTrait< DynamicMatrix<T1,SO>, CRAs... >
   {
      using Type = DynamicMatrix<T1,false>;
   };
   \endcode

// \n \section rowstrait_examples Examples
//
// The following example demonstrates the use of the RowsTrait template, where depending on
// the given matrix type the resulting rows type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the rows type of a row-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,rowMajor>;
   using ResultType1 = typename blaze::RowsTrait<MatrixType1>::Type;

   // Definition of the rows type for the 1st and 3rd row of a column-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,4UL,3UL,columnMajor>;
   using ResultType2 = typename blaze::RowsTrait<MatrixType2,1UL,3UL>::Type;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CRAs >  // Compile time row arguments
struct RowsTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { using Type = INVALID_TYPE; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<MT>, IsVolatile<MT>, IsReference<MT> >
                            , RowsTrait< Decay_<MT>, CRAs... >
                            , Failure >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the RowsTrait type trait.
// \ingroup math_traits
//
// The RowsTrait_ alias declaration provides a convenient shortcut to access the nested
// \a Type of the RowsTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename RowsTrait<MT>::Type;
   using Type2 = RowsTrait_<MT>;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CRAs >  // Compile time row arguments
using RowsTrait_ = typename RowsTrait<MT,CRAs...>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
