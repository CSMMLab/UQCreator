//=================================================================================================
/*!
//  \file blaze/math/traits/CrossTrait.h
//  \brief Header file for the cross product trait
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

#ifndef _BLAZE_MATH_TRAITS_CROSSTRAIT_H_
#define _BLAZE_MATH_TRAITS_CROSSTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsCustom.h>
#include <blaze/math/typetraits/IsInitializer.h>
#include <blaze/math/typetraits/IsView.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Not.h>
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
/*!\brief Base template for the CrossTrait class.
// \ingroup math_traits
//
// \section crosstrait_general General
//
// The CrossTrait class template offers the possibility to select the resulting data type of
// a generic cross product operation between the two given types \a T1 and \a T2. CrossTrait
// defines the nested type \a Type, which represents the resulting data type of the cross
// product. In case \a T1 and \a T2 cannot be combined in a cross product, the resulting data
// type \a Type is set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and
// reference modifiers are generally ignored.
//
//
// \n \section crosstrait_specializations Creating custom specializations
//
// Per default, CrossTrait supports all vector types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the CrossTrait template. The
// following example shows the according specialization for the cross product between two static
// column vectors:

   \code
   template< typename T1, typename T2 >
   struct CrossTrait< StaticVector<T1,3UL,false>, StaticVector<T2,3UL,false> >
   {
      using Type = StaticVector< typename SubTrait< typename MultTrait<T1,T2>::Type
                                                  , typename MultTrait<T1,T2>::Type >::Type, 3UL, false >;
   };
   \endcode

// \n \section crosstrait_examples Examples
//
// The following example demonstrates the use of the CrossTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >  // The two generic types
   typename CrossTrait<T1,T2>::Type      // The resulting generic return type
   cross( T1 t1, T2 t2 )                 //
   {                                     // The function 'cross' returns the cross
      return t1 % t2;                    // product of the two given values
   }                                     //
   \endcode
*/
template< typename T1        // Type of the left-hand side operand
        , typename T2        // Type of the right-hand side operand
        , typename = void >  // Restricting condition
struct CrossTrait
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
   using Type = typename If_< Or< IsConst<T1>, IsVolatile<T1>, IsReference<T1>
                                , IsConst<T2>, IsVolatile<T2>, IsReference<T2> >
                            , CrossTrait< Decay_<T1>, Decay_<T2> >
                            , Failure >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the CrossTrait class template for the left operand being a custom or
//        view type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct CrossTrait< T1, T2
                 , EnableIf_< And< Or< IsCustom<T1>, IsInitializer<T1>, IsView<T1> >
                                 , Not< Or< IsCustom<T2>, IsInitializer<T2>, IsView<T2> > > > > >
{
 public:
   //**********************************************************************************************
   using Type = typename CrossTrait< typename T1::ResultType, T2 >::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the CrossTrait class template for the right operand being a custom or
//        view type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct CrossTrait< T1, T2
                 , EnableIf_< And< Not< Or< IsCustom<T1>, IsInitializer<T1>, IsView<T1> > >
                                 , Or< IsCustom<T2>, IsInitializer<T2>, IsView<T2> > > > >
{
 public:
   //**********************************************************************************************
   using Type = typename CrossTrait< T1, typename T2::ResultType >::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the CrossTrait class template for the both operands being custom or
//        view types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct CrossTrait< T1, T2
                 , EnableIf_< And< Or< IsCustom<T1>, IsInitializer<T1>, IsView<T1> >
                                 , Or< IsCustom<T2>, IsInitializer<T2>, IsView<T2> > > > >
{
 public:
   //**********************************************************************************************
   using Type = typename CrossTrait< typename T1::ResultType, typename T2::ResultType >::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the CrossTrait class template.
// \ingroup math_traits
//
// The CrossTrait_ alias declaration provides a convenient shortcut to access the nested \a Type
// of the CrossTrait class template. For instance, given the types \a T1 and \a T2 the following
// two type definitions are identical:

   \code
   using Type1 = typename CrossTrait<T1,T2>::Type;
   using Type2 = CrossTrait_<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
using CrossTrait_ = typename CrossTrait<T1,T2>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
