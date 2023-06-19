
///////////////////////////////////////////////////////////////////////////////////
//    ReferenceImplementations.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef ReferenceImplementations_hh
#define ReferenceImplementations_hh 1

#include "Operator.hh"


namespace ReferenceImplementations
{


  void comm110ss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void comm220ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm111ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm121ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm221ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm122ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm222_pp_hhss( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_phss( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_pp_hh_221ss( const Operator& X, const Operator& Y, Operator& Z) ;

  void comm330ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm331ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm231ss( const Operator& X, const Operator& Y, Operator& Z ) ;           

  void comm132ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm232ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z ) ;      
  
  void comm133ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm223ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm233_pp_hhss( const Operator& X, const Operator& Y, Operator& Z ) ;     
  void comm233_phss( const Operator& X, const Operator& Y, Operator& Z ) ;        

  void comm333_ppp_hhhss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z ) ;          


  // scalar-tensor commutators
  void comm222_phst( const Operator& X, const Operator& Y, Operator& Z) ;


}// namespace ReferenceImplementations

#endif
