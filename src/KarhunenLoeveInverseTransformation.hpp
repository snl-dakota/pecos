/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

#ifndef KARHUNEN_LOEVE_INVERSE_TRANSFORMATION_HPP
#define KARHUNEN_LOEVE_INVERSE_TRANSFORMATION_HPP

#include "DataTransformation.hpp"


namespace Pecos {


/// Class for KL data transformation.

/** The KarhunenLoeveInverseTransformation employs an KL decomposition
    to map from the frequency domain to the time domain. */

class KarhunenLoeveInverseTransformation: public DataTransformation
{
public:

  //
  //- Heading: Constructors and destructor
  //

  KarhunenLoeveInverseTransformation();  ///< constructor
  ~KarhunenLoeveInverseTransformation(); ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //


private:

  //
  //- Heading: Utility routines
  //

};


inline KarhunenLoeveInverseTransformation::KarhunenLoeveInverseTransformation():
  DataTransformation(BaseConstructor())
{ }


inline KarhunenLoeveInverseTransformation::~KarhunenLoeveInverseTransformation()
{ }

} // namespace Pecos

#endif
