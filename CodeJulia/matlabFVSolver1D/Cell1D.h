/** Class realizes a 1D grid cell
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 25.04.2013
 *
 *
 *  @version 1
 *  class for 1D grid cells (first definition)
 */

#ifndef CELL1D_H
#define CELL1D_H

#include "Cell.h"

class Cell1D : public Cell
{
  public:
    Cell1D( int dim );    ///< constructor
    virtual ~Cell1D();    ///< destructor
  private:
    Cell1D();    ///< invalid constructor
};

#endif    // CELL1D_H
