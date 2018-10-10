/** Class realizes a 1D mesh
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 25.04.2013
 *
 *  @version 2
 *  GetLeftCell and GetRightCell added. Attention: Cell[0] and Cell[n] are returned as left and right neighbor of Cell[0] and Cell[n].
 *  They are therefore their own ghostcells. However one has to be aware that changing the values on those cells immediately changes the
 *  values of ghost cells and the other way round!
 *  @version 1
 *  class for 1D mesh consisting of 1D cells (first definition)
 */

#ifndef MESH1D_H
#define MESH1D_H

#include "Cell1D.h"
#include "Mesh.h"
#include <vector>

class Mesh1D : public Mesh
{
  private:
    /**
     * left and right boundaries in x direction
     */
    double _xLeft, _xRight;

    /**
     * stores dimension in x direction
     */
    unsigned int _dimX;

    /**
     * vector that contains the cells of the mesh
     */
    std::vector<Cell*> _cells;

    /**
     * total length of interval in x direction
     */
    double _lengthX;

    /**
     * cell width dx
     */
    double _dx;

  public:
    Mesh1D( double xLeft, double xRight, unsigned int dimX, PhysicalProblem* physProb, Settings* settings );    ///< constructor
    virtual ~Mesh1D();                                                                                          ///< destructor

    /**
     * get access to a specific cell
     * @param[in] i, cell id (range: 0 to _dimX-1)
     */
    Cell* GetCell( int i ) const;
    /**
     * Getter left cell. Left neighbor of cell 0 is cell 0!
     * @param[in] i cell index
     * @return cell left to cell i
     */
    virtual Cell* GetLeftCell( int i ) const;
    /**
     * Getter right cell. Right neighbor of cell n is cell n!
     * @param[in] i cell index
     * @return cell right to cell i
     */
    virtual Cell* GetRightCell( int i ) const;
    /**
     * Getter _dx
     * @return Cell width _dx
     */
    double GetDx() const;
    /**
     * Getter xLeft
     * @return left boundary position
     */
    double GetXLeft() const;
    /**
     * Getter _xRight
     * @return right boundary position
     */
    double GetXRight() const;
    /**
     * Getter _dimX
     * @return Dimension in x direction
     */
    int GetDimX() const;

  private:
    Mesh1D() {}
};

#endif    // MESH1D_H
