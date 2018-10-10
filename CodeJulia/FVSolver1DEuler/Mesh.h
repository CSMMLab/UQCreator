/** Class realizes a mesh
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 25.04.2013
 *
 *
 *  @version 2
 *  GetLeftCell and GetRightCell added. Attention: Cell[0] and Cell[n] are returned as left and right neighbor of Cell[0] and Cell[n].
 *  They are therefore their own ghostcells. However one has to be aware that changing the values on those cells immediately changes the
 *  values of ghost cells and the other way round!
 *  @version 1
 *  abstract class for mesh consisting of cells (first definition)
 */

#ifndef MESH_H
#define MESH_H

#include "Cell.h"
#include "PhysicalProblem.h"
#include "Settings.h"

class Mesh
{
  private:
  protected:
    int _nCells;    ///< number of cells
    Mesh() {}
    Settings* _settings;

  public:
    Mesh( int nCells, Settings* settings );      ///< constructor
    virtual ~Mesh();                             ///< destructor
    int GetNumberOfCells() const;                ///< getter for number of cells
    virtual Cell* GetCell( int i ) const = 0;    ///< getter for cell with index i
    /**
     * Getter left cell. Left neighbor of cell 0 is cell 0!
     * @param[in] i cell index
     * @return cell left to cell i
     */
    virtual Cell* GetLeftCell( int i ) const = 0;
    /**
     * Getter right cell. Right neighbor of cell n is cell n!
     * @param[in] i cell index
     * @return cell right to cell i
     */
    virtual Cell* GetRightCell( int i ) const = 0;
    virtual double GetDx() const              = 0;    ///< getter for cell width
    /**
     * Getter xLeft
     * @return left boundary position
     */
    virtual double GetXLeft() const = 0;
    /**
     * Getter _xRight
     * @return right boundary position
     */
    virtual double GetXRight() const = 0;
    /**
     * Getter _dimX
     * @return Dimension in x direction
     */
    virtual int GetDimX() const = 0;
};
#endif    // MESH_H
