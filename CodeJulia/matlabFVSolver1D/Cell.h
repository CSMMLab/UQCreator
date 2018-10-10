/** Class realizes a general grid cell
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 25.04.2013
 *
 *
 *  @version 1
 *  abstract class for grid cells (first definition)
 */

#ifndef cell_inc
#define cell_inc

#include "Vector.h"

class Cell
{
  protected:
    Vector* _state;    ///< stores state of current cell
    /**
     * @li _north cell above current cell
     * @li _east cell right to current cell
     * @li _south cell beneath current cell
     * @li _west cell left to current cell
     */
    // Cell* _north, _east, _south, _west;
  public:
    Cell( int dim );    ///< constructor
    virtual ~Cell();    ///< destructor

    /**  getter for vector '_state'
     *  @return state vector of cell averages
     */
    Vector* GetState();

    /**  setter for vector '_state'
     *  @param[in] vec, new values for vector '_state'
     */
    void SetState( Vector* vec );
    /**
     *  Getter _north
     *  @return cell above current cell
     */
    // Cell* North()const;
    /**
     *  Getter _east
     *  @return cell right to current cell
     */
    // Cell* East()const;
    /**
     *  Getter _south
     *  @return cell beneath current cell
     */
    // Cell* South()const;
    /**
     *  Getter _west
     *  @return cell left to current cell
     */
    // Cell* West()const;
    /**
     *  Setter _north
     *  @return cell above current cell
     */
    // void SetNorth(Cell* in);
    /**
     *  Setter _east
     *  @return cell right to current cell
     */
    // void SetEast(Cell* in);
    /**
     *  Setter _south
     *  @return cell beneath current cell
     */
    // void SetSouth(Cell* in);
    /**
     *  Setter _west
     *  @return cell left to current cell
     */
    // void SetWest(Cell* in);

  private:
    Cell();    ///< invalid constructor
};

#endif
