/*****************************************************************************

 This file is part of QSS Solver.

 QSS Solver is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 QSS Solver is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with QSS Solver.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#ifndef MMO_DECL_H_
#define MMO_DECL_H_

#include "mmo_exp.h"
#include "mmo_math.h"
#include "mmo_visitor.h"

#include <string>

using namespace std;

/**
 *
 */
typedef enum
{
  constant,           //!< constant
  parameter,          //!< parameter
  state,              //!< state
  discrete,           //!< discrete
  algebraic,          //!< algebraic
  derivative,         //!< derivative
  assignment,         //!< assignment
  initial_assignment, //!< initial_assignment
  zc_relation,        //!< zc_relation
  zc_oposite_relation,//!< zc_oposite_relation
  algebraic_equation, //!< algebraic_equation
  implicit_equation,  //!< implicit_equation
  reinit,             //!< reinit
  condition,          //!< condition
  function_input,     //!< function_input
  function_output,    //!< function_output
  function_definition,//!< function_definition
  function_formula,   //!< function_formula
  import              //!< import
} MMODeclType;

/**
 *
 */
class MMODecl : public MMOExp
{
public:
  /**
   *
   * @param id
   * @param value
   * @param type
   */
  MMODecl (string id, double value, MMODeclType type);
  /**
   *
   * @param id
   * @param exp
   * @param type
   */
  MMODecl (string id, string exp, MMODeclType type);
  /**
   *
   * @param exp
   * @param type
   */
  MMODecl (string exp, MMODeclType type);
  /**
   *
   */
  MMODecl ();
  /**
   *
   */
  ~MMODecl ();
  /**
   *
   * @param visitor
   */
  void
  accept (MMOVisitor *visitor);
  /**
   *
   * @return
   */
  string
  getId ();
  /**
   *
   * @param i
   */
  void
  id (string i);
  /**
   *
   * @param type
   */
  void
  setType (MMODeclType type);
  /**
   *
   * @return
   */
  string
  getExp ();
  /**
   *
   * @param i
   */
  void
  exp (string i);
  /**
   *
   * @param i
   */
  void
  value (double i);
  /**
   *
   * @return
   */
  double
  getValue ();
  /**
   *
   * @return
   */
  bool
  hasExp ();
  /**
   *
   * @return
   */
  bool
  hasValue ();
  /**
   *
   * @return
   */
  bool
  isAlgebraicEquation ();
  /**
   *
   * @return
   */
  bool
  isInitialAssignment ();
  /**
   *
   * @return
   */
  bool
  isAssignment ();
  /**
   *
   * @return
   */
  bool
  isZeroCrossing ();
  /**
   *
   * @return
   */
  bool
  isOpositeZeroCrossing ();
  /**
   *
   * @return
   */
  bool
  isDerivative ();
  /**
   *
   * @return
   */
  bool
  isParameter ();
  /**
   *
   * @return
   */
  bool
  isConstant ();
  /**
   *
   * @return
   */
  bool
  isState ();
  /**
   *
   * @return
   */
  bool
  isDiscrete ();
  /**
   *
   * @return
   */
  bool
  isAlgebraic ();
  /**
   *
   * @return
   */
  bool
  isCondition ();
  /**
   *
   * @return
   */
  bool
  isFunctionInput ();
  /**
   *
   * @return
   */
  bool
  isFunctionOutput ();
  /**
   *
   * @return
   */
  bool
  isFunctionDefinition ();
  /**
   *
   * @return
   */
  bool
  isFunctionFormula ();
  /**
   *
   * @return
   */
  bool
  isReinit ();
  /**
   *
   * @return
   */
  bool
  isImport();
  /**
   *
   * @return
   */
  bool
  isImplicit ();
private:
  string _id;
  string _exp;
  double _value;
  MMODeclType _type;
  bool _init;
};

#endif  /* MMO_DECL_H_ */
