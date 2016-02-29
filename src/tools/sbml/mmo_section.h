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

#ifndef MMO_SECTION_H_
#define MMO_SECTION_H_

#include <map>

#include <mmo_decl.h>
#include <mmo_exp.h>

using namespace std;

/**
 *
 */
typedef enum
{
  equation,          //!< equation
  algorithm,         //!< algorithm
  initial_algorithm, //!< initial_algorithm
  external_functions, //!< external_functions
  declarations,      //!< declarations
  imports            //!< imports
} MMOSectionType;

/**
 *
 */
class MMOSection : public MMOExp
{
public:
  /**
   *
   */
  MMOSection () :
    _type()
  {
  }
  ;
  /**
   *
   * @param type
   */
  MMOSection (MMOSectionType type);
  /**
   *
   */
  ~MMOSection ();
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
  getId ()
  {
    return (_id);
  }
  ;
  /**
   *
   * @param exp
   */
  void
  add (MMOExp *exp);
  /**
   *
   * @param id
   * @return
   */
  MMOExp*
  find (string id);
  /**
   *
   * @param id
   * @return
   */
  MMODecl*
  findDec (string id);
  /**
   *
   * @return
   */
  bool
  isEquation ()
  {
    return (_type == equation);
  }
  ;
  /**
   *
   * @return
   */
  bool
  isAlgorithm ()
  {
    return (_type == algorithm);
  }
  ;
  /**
   *
   * @return
   */
  bool
  isInitialAlgorithm ()
  {
    return (_type == initial_algorithm);
  }
  ;
  /**
   *
   * @return
   */
  bool
  isExternalFunctions ()
  {
    return (_type == external_functions);
  }
  ;
  /**
   *
   * @return
   */
  bool
  isDeclarations ()
  {
    return (_type == declarations);
  }
  ;
  /**
   *
   * @return
   */
  bool
  isEmpty ()
  {
    return (_exps.size () == 0);
  }
  ;
  /**
   *
   * @param algebraics
   */
  void
  setAlgebraics (list<string> algebraics);
  /**
   *
   * @return
   */
  MMODecl *
  first();
  /**
   *
   * @return
   */
  MMODecl *
  next();
  /**
   *
   * @return
   */
  bool
  end();
  /**
   *
   * @param iv
   */
  void
  setInitialValues (bool iv);
private:
  unsigned int _expIt;
  string _id;
  MMOSectionType _type;
  map<string, MMOExp*> _exps;
  list<string> _algebraics;
  bool _setInitialValues;
};

#endif  /* MMO_SECTION_H_ */
