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

#ifndef MMO_CONVERT_H_
#define MMO_CONVERT_H_

#include "mmo_model.h"
#include "mmo_writer.h"

#include <sbml/SBMLVisitor.h>
#include <string> 

using namespace std;

class MMOConvert : public SBMLVisitor
{
public:
  /**
   *
   * @param name
   * @param replace
   * @param type
   */
  MMOConvert (string prefix, bool replace = false, WR_Type type = WR_MMO);
  /**
   *
   */
  ~MMOConvert ();
  /**
   *
   * @return
   */
  MMOModel
  model ()
  {
    return (_model);
  };
  /**
   *
   * @param x
   */
  void
  visit (const SBMLDocument &x);
  /**
   *
   * @param x
   */
  void
  visit (const Model &x);
  /**
   *
   * @param x
   */
  void
  visit (const KineticLaw &x);
  /**
   *
   * @param x
   * @param type
   */
  void
  visit (const ListOf &x, int type);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const SBase &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const FunctionDefinition &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const UnitDefinition &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Unit &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const CompartmentType &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const SpeciesType &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Compartment &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Species &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Parameter &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const InitialAssignment &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Rule &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const AlgebraicRule &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const AssignmentRule &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const RateRule &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Constraint &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Reaction &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const SimpleSpeciesReference &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const SpeciesReference &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const ModifierSpeciesReference &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Event &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const EventAssignment &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Trigger &x);
  /**
   *
   * @param x
   * @return
   */
  bool
  visit (const Delay &x);
  /**
   *
   * @param x
   */
  void
  visit (const Priority &x);
  /**
   *
   * @param x
   */
  void
  leave (const SBMLDocument &x);
  /**
   *
   * @param x
   */
  void
  leave (const Model &x);
  /**
   *
   * @param x
   */
  void
  leave (const KineticLaw &x);
  /**
   *
   * @param x
   */
  void
  leave (const Priority &x);
  /**
   *
   * @param x
   */
  void
  leave (const Reaction &x);
  /**
   *
   * @param x
   */
  void
  leave (const SBase &x);
  /**
   *
   * @param x
   * @param type
   */
  void
  leave (const ListOf &x, int type);
private:
  map<string, ASTNode *> _reaction_def;
  list<Event> _events;
  bool _in_reaction;
  WR_Type _type;
  string _reaction_name;
  MMOModel _model;
};

#endif /* MMO_CONVERT_H_ */
