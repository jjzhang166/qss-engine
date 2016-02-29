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

#include <sbml/math/ASTNode.h>
#include <sbml/math/FormulaFormatter.h>
#include <sbml/math/L3Parser.h>
#include "mmo_utils.h"
#include "mmo_decl.h"
#include "mmo_math.h"

MMODecl::MMODecl () :
_value(),
_type(constant),
_init()
{
}

MMODecl::MMODecl (string id, double value, MMODeclType type)
{
  _id = id;
  _exp = "";
  _type = type;
  _value = value;
  _init = false;
}

MMODecl::MMODecl (string id, string exp, MMODeclType type)
{
  _id = id;
  _exp = exp;
  _type = type;
  _value = 0;
  _init = false;
}

MMODecl::MMODecl (string id, MMODeclType type)
{
  if (type == zc_relation || type == condition)
    {
      _id = "";
      _exp = id;
    }
  else
    {
      _id = id;
      _exp = "";
    }
  _type = type;
  _value = 0;
  _init = false;
}

MMODecl::~MMODecl ()
{
}

void
MMODecl::accept (MMOVisitor *visitor)
{
  visitor->visit (this);
  visitor->leave (this);
}

string
MMODecl::getId ()
{
  return (_id);
}

void
MMODecl::id (string i)
{
  _id = i;
}

void
MMODecl::setType (MMODeclType type)
{
  _type = type;
}

string
MMODecl::getExp ()
{
  return (_exp);
}

void
MMODecl::exp (string i)
{
  _exp = i;
}

void
MMODecl::value (double i)
{
  _value = i;
}

double
MMODecl::getValue ()
{
  return (_value);
}

bool
MMODecl::hasExp ()
{
  return (_exp != "");
}

bool
MMODecl::hasValue ()
{
  return (_exp == "");
}

bool
MMODecl::isAlgebraicEquation ()
{
  return (_type == algebraic_equation);
}

bool
MMODecl::isInitialAssignment ()
{
  return (_type == initial_assignment);
}

bool
MMODecl::isAssignment ()
{
  return (_type == assignment);
}

bool
MMODecl::isZeroCrossing ()
{
  return (_type == zc_relation);
}

bool
MMODecl::isOpositeZeroCrossing ()
{
  return (_type == zc_oposite_relation);
}
bool
MMODecl::isDerivative ()
{
  return (_type == derivative);
}

bool
MMODecl::isParameter ()
{
  return (_type == parameter);
}

bool
MMODecl::isConstant ()
{
  return (_type == constant);
}

bool
MMODecl::isState ()
{
  return (_type == state);
}

bool
MMODecl::isDiscrete ()
{
  return (_type == discrete);
}

bool
MMODecl::isAlgebraic ()
{
  return (_type == algebraic);
}

bool
MMODecl::isCondition ()
{
  return (_type == condition);
}

bool
MMODecl::isFunctionInput ()
{
  return (_type == function_input);
}

bool
MMODecl::isFunctionOutput ()
{
  return (_type == function_output);
}

bool
MMODecl::isFunctionDefinition ()
{
  return (_type == function_definition);
}

bool
MMODecl::isFunctionFormula ()
{
  return (_type == function_formula);
}

bool
MMODecl::isReinit ()
{
  return (_type == reinit);
}

bool
MMODecl::isImport ()
{
  return (_type == import);
}

bool
MMODecl::isImplicit ()
{
  return (_type == implicit_equation);
}
