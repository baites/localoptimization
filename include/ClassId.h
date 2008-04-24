/*
 *  ClassId.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/26/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CLASSID_H
#define CLASSID_H

#include<string>

//! Local optimization library.
namespace LO { 

//! Implementation of a primitive class id. 
class ClassId {

public:

  //! Compare if two objects are the same class.
  bool sameClass(const ClassId & object) const
  {
    return id() == object.id();
  } 

  //! Return an string with object type.
  virtual std::string const id() const = 0;
};

#define DefineClassIdentity(name) std::string const id() const { return name; }

}

#endif
	