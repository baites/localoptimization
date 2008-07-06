/*
 *  ErrorHandler.h
 *  lo++
 *
 *  Created by Victor Eduardo Bazterra on 5/26/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ERROR_HANDLER_H
#define ERROR_HANDLER_H

#include <string>
#include <iostream>

namespace LO {

static std::ostream* stream = &std::cerr;

//! Primitive implementation for locating warnings and errors in the code.
class SourceLocator {

public:

  //! Line number of a given location in the code.
  long line;
  
  //! File name of a given location in the code.
  std::string file;
  
  //! Constructor by file and line number of a given location in the code.
  SourceLocator(std::string const & _file, long _line) : file(_file), line(_line) {};

};

#define LOCATION SourceLocator(__FILE__, __LINE__)

//! This function is called when a error in a function occurs.
/*
  The user needs to call this function when an error ocurrs inside of a function
  that do not beling to any class.
  \param location MACRO to mark the location in the code.
  \param function name.
  \param error message.  
*/
void FunctionError (
  SourceLocator const location,
  std::string const & function,
  std::string const & message
)
{
  *stream << function << ":" << std::endl;
  *stream << message << std::endl;
  *stream << location.file << " : " << location.line << std::endl;   
}

//! This function is called when a error in a class function occurs.
/*
  The user needs to call this function when an error ocurrs inside of a function
  that do belong to a class.
  \param location MACRO to mark the location in the code.
  \param class name.
  \param function name.
  \param error message.  
*/
void ClassError(
  SourceLocator const location,
  std::string const & classname,
  std::string const & function,
  std::string const & message
)
{
  *stream << classname << "::" << function << ":" << std::endl;
  *stream << message << std::endl;
  *stream << location.file << " : " << location.line << std::endl;   
}

}

#endif
