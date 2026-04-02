// trace.C
// implementation of class Trace
//
// The implementation is trivial.
//
static char Trace_C[] = "$Id: trace.cpp,v 1.2 2015/10/29 04:53:35 sbrozell Exp $" ;

#include "trace.h"

// static initialization
std::string Trace::indentation( "" ) ;
std::string const Trace::indentation_spacing( " " ) ;
std::string const Trace::preamble( "Trace:" ) ;
bool Trace::traceIsActive( false ) ;

