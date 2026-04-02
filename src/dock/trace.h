// 
// definition of class Trace
// 
// This class logs diagnostic information.
// 
#ifndef TRACE_H
#define TRACE_H
static char Trace_header[] = "$Id: trace.h,v 1.7 2016/04/26 20:05:09 sbrozell Exp $" ;

#include <string>

class Trace {

public:

    Trace( char const * name ) ;
    Trace( std::string const & name ) ;
    // No default constructor
    // No copy constructor:  Do Not define
    Trace( Trace const & ) ;
    ~Trace() ;
    // No assignment operator:  Do Not define
    Trace & operator=( Trace const & ) ;

    void note( char const * message ) const ;
    void note( std::string const & message ) const ;
    void boolean( char const * message, bool b ) const ;
    void integer( char const * message, int i ) const ;
    void identity( void const * address ) const ;

    static bool isTracingOn() ;
    static void traceOff() ;
    static void traceOn() ;
    static void tracing( bool on ) ;

private:

    static std::string indentation ;
    static std::string const indentation_spacing ;
    static std::string const preamble ;
    static bool traceIsActive ;

    std::string * theFunctionName ;
};

#include "trace.inl"

#endif  // TRACE_H
