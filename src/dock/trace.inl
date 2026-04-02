// trace.inl
// implementation of class Trace inline methods
// 
static char Trace_inline[] = "$Id: trace.inl,v 1.6 2016/04/26 20:05:10 sbrozell Exp $" ;

#include <iostream>

inline
bool
Trace::isTracingOn( )
{
    return traceIsActive ;
}

inline
Trace::Trace( char const * name )
: theFunctionName( 0 )
{
    if ( isTracingOn() && name != 0 ) {
        indentation += indentation_spacing ;
        std::clog << preamble << indentation << "Enter " << name << std::endl ;
        theFunctionName = new std::string( name ) ;
    }
}

inline
Trace::Trace( std::string const & name )
: theFunctionName( 0 )
{
    if ( isTracingOn() ) {
        indentation += indentation_spacing ;
        std::clog << preamble << indentation << "Enter " << name << std::endl ;
        theFunctionName = new std::string( name ) ;
    }
}

inline
Trace::~Trace()
{
    if ( isTracingOn() && theFunctionName != 0 ) {
        std::clog << preamble << indentation << "Exit  " << * theFunctionName << std::endl ;
        indentation.erase( indentation.size() - indentation_spacing.size() ) ;
        delete theFunctionName ;
    }
}

inline
void
Trace::boolean( char const * message, bool b )
const
{
    if ( isTracingOn() ) {
        std::clog << preamble << indentation << message << " = " << b << std::endl ;
    }
}

inline
void
Trace::note( char const * message )
const
{
    if ( isTracingOn() && message != 0) {
        std::clog << preamble << indentation << message << std::endl ;
    }
}

inline
void
Trace::note( std::string const & message )
const
{
    if ( isTracingOn() ) {
        std::clog << preamble << indentation << message << std::endl ;
    }
}

inline
void
Trace::integer( char const * message, int i )
const
{
    if ( isTracingOn() ) {
        std::clog << preamble << indentation << message << " = " << i << std::endl ;
    }
}

inline
void
Trace::identity( void const * address )
const
{
    if ( isTracingOn() && address != 0) {
        std::clog << preamble << indentation << "My address is " << address << std::endl ;
    }
}

inline
void
Trace::traceOff( )
{
    traceIsActive = false ;
}

inline
void
Trace::traceOn( )
{
    traceIsActive = true ;
}

inline
void
Trace::tracing( bool on )
{
    traceIsActive = on ;
}

