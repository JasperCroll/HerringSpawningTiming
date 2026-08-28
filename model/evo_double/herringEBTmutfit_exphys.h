/***
   This file contains model specific settings for corresponding to the C file with the same name.
***/

/*
 *===========================================================================
 *
 * #define POPULATION_NR	(Required)
 * #define I_STATE_DIM		(Required)
 * #define I_CONST_DIM		(Required)
 * 	Defines the number of structured populations in the problem, the
 *	dimension of the i-state and the number of constant variables that
 *	characterize a cohort. These dimensions should be the same for all the
 *	populations. 
 *
 * #define ENVIRON_DIM		(Required)
 * 	Defines the number of variables characterizing the environment.
 *
 * #define OUTPUT_VAR_NR	(Required)
 * 	Defines the number of quantities that have to be written to the output
 *	file each time that output should be generated.
 *
 * #define PARAMETER_NR		(Optional)
 * 	Define the number of free parameters in the problem. These parameters
 *	can be changed between various runs without compilation of the
 *	program. Parameters that are fixed in the problem can be defined as
 *	constants in the problem-specific program file written by the
 *	user. Changing of these constants requires a new compilation of the
 *	program before use. 
 *
 * #define TIME_METHOD		(Optional)
 *
 *	Possible values: RK2, RK4, RKF45, RKCK, RKTSI5, DOPRI5, DOPRI8 or RADAU5
 *
 * 	The control macro TIME_METHOD can be set to either RK2, RK4, RKF45,
 *	RKCK (which is the default value), RKTSI5, DOPRI5, DOPRI8, or RADAU5
 *	indicating the time integration method to use:  
 *
 *    - the 2nd order Runge-Kutta integration method with fixed step size,
 *    - the 4th order Runge-Kutta integration method with fixed step size,
 *    - the Runge-Kutta-Fehlberg 4/5th order integration method with an
 * 	    adaptive step size,
 *    - an explicit Runge-Kutta integration 4(5) order integration method 
 *      due to Tsitouras (2011) with step size control and dense output.
 *    - the Cash-Karp Runga-Kutta integration method with adaptive step size, 
 *    - DOPRI5: an explicit Runge-Kutta method of order (4)5 due to Dormand
 *	    & Prince with step size control and dense output.
 *    - DOPRI8: an explicit Runge-Kutta method of order 8(5,3) due to Dormand
 *	    & Prince with step size control and dense output.
 *    - RADAU5: an implicit Runge-Kutta method of order 5 with step size
 *	    control and dense output.
 *
 * 	The RKTSI5, DOPRI5, DOPRI8 and RADAU5 methods can detect and locate
 *	discontinuities or events. These events are signalled by the routine
 * 	EventLocation() in the program definition file. Integration will be
 *	carried out exactly up to the moment that the event takes place and
 *	will be restarted subsequently.
 *
 * #define EVENT_NR		(Optional)
 *	Defines the number of different event indicators to track. Events are
 *	signalled by indicators, the value of which is set by the routine
 *	EventLocation(). A zero (0) value of an event indicator flags the
 *	timing of an event. 
 *
 * 	NB: For event location the RKTSI5, DOPRI5, DOPRI8 or RADAU5 method is needed.
 *	    If EVENT_NR is set to a non-zero, positive value, TIME_METHOD will
 *	    automatically be changed to RKTSI5, if TIME_METHOD is not already
 *	    set to RKTSI5, DOPRI5, DOPRI8 or RADAU5. 
 *
 * #define DYNAMIC_COHORTS	(Optional)  (0|1)
 * 	The macro DYNAMIC_COHORTS will determine when new cohorts are to be
 *	generated. After a discontinuity/event has been detected on the basis
 *	of the routine EventLocation(), the routine ForceCohortEnd() will be
 *	called. This routine should return a value of 0 or 1, indicating
 *	whether the current boundary cohort should be closed at the occurrence
 *	of the current event or not. If DYNAMIC_COHORTS equals 0, new cohorts
 *	are formed whenever the return value of the function ForceCohortEnd()
 *	equals 1, IN ADDITION TO the usual formation of new cohorts at
 *	equidistant moments in time. If DYNAMIC_COHORTS equals 1, new cohorts
 *	are only formed whenever the return value of the function
 *	ForceCohortEnd() equals 1. 
 *
 * 	NB: For dynamic cohort closure using the routine ForceCohortEnd() the
 *	    DOPRI5, DOPRI8 or RADAU5 method is needed. If DYNAMIC_COHORTS is
 *	    set to a non-zero, positive value, TIME_METHOD will automatically
 *	    be changed to DOPRI5, if TIME_METHOD is not already set to DOPRI5,
 *          DOPRI8 or RADAU5. 
 *
 * #define BIFURCATION		(Optional)  (0|1)
 *      If the macro BIFURCATION is set to 1 the program will carry out an
 *      integration run over a very long time, during which it will stepwise
 *      increase one parameter at regular intervals. The duration of these
 *      intervals, during which the bifurcation parameter is constant, equals
 *      the maximum integration time that is specified in the CVF file. Because
 *      the bifurcation parameter is increased during a single integration run,
 *      the initial state of the population for a particular value of the
 *      bifurcation parameter is the final state that was reached for the
 *      previous value of the parameter.
 *
 *	To make use of this bifurcation feature of the program, the CVF file
 *	has to contain 6 additional lines with control variables following the
 *	specification of all model parameters. These lines should, for example,
 *	look like:
 *
 *	"Index of the bifurcation parameter"	       			3.0
 *	"Step size in the bifurcation parameter"			0.1
 *	"Final value of the bifurcation parameter"		       10.0
 *	"Linear (0) or logarithmic (1) parameter change"        	0.0
 *	"Period of output generation in bifurcation"		      100.0
 *	"Period of state output in bifurcation"				0.0
 *
 *	The first additional control variable specifies the index number of the
 *	parameter that should be changed, while the second and third determine
 *	the stepsize and the final value of the parameter. The fourth control
 *	variable determines whether parameter steps should be additive (on a
 *	linear scale) or multiplicative (on a logarithmic scale).
 *	The fifth and sixth additional control variable determine the time
 *	period within the integration interval for 1 particular parameter
 *	value, during which normal and state output is generated. For example,
 *	if the parameter is changed every 250 time units, the above value
 *	determine that normal output is discarded for the first 150 of these
 *	250 time units, but written to file for the last 100 time units, while
 *	state output is only generated at the very end of the 250 time units.
 *
 */

#define POPULATION_NR	4
#define I_STATE_DIM	    2
#define I_CONST_DIM	    9
#define ENVIRON_DIM	    4
#define OUTPUT_VAR_NR	19
#define PARAMETER_NR	77
#define TIME_METHOD	    RKCK
#define EVENT_NR       	0
#define DYNAMIC_COHORTS 0
#define BIFURCATION	    1



/*==========================================================================*/
