/***
   This file contains model specific code for a structured herring population feeding on two unstructured resources.

   This model implementation in tailored to calculate the invasion fitness
   of a mutant individual with a spawning startdate one day before, one day after or at the same day as a resident population.

   In this specific code deals with the case in which temperature and light only affect prey dynamics and herring physiology.
***/

/*==========================================================================
 *
 *			INCLUDING THE HEADER FILE
 *
 * Include the header file for standard library prototyping and definition
 * of the data-structures.
 *
 *==========================================================================
 */
#include  "escbox.h"
#include  "stdbool.h"



/*
 *==========================================================================
 *
 *		LABELLING ENVIRONMENT AND I-STATE VARIABLES
 *
 * For convenience it is possible to label the environment and i-state
 * variables with a more meaningful name by defining, for instance: 
 *
 *		#define time	env[0]
 *		#define length	i_state(0)
 *
 * Note in this respect again the zero based array indexing in the
 * C-language.
 *
 *==========================================================================
 */

#define time		env[0]
#define zoopL		env[1]
#define bent		env[2]
#define eggbuffer   env[3]

#define herring		0


#define	indx		i_state(0)
#define indy		i_state(1)

#define	IDage		i_const( 0)
#define	IDmass		i_const( 1)
#define	IDlen		i_const( 2)
#define	IDcond 		i_const( 3)
#define	IDbtime		i_const( 4)
#define	IDbage 		i_const( 5)
#define	IDmatage 	i_const( 6)
#define	IDinit 		i_const( 7)
#define	IDmutindex  i_const( 8)

#define kel 		273.15
#define kB			8.6173E-5

#define TINY		1E-7
#define MAXAGE		7300
#define SCALE	1E6


/*
 *==========================================================================
 *
 *		DEFINING AND LABELLING CONSTANTS AND PARAMETERS
 *
 * Define the constant names and values used in the user-specified
 * routines. Most parameters in the problem can be treated as constants for
 * this file. This seems the most easy way to specify the parameters. The
 * parameters that are free to change will be in a vector called
 * "parameter[]", defined elsewhere. As with the i-state variables it is 
 * possible to label these parameters with a more meaningful name, e.g.:
 *
 *		#define	   ALPHA    parameter[0]
 *
 *==========================================================================
 */

// Resource
#define PZ		parameter[ 0] // Productivity zooplankton
#define RZ		parameter[ 1] // Turnover zooplankton
#define PB		parameter[ 2] // Productivity Benthos
#define RB		parameter[ 3] // Turnover benthos

// Allometric scaling
#define LAMBDA1	parameter[ 4] // allometrix size scaling
#define LAMBDA2	parameter[ 5] // alometric size scaling

// Development
#define EPSMIN	parameter[ 6] // Minimum egg development time
#define EPSP	parameter[ 7] // Temperature dependent development time scalar
#define WB		parameter[ 8] // Weight at hatching
#define AF		parameter[ 9] // Age first feeding
#define LM		parameter[10] // Length at maturation

// Condition
#define QJ		parameter[11] // Maximum juvenile condition
#define QA		parameter[12] // Maximum adult condition
#define QR		parameter[13] // Minimum conditon for reproduction
#define QS		parameter[14] // Minimum condition before starvation

// Ingestion
#define D1		parameter[15] // Allometric scalar digestion
#define D2		parameter[16] // Allometric scalar digestion
#define AMAX	parameter[17] // Maximum attack rate
#define ALPHA	parameter[18] // Allometric scalar zooplankton attack rate
#define B1		parameter[19] // Allometric scalar benthic attack rate
#define B2		parameter[20] // Allometric scalar benthic attack rate

#define BMAX	parameter[21] // Maximum benthic consumption

#define WOPT	parameter[22] // Optimum zooplankton consumption mass
#define LH		parameter[23] // size at half the maximum benthic consumption

#define SB		parameter[24] // feeding change slope scalar

// Energy allocation
#define KE		parameter[25] // Food energy conversion
#define RHO1	parameter[26] // Allometric scalar metabolism
#define RHO2	parameter[27] // Allometric scalar metabolism

// Reproduction
#define KR		parameter[28] // Gonal conversion efficiency
#define SM		parameter[29] // Time with maximum spawning after start
#define SD		parameter[30] // Duration of spawning period
#define SSTART	parameter[31] // Start of spawning period

// Mortality
#define MUL		parameter[32] // feeding change slope scalar
#define MUS		parameter[33] // Size dpeendent mortality constant 
#define SC		parameter[34] // Starvation rate scalar
#define MUC		parameter[35] // Size dependent mortality scalar
#define MU0		parameter[36] // Background mortality

// Fisheries
#define MUF		parameter[37] // Fisheries mortality scalar
#define LHF		parameter[38] // Size half maximum fishing selectivity
#define SF		parameter[39] // Fisheries selectivity slope
#define LFMIN	parameter[40] // minimum fishing size

// Temperature and light dependence
#define GMZ1	parameter[41] // Light dependence zooplankton
#define GMZ2	parameter[42] // Light dependence slope zooplankton
#define LHZ		parameter[43] // Light intensity with half maximum zooplankton productivity

#define TR		parameter[44] // Reference temperature

#define EAR		parameter[45] // Activation energy resource turnover rate

#define EAP		parameter[46] // Activation energy resource production
#define EDP		parameter[47] // Deactivation energy resource production
#define THP		parameter[48] // Temperature half deactivation resource production

#define EAD		parameter[49] // Activation energy digestion
#define EDD		parameter[50] // Deactivation energy digestion
#define THD		parameter[51] // Temperature half deactivation digestion

#define EAM		parameter[52] // Activation energy metabolic rate
#define EI		parameter[53] // Size temperature interaction metabolic rate

#define EAS		parameter[54] // Activation energy egg survival
#define EDS		parameter[55] // Deactivation energy egg survival
#define THS		parameter[56] // Temperature half deactivation egg survival

#define EAE		parameter[57] // Activatio nenergy egg duration

#define TAC		parameter[58] // Temperature amplitude center
#define TAA		parameter[59] // Temperature amplitude amplitude
#define TAS		parameter[60] // Temperature amplitude shift
#define TSC		parameter[61] // Temperature shift center
#define TSA		parameter[62] // Temperature shift amplitude
#define TSS		parameter[63] // Temperature shift shift
#define TTC		parameter[64] // Temperature center

#define LC		parameter[65] // Light center
#define LA		parameter[66] // Light amplitude
#define LS		parameter[67] // Light shift

#define WPROD	parameter[68] // fraction productivity in winter without light.

// scalars
#define VB		parameter[69] // Water srufance volume scalar
#define HWSTART		parameter[70] // Heatwave start of max increase period
#define HWDUR		parameter[71] //Heatwave duration of max increase period
#define HWINC		parameter[72] // Heatwave increase
#define HWSLOPE		parameter[73] // heatwave start and end duration

#define LHY		parameter[74] // Size half maximum fishing selectivity
#define SY		parameter[75] // Fisheries selectivity slope
#define LYMIN	parameter[76] // minimum fishing size

#define MUTINIT 1E9

#define MUTYEAR 90
#define MUTPERIOD 36
#define BIFPERIOD 150

double LROcoho[3][MUTPERIOD];
double gentime[3][MUTPERIOD];

#define MEANTEMP			7.3615
#define MEANLIGHTPROD		0.6061





void UpdateIDcards(double* env, population* pop);



double gettemp(double reltime){ 
	double curtime = reltime + SSTART;
	double Ta = TAC + TAA * sin(2 * M_PI / 365 * (curtime - TAS));
	double Ts = TSC + TSA * sin(2 * M_PI / 365 * (curtime - TSS));
	double settemp = TTC + Ta*sin(2 * M_PI / 365 * (curtime - Ts));
	//double settemp = 6.482;

	return settemp; 
}

double getlight(double reltime){
	double curtime = reltime + SSTART;
	double setlight = LC + LA * cos(2 * M_PI / 365 * (curtime - LS));
	return setlight;
}

double getheatwave(double reltime) {
	double curtime = reltime + SSTART;



	if (curtime <= HWSTART - HWSLOPE) return 0;
	if (curtime >= HWSTART + HWDUR + HWSLOPE) return 0;

	if (curtime >= HWSTART - HWSLOPE && curtime <= HWSTART) return HWINC - (HWSTART - curtime) / HWSLOPE * HWINC;
	if (curtime >= HWSTART && curtime <= HWSTART + HWDUR) return HWINC;
	if (curtime >= HWSTART + HWDUR && curtime <= HWSTART + HWDUR + HWSLOPE) return (HWSLOPE - (curtime - HWSTART - HWDUR)) / HWSLOPE * HWINC;
	



	return 0;

}

double gam(double light, double gamma1, double gamma2, double lH){
	double val= WPROD + (1-WPROD)/(1+gamma1*exp(-gamma2*(light-lH)));

	return val;
}

double phiarr(double T, double Ea, double Tr){
	double val=exp(-Ea/kB*(1/(T+kel)-(double)1/(Tr+kel)));
	return val;
}

double phiSS(double T, double Ea, double Ed, double Tr, double Td){
	double val = exp(-Ea/kB*(1/(T+kel)-(double)1/(Tr+kel)))/((double)1+exp(-Ed/kB*(1/(T+kel)-(double)1/(Td+kel))));
	return val;
}



/*
 *==========================================================================
 *
 * USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *==========================================================================
 */

void	UserInit(int argc, char** argv, double* env, population* pop)

/*
 * The function header above should not be changed.
 *
 * In the routine below the user can specify some operations that are to
 * be performed before the actual integration starts. Initialization of
 * user defined variables can, for instance, be carried out here.
 *
 *     NB : -	The integer array "cohort_no[]" is globally available and
 *		denotes the number of internal cohorts present in each
 *		structured population.
 */
{

	switch (argc) {
	case 5:
		HWINC = strtoll(argv[4], NULL, 10);
		printf("HWINC set to %f\n", HWINC);
	case 4:
		HWDUR = strtoll(argv[3], NULL, 10);
		printf("HWDUR set to %f\n", HWDUR);
	case 3:
		HWSTART = strtoll(argv[2], NULL, 10);
		printf("HWSTART set to %f\n", HWSTART);
	}


	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < cohort_no[i]; j++) {
			popIDcard[i][j][IDbtime] = time - popIDcard[i][j][IDage] - popIDcard[i][j][IDbage];
		}
	}




	

	return;
}



/*
 *==========================================================================
 *
 *	SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *==========================================================================
 */

void	SetBpointNo(double *env, population *pop, int *bpoint_no)
  /*
   * This function can also be defined as:

   void	CycleStart(double *env, population *pop, int *bpoint_no)

   * The function names SetBpointNo and CycleStart are synonyms.
   *
   * Apart from this, the function header above should not be changed.
   *
   * Specify below the number of boundary cohorts that should be created at
   * the start of the next cohort integration cycle. Fill the array
   * "bpoint_no[]" with the appropriate integer values. The length of the
   * array is "POPULATION_NR".
   * The state of the environment and the population can be used to adapt
   * the number of the fixed points on the boundary to the current state.
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   */
  
{
	bpoint_no[0]=0;
	bpoint_no[1] = 0;	
	bpoint_no[2] = 0;
	bpoint_no[3] = 0;

	
	return;
}



/*==========================================================================*/

void	SetBpoints(double *env, population *pop, population *bpoints)
  
  /*
   * The function header above should not be changed.
   *
   * Define below the fixed points at the boundary of the state space that
   * will characterize the boundary cohorts in the subsequent cohort cycle.
   * These boundary points should be addressed in the same way as the
   * population cohorts and the offspring cohorts, i.e.:
   *
   *		      bpoints[1][4][i_state(2)]
   *
   * to denote the fixed value of the 3rd i_state variable of the 5th
   * boundary point in the 2nd population (Note again the base of 0 of an
   * array in the C-language). 
   * Set the fixed value for all the i_state variables in all the boundary
   * cohorts for all populations. The number of populations is given by
   * "POPULATION_NR", the number of boundary cohorts is given by the
   * globally available vector "bpoint_no[]" of length "POPULATION_NR" and
   * the number of i_state variables is obviously determined by
   * "I_STATE_DIM". 
   * The state of the environment and the population can be used to adapt
   * the values of the fixed points on the boundary to the current state.
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   */
{
	return;
}



/*
 *==========================================================================
 *
 *			SPECIFICATION OF DERIVATIVES
 *
 *==========================================================================
 */

void	Gradient(double *env,     population *pop,     population *ofs,
		 double *envgrad, population *popgrad, population *ofsgrad,
		 population *bpoints)
  
  /*
   * The function header above should not be changed.
   *
   * Define the derivatives of the various environment variables,
   * which have to be returned to the main program in the array
   * "envgrad[]". 
   * Define also the derivatives of the various population variables,
   * which have to be returned to the main program in the matrix of cohort
   * variables "popgrad[][][]" for each population.
   * Finally define the derivatives of the various offspring variables,
   * which have to be returned to the main program in the matrix of cohort
   * variables "ofsgrad[][][]" for each population.
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   *	      - The integer array "bpoint_no[]" is globally available and
   *		denotes the number of boundary cohorts present in each
   *		structured population during the current cohort cycle.
   *	      - The offspring number can be zero. If used in a division, check
   *		for this equality to zero!! 
   */

{

	UpdateIDcards(env, pop);

	register int i;
	double cond, mass, len;
	double mu_size, mu_starv, mu_fish, mu;
	double fall, Eg, gx, gy, Ea, Em;
	double IzS, IzL, Ib, BlB, BlZ, hand, etab, etazS, etazL, Az, Ab;

	double consZS = 0;
	double consZL = 0;
	double consB = 0;

	double temp = gettemp(time) + getheatwave(time);
	double light = getlight(time);



	for (int j = 0; j < 4; j++) {
		for (i = 0; i < cohort_no[j]; i++) {

			// default
			popgrad[j][i][number] = 0.0;

			popgrad[j][i][indx] = 0.0;
			popgrad[j][i][indy] = 0.0;

			if (popIDcard[j][i][IDage] < -popIDcard[j][i][IDbage]) {
				continue;
			}

			if (popIDcard[j][i][IDage] < 0) { //nonfeeding stage
				popgrad[j][i][number] = -MUL * pop[j][i][number];
				continue;
			}

			if (pop[j][i][indy] < 0) continue;


			mass = popIDcard[j][i][IDmass];
			len = popIDcard[j][i][IDlen];
			cond = popIDcard[j][i][IDcond];

			// mortality
			mu_size = MUS * exp(-mass / MUC);
			mu_starv = 0;
			if (len > LFMIN) mu_fish = MUF / (1 + exp(-SF * (len - LHF)));
			else mu_fish = 0;

			if (cond < QS) mu_starv = SC * (QS / cond - 1);

			mu = MU0 + mu_starv + mu_size + mu_fish;


			//ingestion
			Az = phiarr(temp, EAM, TR) * pow(mass, EI / kB * (1 / (temp + kel) - 1 / (TR + kel))) * AMAX * pow((mass / WOPT * exp(1 - mass / WOPT)), ALPHA);
			Ab = phiarr(temp, EAM, TR) * B1 * pow(mass, B2 + EI / kB * (1 / (temp + kel) - 1 / (TR + kel)));


			etazL = max(0, Az * zoopL); // netagive value protection)
			etab = max(0, Ab * bent);

			hand = D1 * pow(mass, D2) / phiSS(temp, EAD, EDD, TR, THD);
			BlB = BMAX / (1 + exp(-SB * (len - LH)));

			IzL = (1 - BlB) * etazL / (1 + hand * ((1 - BlB) * etazL + BlB * etab));
			Ib = BlB * etab / (1 + hand * ((1 - BlB) * etazL + BlB * etab));

			if (j == 0) {
				consZL += max(0, IzL * pop[j][i][number] * zoopL);
				consB += max(0, Ib * pop[j][i][number] * bent);
			}

			Ea = KE * (IzL + Ib);


			Em = phiarr(temp, EAM, TR) * pow(mass, EI / kB * ((double)1 / (temp + kel) - (double)1 / (TR + kel))) * RHO1 * pow(pop[j][i][indx] + pop[j][i][indy], RHO2);
			Eg = Ea - Em;


			if (len < LM && cond < QJ) fall = cond * cond / ((1 + QJ) * QJ * QJ);
			if (len < LM && cond >= QJ) fall = (double)1 / (1 + QJ);
			if (len >= LM && cond < QA) fall = cond * cond / ((1 + QA) * QA * QA);
			if (len >= LM && cond >= QA) fall = (double)1 / (1 + QA);

			gx = fall * Eg;
			gy = (1 - fall) * Eg;


			if (Eg < 0) {
				gx = 0;
				gy = Eg;
			}

			// final gradients
			popgrad[j][i][number] = -mu * pop[j][i][number];

			popgrad[j][i][indx] = gx;
			popgrad[j][i][indy] = gy;



		}
	}


	envgrad[0] = 1;//time
	envgrad[1] = max( 0, gam(light, GMZ1, GMZ2, LHZ) * phiSS(temp, EAP, EDP, TR, THP) * PZ) -  phiarr(temp, EAR, TR) * RZ * zoopL - consZL/SCALE; // zooplankton
	envgrad[2] = max(0, phiSS(temp, EAP, EDP, TR, THP) * PB ) - phiarr(temp, EAR, TR) * RB * bent - consB * VB/SCALE; // benthos
	envgrad[3] = 0; // egg buffer
	 
	return;
}



/*
 *==========================================================================
 *
 *	SPECIFICATION OF EVENT LOCATION AND DYNAMIC COHORT CLOSURE
 *
 *==========================================================================
 */

void	EventLocation(double *env, population *pop, population *ofs,
		      population *bpoints, double *events)
  
  /*
   * The function header above should not be changed. 
   *
   * When the value of EVENT_NR is larger than 0, the current routine is called
   * to check whether a discontinuity in the time derivatives or an event has
   * occurred. Such a discontinuity or event is indicated by a zero (0) value
   * of one of the elements of the array events[]. Hence, every time the
   * value of one of these array elements changes from positive to negative or
   * vice versa within a time integration step, this is interpreted as an event
   * occurrence. If an event occurs, the integration will not step over the
   * time point of the event, but will localize the moment of its occurrence by
   * means of a root finding procedure. This localization will stop when it has
   * found a time, for which the value of the particular array element that
   * changed sign is within the tolerance bounds, distinguishing values from
   * zero, as set in the CVF file. After the event time has been localized,
   * integration will be carried out exactly up to the moment of event
   * occurrence. Subsequently the routine ForceCohortEnd() (see below) is
   * called and the integration is started from the moment of the event
   * location onwards. If no event location is required, this routine should
   * not change the value of the elements in events[] or the macro EVENT_NR
   * should be set to 0. 
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   *	      - The integer array "bpoint_no[]" is globally available and
   *		denotes the number of boundary cohorts present in each
   *		structured population during the current cohort cycle.
   *	      - The offspring number can be zero. If used in a division, check
   *		for this equality to zero!! 
   */
  
{
	return;
}




/*==============================================================================*/

int	ForceCohortEnd(double *env, population *pop, population *ofs,
		       population *bpoints)
  
  /*
   * The function header above should not be changed. 
   *
   * When the value of EVENT_NR is larger than 0, the current routine is called
   * after an event has been localized on the basis of the EventLocation()
   * routine defined above. This routine should then return either COHORT_END
   * (=1) or NO_COHORT_END (=0), indicating whether a cohort closure is
   * currently required or not. If the macro DYNAMIC_COHORTS in the program
   * header file is set to 0 (its default value) closing off boundary cohorts
   * and transforming them into internal cohorts, takes place both at constant
   * time intervals, as specified by the cohort cycle limit variable in the
   * .CVF file and whenever the current routine returns COHORT_END. This allows
   * for dynamic closure of boundary cohorts, for instance, to ensure a limited
   * number of individuals per cohort. If DYNAMIC_COHORTS is set to 1 in the
   * program header file, closing off of boundary cohorts is exclusively
   * performed when the current routine returns COHORT_END.
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   *	      - The integer array "bpoint_no[]" is globally available and
   *		denotes the number of boundary cohorts present in each
   *		structured population during the current cohort cycle.
   *	      - The offspring number can be zero. If used in a division, check
   *		for this equality to zero!! 
   */
  
{
  

  /* return COHORT_END;    */
  return NO_COHORT_END; 
}



/*
 *==========================================================================
 *
 *		SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *==========================================================================
 */

void	InstantDynamics(double *env, population *pop, population *ofs)
  
  /*
   * This function can also be defined as:

   void	CycleEnd(double *env, population *pop, population *ofs)

   * The function names SetBpointNo and CycleStart are synonyms.
   *
   * Apart from this, the function header above should not be changed.
   *
   * This routine is called at the end of each cohort cycle. The routine
   * can be used to implement any type of instantaneous dynamics, occurring
   * between two subsequent cohort cycles. It can, for instance, be used
   * for a pulsed, instantaneous reproduction process or to set the number
   * of individuals in cohorts that have reached their maximum lifespan to
   * 0. 
   * Note that the transformation of boundary cohorts into internal cohorts
   * has already been performed and that the boundary cohorts are hence
   * characterized by the number of individuals and their transformed
   * moments (usually denoted by the symbol mu).
   * In an instantaneous reproduction process the i_state of the offspring
   * can therefore be simply specified in terms of the mean i_state!
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   *	      - The integer array "bpoint_no[]" is globally available and
   *		denotes the number of boundary cohorts present in each
   *		structured population during the current cohort cycle.
   *	      - The offspring number can be zero. If used in a division, check
   *		for this equality to zero!! 
   */

{


	if (fmod(time,  BIFPERIOD * 365) < 1 + TINY && fmod(time, BIFPERIOD * 365) > 1 - TINY) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < MUTPERIOD; j++) {
				LROcoho[i][j] = 0;
				gentime[i][j] = 0;
			}

		}
	}


	UpdateIDcards(env, pop);

	

	register int i;

	double mass, len, temp;

	// prevent negative environment
	if (env[1] < 0) env[1] = 0;
	if (env[2] < 0) env[2] = 0;
	if (env[3] < 0) env[3] = 0;

	// set some cohorts to zero
	for (int j = 0; j < 4; j++) {
		for (i = 0; i < cohort_no[j]; i++) {
			if (popIDcard[j][i][IDage] > MAXAGE || popIDcard[j][i][IDcond] < 0) pop[j][i][number] = -1;
		}
	}



	if (fmod(time, 365) > 364 - TINY && fmod(time, 365) < 364 + TINY) { // early mutant
		temp = gettemp(time) + getheatwave(time);


		// collect eggs mutant
		for (i = 0; i < cohort_no[2]; i++) {
			if (pop[2][i][indy] / pop[2][i][indx] > QR && popIDcard[2][i][IDlen] > LM && pop[2][i][number] > 0) {

				double yearrep = KR * (pop[2][i][indy] - QR * pop[2][i][indx]) / WB * pop[2][i][number] / MUTINIT;
				int cohoindex = (int)round(popIDcard[2][i][IDmutindex]);

				LROcoho[1][cohoindex] += yearrep;
				gentime[1][cohoindex] += popIDcard[2][i][IDage] * yearrep;

				pop[2][i][indy] = QR * pop[2][i][indx];

			}
		}

		// new cohorts

		if (fmod(time+1, BIFPERIOD * 365) > MUTYEAR * 365 - TINY && fmod(time+1, BIFPERIOD * 365) < (MUTYEAR + MUTPERIOD - 1) * 365 + TINY) {
			// cohorts

			int NewCohort = AddCohorts(pop, 2, SD + 1);

			for (i = 0; i < SD + 1; i++) {
				temp = gettemp(time + i) + getheatwave(time + i);

				if (i == 0)			pop[2][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * 1.0 / (1 + exp(-(0 - SM))) * MUTINIT;// first reproduction day to make distribution add to 1
				else if (i == SD) 		pop[2][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 - 1.0 / (1 + exp(-(SD - 1 - SM)))) * MUTINIT; // last reproduction day to make distribution add to 1
				else 				pop[2][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS)* (1.0 / (1 + exp(-(i - SM))) - 1.0 / (1 + exp(-(i - 1 - SM))))* MUTINIT;

				popIDcard[2][NewCohort + i][IDbtime] = time + i; // egg release time
				popIDcard[2][NewCohort + i][IDbage] = (AF + EPSMIN + EPSP * phiarr(MEANTEMP, EAE, TR)); // egg and larvea duration

				popIDcard[2][NewCohort + i][IDage] = -i - popIDcard[2][NewCohort + i][IDbage]; // the age is zero as soon as individuals turn into juveniles
				pop[2][NewCohort + i][indx] = WB;
				pop[2][NewCohort + i][indy] = QJ * WB;

				popIDcard[2][NewCohort + i][IDmass] = (1 + QR) * pop[2][NewCohort + i][indx];
				popIDcard[2][NewCohort + i][IDlen] = LAMBDA1 * pow(popIDcard[2][NewCohort + i][IDmass], LAMBDA2);
				popIDcard[2][NewCohort + i][IDcond] = pop[2][NewCohort + i][indy] / pop[2][NewCohort + i][indx];

				popIDcard[2][NewCohort + i][IDmatage] = MAXAGE;

				popIDcard[2][NewCohort + i][IDinit] = -1;

				popIDcard[2][NewCohort + i][IDmutindex] = round(fmod((time + 1) / 365, BIFPERIOD) - MUTYEAR);
			}
		}
	}



	if (fmod(time, 365) > 365 - TINY || fmod(time, 365) < TINY) {



		double temp = gettemp(time) + getheatwave(time);


		// collect eggs resudent
		eggbuffer = 0;
		for (i = 0; i < cohort_no[0]; i++) {
			if (pop[herring][i][indy] / pop[herring][i][indx] > QR && popIDcard[herring][i][IDlen] > LM && pop[0][i][number] > 0) {
				eggbuffer += KR * (pop[herring][i][indy] - QR * pop[herring][i][indx]) / WB * pop[herring][i][number];
				pop[herring][i][indy] = QR * pop[herring][i][indx];
			}
		}


		// collect eggs mutant
		for (i = 0; i < cohort_no[1]; i++) {
			if (pop[1][i][indy] / pop[1][i][indx] > QR && popIDcard[1][i][IDlen] > LM && pop[1][i][number] > 0) {

				double yearrep = KR * (pop[1][i][indy] - QR * pop[1][i][indx]) / WB * pop[1][i][number] / MUTINIT;
				int cohoindex = (int)round(popIDcard[1][i][IDmutindex]);

				LROcoho[0][cohoindex] += yearrep;
				gentime[0][cohoindex] += popIDcard[1][i][IDage] * yearrep;

				pop[1][i][indy] = QR * pop[1][i][indx];

			}
		}

		// cohorts
		int NewCohort = AddCohorts(pop, herring, SD + 1);

		for (i = 0; i < SD + 1; i++) {
			temp = gettemp(time + i) + getheatwave(time + i);

			if (i == 0)			pop[herring][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * 1.0 / (1 + exp(-(0 - SM))) * eggbuffer;// first reproduction day to make distribution add to 1
			else if (i == SD) 		pop[herring][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 - 1.0 / (1 + exp(-(SD - 1 - SM)))) * eggbuffer; // last reproduction day to make distribution add to 1
			else 				pop[herring][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 / (1 + exp(-(i - SM))) - 1.0 / (1 + exp(-(i - 1 - SM)))) * eggbuffer;

			popIDcard[herring][NewCohort + i][IDbtime] = time + i; // egg release time
			popIDcard[herring][NewCohort + i][IDbage] = (AF + EPSMIN + EPSP * phiarr(MEANTEMP, EAE, TR)); // egg and larvea duration

			popIDcard[herring][NewCohort + i][IDage] = -i - popIDcard[herring][NewCohort + i][IDbage]; // the age is zero as soon as individuals turn into juveniles
			pop[herring][NewCohort + i][indx] = WB;
			pop[herring][NewCohort + i][indy] = QJ * WB;

			popIDcard[herring][NewCohort + i][IDmass] = (1 + QR) * pop[herring][NewCohort + i][indx];
			popIDcard[herring][NewCohort + i][IDlen] = LAMBDA1 * pow(popIDcard[herring][NewCohort + i][IDmass], LAMBDA2);
			popIDcard[herring][NewCohort + i][IDcond] = pop[herring][NewCohort + i][indy] / pop[herring][NewCohort + i][indx];

			popIDcard[herring][NewCohort + i][IDmatage] = MAXAGE;
		}



		if (fmod(time, BIFPERIOD * 365) > MUTYEAR * 365 - TINY && fmod(time, BIFPERIOD * 365) < (MUTYEAR + MUTPERIOD - 1) * 365 + TINY) {
			// cohorts

			int NewCohort = AddCohorts(pop, 1, SD + 1);


			for (i = 0; i < SD + 1; i++) {
				temp = gettemp(time + i) + getheatwave(time + i);

				if (i == 0)			pop[1][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * 1.0 / (1 + exp(-(0 - SM))) * MUTINIT;// first reproduction day to make distribution add to 1
				else if (i == SD) 		pop[1][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 - 1.0 / (1 + exp(-(SD - 1 - SM)))) * MUTINIT; // last reproduction day to make distribution add to 1
				else 				pop[1][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 / (1 + exp(-(i - SM))) - 1.0 / (1 + exp(-(i - 1 - SM)))) * MUTINIT;

				popIDcard[1][NewCohort + i][IDbtime] = time + i; // egg release time
				popIDcard[1][NewCohort + i][IDbage] = (AF + EPSMIN + EPSP * phiarr(MEANTEMP, EAE, TR)); // egg and larvea duration

				popIDcard[1][NewCohort + i][IDage] = -i - popIDcard[1][NewCohort + i][IDbage]; // the age is zero as soon as individuals turn into juveniles
				pop[1][NewCohort + i][indx] = WB;
				pop[1][NewCohort + i][indy] = QJ * WB;

				popIDcard[1][NewCohort + i][IDmass] = (1 + QR) * pop[1][NewCohort + i][indx];
				popIDcard[1][NewCohort + i][IDlen] = LAMBDA1 * pow(popIDcard[1][NewCohort + i][IDmass], LAMBDA2);
				popIDcard[1][NewCohort + i][IDcond] = pop[1][NewCohort + i][indy] / pop[1][NewCohort + i][indx];

				popIDcard[1][NewCohort + i][IDmatage] = MAXAGE;

				popIDcard[1][NewCohort + i][IDinit] = eggbuffer;

				popIDcard[1][NewCohort + i][IDmutindex] = round(fmod((time) / 365, BIFPERIOD) - MUTYEAR);

			}

		}

		// IDinit of early mutants need to be set but cannot be done earlier
		for (i = 0; i < cohort_no[2]; i++) {
			if (popIDcard[2][i][IDinit] == -1) {
				popIDcard[2][i][IDinit] = eggbuffer;
			}
		}

	}

	if (fmod(time, 365) > 1 - TINY && fmod(time, 365) < 1 + TINY) { // late mutant

		temp = gettemp(time) + getheatwave(time);


		// collect eggs mutant
		for (i = 0; i < cohort_no[3]; i++) {
			if (pop[3][i][indy] / pop[3][i][indx] > QR && popIDcard[3][i][IDlen] > LM && pop[3][i][number] > 0) {

				double yearrep = KR * (pop[3][i][indy] - QR * pop[3][i][indx]) / WB * pop[3][i][number] / MUTINIT;
				int cohoindex = (int)round(popIDcard[3][i][IDmutindex]);

				LROcoho[2][cohoindex] += yearrep;
				gentime[2][cohoindex] += popIDcard[3][i][IDage] * yearrep;

				pop[3][i][indy] = QR * pop[3][i][indx];
			}
		}

		// new cohorts
		if (fmod(time -1, BIFPERIOD * 365) > MUTYEAR * 365 - TINY && fmod(time - 1, BIFPERIOD * 365) < (MUTYEAR + MUTPERIOD - 1) * 365 + TINY) {
			// cohorts
			int NewCohort = AddCohorts(pop, 3, SD + 1);

			for (i = 0; i < SD + 1; i++) {
				temp = gettemp(time + i) + getheatwave(time + i);

				if (i == 0)			pop[3][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * 1.0 / (1 + exp(-(0 - SM))) * MUTINIT;// first reproduction day to make distribution add to 1
				else if (i == SD) 		pop[3][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 - 1.0 / (1 + exp(-(SD - 1 - SM)))) * MUTINIT; // last reproduction day to make distribution add to 1
				else 				pop[3][NewCohort + i][number] = phiSS(MEANTEMP, EAS, EDS, TR, THS) * (1.0 / (1 + exp(-(i - SM))) - 1.0 / (1 + exp(-(i - 1 - SM)))) * MUTINIT;

				popIDcard[3][NewCohort + i][IDbtime] = time + i; // egg release time
				popIDcard[3][NewCohort + i][IDbage] = (AF + EPSMIN + EPSP * phiarr(MEANTEMP, EAE, TR)); // egg and larvea duration

				popIDcard[3][NewCohort + i][IDage] = -i - popIDcard[3][NewCohort + i][IDbage]; // the age is zero as soon as individuals turn into juveniles
				pop[3][NewCohort + i][indx] = WB;
				pop[3][NewCohort + i][indy] = QJ * WB;

				popIDcard[3][NewCohort + i][IDmass] = (1 + QR) * pop[3][NewCohort + i][indx];
				popIDcard[3][NewCohort + i][IDlen] = LAMBDA1 * pow(popIDcard[3][NewCohort + i][IDmass], LAMBDA2);
				popIDcard[3][NewCohort + i][IDcond] = pop[3][NewCohort + i][indy] / pop[3][NewCohort + i][indx];

				popIDcard[3][NewCohort + i][IDmatage] = MAXAGE;

				popIDcard[3][NewCohort + i][IDinit] = eggbuffer;

				popIDcard[3][NewCohort + i][IDmutindex] = round(fmod((time - 1) / 365, BIFPERIOD) - MUTYEAR);
			}
		}

	}



	return;
}



/*
 *==========================================================================
 *
 *			SPECIFICATION OF OUTPUT VARIABLES
 *
 *==========================================================================
 */

void	DefineOutput(double *env, population *pop, double *output)
  
  /*
   * The function header above should not be changed.
   *
   * Define below the values of the output variables in terms of the
   * population and environment statistics. These values have to be
   * returned to the main program in the array "output[]".
   *
   *     NB : -	The integer array "cohort_no[]" is globally available and
   *		denotes the number of internal cohorts present in each
   *		structured population. 
   */

{						/* OUTPUT VARIABLES:        */


	UpdateIDcards(env, pop);

	register int		i;	

	// environment
	output[0] = env[0];
	output[1] = env[1];
	output[2] = env[2];
	output[3] = env[3];
  
	// pop starts
	output[4] = 0;
	output[5] = 0;
	output[6] = 0;
	output[7] = MAXAGE;
	output[8] = 0;

	output[9] = gettemp(time) + getheatwave(time);

	output[10] = 0;
	output[11] = 0;
	output[12] = 0;

	output[13] = 0;
	output[14] = 0;
	output[15] = 0;

	output[16] = 0;
	output[17] = 0;
	output[18] = 0;

	for (int i = 0; i < MUTPERIOD; i++) {


		if (LROcoho[0][i] > 0 && gentime[0][i] > 0) output[10] += log(LROcoho[0][i]) / (gentime[0][i] / LROcoho[0][i]) * 365 / MUTPERIOD;
		if (LROcoho[1][i] > 0 && gentime[1][i] > 0) output[11] += log(LROcoho[1][i]) / (gentime[1][i] / LROcoho[1][i]) * 365 / MUTPERIOD;
		if (LROcoho[2][i] > 0 && gentime[2][i] > 0) output[12] += log(LROcoho[2][i]) / (gentime[2][i] / LROcoho[2][i]) * 365 / MUTPERIOD;

		if (LROcoho[0][i] > 0) output[13] += log(LROcoho[0][i]) / MUTPERIOD;
		if (LROcoho[1][i] > 0) output[14] += log(LROcoho[1][i]) / MUTPERIOD;
		if (LROcoho[2][i] > 0) output[15] += log(LROcoho[2][i]) / MUTPERIOD;

		if (LROcoho[0][i] > 0) output[16] = gentime[0][i] / LROcoho[0][i] / 365;
		if (LROcoho[1][i] > 0) output[17] = gentime[1][i] / LROcoho[1][i] / 365;
		if (LROcoho[2][i] > 0) output[18] = gentime[2][i] / LROcoho[2][i] / 365;
		
	}


	for(i=0; i<cohort_no[0]; i++){	
		if (popIDcard[herring][i][IDage] < -popIDcard[herring][i][IDbage]) continue;

		if (popIDcard[herring][i][IDage] < 0 ) output[4] += pop[0][i][number];
		else if (popIDcard[herring][i][IDlen] < LM) output[5] += pop[0][i][number];
		else output[6] += pop[0][i][number];

		if (popIDcard[herring][i][IDmatage] < output[7]) output[7] = popIDcard[herring][i][IDmatage];
		if (popIDcard[herring][i][IDlen] > output[8]) output[8] = popIDcard[herring][i][IDlen];
    }

	return;



}

void UpdateIDcards(double* env, population* pop)
{
	register int		i;


	
	// set cohort values
	for (int j = 0; j < 4; j++) {
		for (i = 0; i < cohort_no[j]; i++)
		{
			popIDcard[j][i][IDage] = time - popIDcard[j][i][IDbtime] - popIDcard[j][i][IDbage];
			popIDcard[j][i][IDmass] = (1 + QR) * pop[j][i][indx];
			popIDcard[j][i][IDlen] = LAMBDA1 * pow(popIDcard[j][i][IDmass], LAMBDA2);
			popIDcard[j][i][IDcond] = pop[j][i][indy] / pop[j][i][indx];

			if (popIDcard[j][i][IDmatage] > popIDcard[j][i][IDage] && popIDcard[j][i][IDlen] > LM) popIDcard[j][i][IDmatage] = popIDcard[j][i][IDage];
		}
	}


	return;
}



/*==========================================================================*/
