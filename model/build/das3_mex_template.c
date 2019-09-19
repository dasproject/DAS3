// Compilation settings
#define VERBOSE 0					// set this to 1 for debugging

// size of the model (some other size constants are in the header file)
#define MAXMUSCLES 200						// maximum number of muscles in the model
#define NSEG ??								// number of segments
#define MAXSTATES (2*NDOF+2*MAXMUSCLES)		// max number of system state variables
#define MAXMUSDOF 	11              		// max number of DOF between a muscle's origin and insertion 
#define MAXPOLTERMS 70              		// max number of polynomial terms in muscle geometry 
#define NAMLEN  	40              		// max length of names

// macro to extract all mass properties of a segment and store them as scalars in the parameters struct
// S: segment name, see header file.  NUM: index for the segments array that came out of Opensim
#define EXTRACTINERTIAL(S,NUM) {										\
	parameters.S##Mass = extract2(P,"segments",(NUM),"mass",1);		\
	parameters.S##CMx 	= extract2(P,"segments",(NUM),"mass_center",1);	\
	parameters.S##CMy 	= extract2(P,"segments",(NUM),"mass_center",2);	\
	parameters.S##CMz 	= extract2(P,"segments",(NUM),"mass_center",3);	\
	parameters.S##Ixx 	= extract2(P,"segments",(NUM),"inertia",1);		\
	parameters.S##Ixy 	= extract2(P,"segments",(NUM),"inertia",2);		\
	parameters.S##Ixz 	= extract2(P,"segments",(NUM),"inertia",3);		\
	parameters.S##Iyy 	= extract2(P,"segments",(NUM),"inertia",5);		\
	parameters.S##Iyz 	= extract2(P,"segments",(NUM),"inertia",6);		\
	parameters.S##Izz 	= extract2(P,"segments",(NUM),"inertia",9);		\
}

typedef struct {
    char   name[NAMLEN];    				// Name of muscle 
	double Lceopt; 							// Optimal length of CE (m)
    double pennopt;         				// Pennation angle at optimal length of CE, in radians
	double width;							// Width of CE force-length relationship relative to Lceopt
	double Fmax;							// Maximal isometric force of CE (N)
	double Vmax;							// Max. contraction velocity in Lceopt/s
	double Tact, Tdeact;					// Activation and deactivation time constants
	double gmax;							// Maximum eccentric force
	double SEEslack;						// Slack length of the SEE (m)
	double PEEslack;						// Slack length of the PEE, relative to Lceopt
	double mass;						    // Muscle mass, in kg
	double umax;							// Strain of SEE at Fmax load
	double krel;                            // Stiffness of PEE, force/Fmax at elongation of Width*Lceopt
	double HillA;							// Normalized Hill parameter a/Fmax for F-v relationship (usually 0.25)
	double Gmax;							// Maximum eccentric muscle force, relative to Fmax (usually 1.5)
	int    nmusdof;              			// Number of DOFs between origin and insertion 
    int    musdof[MAXMUSDOF];    			// List of DOFs between origin and insertion 
    int    npolterms;            			// Number of terms in polynomial 
    double  polcoeff[MAXPOLTERMS];        	// Polynomial coefficients 
    int    expon[MAXPOLTERMS][MAXMUSDOF];   // Polynomial exponents 
    int    crossGH;            				// 1 if the muscle cross GH, 0 otherwise 
	
    // the following parameters are only relevant if the muscle crosses GH
    int    nFpolterms[3];            		// number of terms in force vector polynomial 
    double  xFpolcoeff[MAXPOLTERMS];        // polynomial coefficients for the x direction of force 
    double  yFpolcoeff[MAXPOLTERMS];        // polynomial coefficients for the y direction of force 
    double  zFpolcoeff[MAXPOLTERMS];        // polynomial coefficients for the z direction of force 
    int    xFexpon[MAXPOLTERMS][MAXMUSDOF]; // polynomial exponents for the x direction of force 
    int    yFexpon[MAXPOLTERMS][MAXMUSDOF]; // polynomial exponents for the y direction of force 
    int    zFexpon[MAXPOLTERMS][MAXMUSDOF]; // polynomial exponents for the z direction of force 

	// the following parameters are derived from other parameters during initialization
	double c3;								// Continuity parameter for eccentric force-velocity relationship
	double kPEE;							// Stiffness parameter of PEE, relative to Fmax/Lceopt^2
	double kSEE;							// Stiffness parameter of SEE, in Fmax/m^2	
} muscleprop;

typedef struct {
	char Name[NAMLEN];     	// Name of joint (= name of DOF)
	double Lowerlim;       	// Lower joint limit
	double Upperlim;       	// Upper joint limit 
	double Midpoint;		// halfway between limits
	double K1;				// linear stiffness (small), in N/rad
	double K2;				// quadratic stiffness, outside range of motion, in N/rad^2
	double B;				// damping, in Ns/rad
}jointprop;

// static global variables are used to preserve model parameters betweem MEX function calls
static int initialized = 0;					// this will be set to 1 when initialization is complete
static param_struct parameters;				// contains model parameters that must be made available to Autolev multibody model
static muscleprop muscles[MAXMUSCLES];		// contains muscle parameters
static jointprop joints[NDOF];				// joint moment parameters
static int NMUS, ncontacts, NSTATES, nf, nstrainterms;		// model size parameters.  nf is number of elements in the model expressions f(x,xdot,u,M)
static mwSize NMUSmx, NSTATESmx;  // NMUS in new Matlab mex type to use 64-bit API
static double zeros[MAXSTATES];				// an array of zeros
static int nonzeromomentarms;			// number of nonzero moment arms (nonzeros in dL/dq), for memory allocation

// get pointer to a MEX input vector
double* getvec(const mxArray *prhs[], int argnum, int len) {
	int nrows, ncols;
	nrows = mxGetM(prhs[argnum-1]);
	ncols = mxGetN(prhs[argnum-1]);
	if (!mxIsDouble(prhs[argnum-1]) || mxIsComplex(prhs[argnum-1]) ) {
		printf("Input %d must be double\n", argnum);
		mexErrMsgTxt("Incorrect type for input.");
	}
	if ((nrows != len) || (ncols != 1) ) {
		printf("Input %d must be a %d x 1 column vector\n", argnum, len);
		mexErrMsgTxt("Incorrect size for input.");
	}
	return mxGetPr(prhs[argnum-1]);
}

//===================================================================================================
// normalize: normalize a vector to length 1, warn if length was not close enough to 1 to begin with
//===================================================================================================
void normalize(double *x, double *y, double *z) {
	double length;
	length = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z));
//	if (fabs(length - 1.0) > 1e-2) {
//		printf("warning: vector was not normalized: %f %f %f -> length %f\n", *x, *y, *z, length);		
//	}
	*x = *x / length;
	*y = *y / length;
	*z = *z / length;
}

//===================================================================================
// MusclePath: calculate muscle-tendon length and its derivatives w.r.t. joint angles
//===================================================================================
double MusclePath(muscleprop *mus, double q[NDOF], double dL_dq[NDOF], double dL_dqdq[NDOF][NDOF]) {

	// returns the length Lm, its gradient dLm/dq, and its Hessian dLm/dqdq for muscle m
	
	int i,j,k,m,kdof,mdof;
	double L, term, dterm;

	// initialize length and all derivatives to 0 
	L = 0.0;
	for (i=0; i<NDOF; i++) {
		dL_dq[i] = 0.0;    
		for (j=0; j<NDOF; j++)		
			dL_dqdq[i][j] = 0.0;       
	}		
	
	// add contributions from each polynomial term
    for(i=0; i < mus->npolterms; i++) {
	
		// add this term's contribution to the muscle length 
		term = mus->polcoeff[i];
		for(j=0; j < mus->nmusdof; j++) {
			mdof = mus->musdof[j];
			for(k=0; k < mus->expon[i][j]; k++)
				term = term * q[mdof];			// this creates polcoeff[i] * product of all q_j to the power expon[i][j] 
		}
		L = L + term;	

		// first derivatives of L with respect to all q's
		for (k=0; k < mus->nmusdof; k++) {
			kdof = mus->musdof[k];
			if ((mus->expon[i][k] > 0) && (q[kdof] != 0.0)) {		// derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
				dterm = mus->expon[i][k]*term/q[kdof];
				dL_dq[kdof] = dL_dq[kdof] + dterm;
			
				// second derivatives dL_dqk_dqm for m not equal to k
				for (m=0; m < k; m++) {									// no need for m > k because the dL_dqdq matrix is symmetric
					mdof = mus->musdof[m];
					if ((mus->expon[i][m] > 0) && (q[mdof] != 0.0))
						dL_dqdq[kdof][mdof] = dL_dqdq[kdof][mdof] + mus->expon[i][m]*dterm/q[mdof];
						dL_dqdq[mdof][kdof] = dL_dqdq[kdof][mdof] ;		// symmetry
				}
				
				// second derivatives dL_dqk_dqk
				if (mus->expon[i][k] > 1)
					dL_dqdq[kdof][kdof] = dL_dqdq[kdof][kdof] + (mus->expon[i][k]-1)*dterm/q[kdof];
					
			}
		}
    }
	
	return L;
}

// ==========================================================================================
// get_force_vector: calculates the muscle's GH force vector at the current state
// 						this is how much 1 N of muscle force would contribute to the GH force
// ==========================================================================================
void get_force_vector(muscleprop *m, double q[], double force_vec[]) {

	int i,j,k;
	double term, qnonzero;

	for (i=0; i<3 ; i++) force_vec[i] = 0.0;            // initialize the force vector to 0

	// x component
	for(i=0; i < m->nFpolterms[0]; i++) {
		// add this term's contribution to the force component
		term = m->xFpolcoeff[i];
		for(j=0; j < m->nmusdof; j++) {
			qnonzero = q[m->musdof[j]];
			for(k=0; k < m->xFexpon[i][j]; k++)
				term = term * qnonzero;	// this creates q to the power expon[][]
		}
		force_vec[0] = force_vec[0] + term;
	}
	
	// y component
	for(i=0; i < m->nFpolterms[1]; i++) {
		// add this term's contribution to the force component
		term = m->yFpolcoeff[i];
		for(j=0; j < m->nmusdof; j++) {
			qnonzero = q[m->musdof[j]];
			for(k=0; k < m->yFexpon[i][j]; k++)
				term = term * qnonzero;	// this creates q to the power expon[][]
		}
		force_vec[1] = force_vec[1] + term;
	}

	// z component	
	for(i=0; i < m->nFpolterms[2]; i++) {
		// add this term's contribution to the force component
		term = m->zFpolcoeff[i];
		for(j=0; j < m->nmusdof; j++) {
			qnonzero = q[m->musdof[j]];
			for(k=0; k < m->zFexpon[i][j]; k++)
				term = term * qnonzero;	// this creates q to the power expon[][]
		}
		force_vec[2] = force_vec[2] + term;
	}
	normalize(&force_vec[0],&force_vec[1],&force_vec[2]);
	return;
}

//===================================================================================================
// MuscleDynamics: the implicit muscle contraction dynamics f()=0, returns f(a,s,sdot,Lm) and its derivatives
// For pennated muscle model, the state variable s is Lce*cos(p)
// state variable s is dimensionless (normalized to Lceopt)
// imbalance f is dimensionless (normalized to Fmax)
//===================================================================================================
double MuscleDynamics(muscleprop *m,									// muscle parameters (input)
	double a, double s, double sdot, double Lm,							// the input variables
	double *df_da, double *df_ds, double *df_dsdot, double *df_dLm, 	// the gradients (output)													
	double *force, double *dforce_dLm, double *dforce_ds){				// muscle force (N) and derivatives (output)
	
	double Lce, cosp, dLce_ds, dcosp_ds, b;
	double Lcedot, dLcedot_dsdot, dLcedot_ds;
	double F1, dF1_dLce, dF1_ds;
	double f,x,k1;
	double F2, dF2_dLcedot, dF2_dsdot, dF2_ds, dF2_da;
	double F3, dF3_dLce, dF3_ds;
	double F4, dF4_ds, dF4_dLm;
	double F5, dF5_dsdot;
	double lambda, dlambda_da;		// factor for activation-dependent maximal shortening velocity
	double c;						// continuity parameter for eccentric force-velocity relationship
			
	// If there is pennation, compute Lce and cos(p) from state s using the constant volume constraint: Lce*sin(p) = Lceopt*sin(popt)
    // Lce is dimensionless (normalized to Lceopt)
    // If pennation is zero, we can't do this because volume is zero, and Lce is equal to s
    if (m->pennopt < 0.01) {
        cosp = 1;
        Lce = s;
        dLce_ds = 1;
        dcosp_ds = 0;
        }
    else {
        b = sin(m->pennopt);            
        Lce = sqrt(s*s + b*b);
        cosp = s/Lce;
        dLce_ds = cosp;
        dcosp_ds = b*b/Lce/Lce/Lce;
	}
	
	// Compute Lcedot and its derivatives with respect to sdot and s
	Lcedot = sdot*cosp;
	dLcedot_dsdot = cosp;
	dLcedot_ds = sdot*dcosp_ds;
		
	// F1 is the normalized isometric force-length relationship at max activation
	x = (Lce - 1.0)/m->width;
	F1 = exp(-x*x);
	dF1_dLce = -2.0*x*F1 / m->width;
	dF1_ds = dF1_dLce * dLce_ds;
		
	// F2 is the normalized force-velocity relationship
//	lambda = 0.5025 + 0.5341*a;				// lambda is the factor for activation dependence of Vmax
//	dlambda_da = 0.5341;
	 lambda = 1.0;
	 dlambda_da = 0.0;
	if (Lcedot < 0) {
		// concentric contraction
		x = lambda * m->Vmax - Lcedot/m->HillA;
		F2 = (lambda * m->Vmax + Lcedot)/x;
		dF2_dLcedot = (1.0 + F2/m->HillA)/x;
		dF2_da = -dlambda_da * m->Vmax * Lcedot * (1.0 + 1.0/m->HillA) / x / x;
	}
	else {
		// eccentric contraction
		c = lambda * m->c3;
		x = Lcedot + c;
		F2 = (m->gmax*Lcedot + c) / x;
		dF2_dLcedot = (m->gmax - F2)/x;
		dF2_da = dlambda_da * m->c3 * Lcedot * (1.0 - m->gmax)/x/x;
	}
	dF2_dsdot =  dF2_dLcedot * dLcedot_dsdot;
	dF2_ds = dF2_dLcedot * dLcedot_ds;
	
	// F3 is the PEE force-length relationship
	k1 = 10.0/m->Fmax*m->Lceopt;	// stiffness of the linear term is 10 N/m, convert to Fmax/Lceopt units	
	x = (Lce - m->PEEslack);		// elongation of PEE, relative to Lceopt
	F3 = k1*x;						// low stiffness linear term
	dF3_dLce = k1;
	if (x>0) {						// add quadratic term for positive elongation						
		F3 = F3 + m->kPEE*x*x;
		dF3_dLce = dF3_dLce + 2*m->kPEE*x;
	}
	dF3_ds = dF3_dLce * dLce_ds;
	
	//  F4 is the SEE force-length relationship
	k1 = 10.0/m->Fmax;			// stiffness of the linear term is 10 N/m, convert to Fmax/meter	
	x = Lm - s * m->Lceopt - m->SEEslack;			// elongation of SEE, in meters
	F4 = k1*x;										// low stiffness linear term
	dF4_ds = -k1*m->Lceopt;
	dF4_dLm = k1;
	if (x>0) {										// add quadratic term for positive deformation
		F4 = F4 + m->kSEE*x*x;
		dF4_ds = dF4_ds - 2 * m->kSEE * m->Lceopt * x;
		dF4_dLm = dF4_dLm + 2 * m->kSEE * x;
	}

	// F5 is viscous damping in the projected CE (0.001 of Fmax at 1 Lceopt/s) to ensure df/dLcedot is never zero
	// this is only really needed if we want to solve Lcedot explicitly for all possible muscle states (including a=0)
	F5 = .001*sdot;
	dF5_dsdot = .001;
		
	// Compute f, the force imbalance in the muscle contraction, and its derivatives
	f = F4 - (a*F1*F2 + F3)*cosp - F5;
	*df_da = -(F1*F2 + a*F1*dF2_da)*cosp;
	*df_ds = dF4_ds - (a*(dF1_ds*F2 + F1*dF2_ds) + dF3_ds)*cosp - (a*F1*F2 + F3)*dcosp_ds;
	*df_dsdot = -a*F1*dF2_dsdot*cosp - dF5_dsdot;
	*df_dLm = dF4_dLm;
	
	// Muscle force is the force in SEE
	*force = m->Fmax*F4;
	*dforce_dLm = m->Fmax * dF4_dLm;
	*dforce_ds = m->Fmax * dF4_ds;
			
	// Return the imbalance
	return f;
}

//=========================================================================
// extract: returns value of P.fieldname1(i1) (double)
//=========================================================================
double extract(const mxArray *P, char *fieldname1, unsigned int i1) {
	
	mxArray *field;

	if ( ((field = mxGetField(P,0,fieldname1))==NULL) || 
	    (mxGetNumberOfDimensions(field)!=2) ||
		!mxIsDouble(field)||mxIsComplex(field) ) {
			printf("field name: %s\n", fieldname1);
			mexErrMsgTxt("Valid field with this name not found in parameters structure.");
	}
	return mxGetPr(field)[i1-1];
}

//=========================================================================
// extract2: returns value of P.fieldname1{i1}.fieldname2(i2) (double) 
//=========================================================================
double extract2(const mxArray *P, char *fieldname1, unsigned int i1, char *fieldname2, unsigned int i2) {
	
	mxArray *field;

	if ( ((field = mxGetField(P,0,fieldname1))==NULL) || 
	    (mxGetNumberOfDimensions(field)!=2) || 
		!mxIsCell(field)||(i1==0)) {
			printf("field name: %s\n", fieldname1);
			mexErrMsgTxt("Initialize: Valid field with this name not found in parameters structure.");
	}
    
    field = mxGetCell(field, i1-1);
    
	return extract(field, fieldname2, i2);
}

// =========================================================================
// mexFunction: this is the actual MEX function interface
// =========================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	char *command;				// command string given as first input to the MEX function

	// get the first argument, make sure that it is a command string
	if (nrhs < 1)
		mexErrMsgTxt("At least one input argument is needed.");
    if (!mxIsChar(prhs[0])) 
		mexErrMsgTxt("First argument must be a string.");
 	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("First argument must be a row vector.");
	command = mxArrayToString(prhs[0]);
	if(command == NULL) 
	  mexErrMsgTxt("First argument was not a command string.");
	  
	// initialize?
	if (strcmp(command, "Initialize") == 0) {
		const mxArray *P;									// pointer to the parameter struct
		printf("*****************************************************\n");
		printf("*                      DAS3MEX                      *\n");  
		printf("*  (c) 2010-2012 Case Western Reserve University    *\n");
		printf("*****************************************************\n");
		printf("Initializing...\n");
		if (nrhs != 2)
			mexErrMsgTxt("Initialize: Two input arguments required.");
			
		// get a pointer P to the second argument and check that it is a 1x1 struct
		P = prhs[1];
		if ((mxGetM(P) != 1) || (mxGetN(P) != 1)) {
			mexErrMsgTxt("Dynamics: Incorrect size for parameters, must be 1 x 1.");
		}
		if (!mxIsStruct(P)) {
			mexErrMsgTxt("Initialize: Incorrect type for parameters, must be struct.");	
		}
		
		// from P, extract the parameters needed by the Autolev generated code, and store them in the C struct "parameters"
		// see header file for information about the fields in the parameters.struct
		
		// gravity
		parameters.GravAxisX = extract(P, "gravity", 1);
		parameters.GravAxisY = extract(P, "gravity", 2);
		parameters.GravAxisZ = extract(P, "gravity", 3);
		
		// scapula-thorax contact force model
		parameters.epscontact = extract(P, "scap_thorax_eps", 1); 	// epsilon for scapula-thorax contact model
		parameters.kcontact = extract(P, "scap_thorax_k", 1);		// stiffness k (N/m) in scapula-thorax contact model
		parameters.Ax  = extract(P, "thorax_radii", 1);				// radii of the thorax ellipsoid
		parameters.Ay  = extract(P, "thorax_radii", 2);
		parameters.Az  = extract(P, "thorax_radii", 3);
		parameters.Mx  = extract(P, "thorax_center", 1);			// center of the thorax ellipsoid
		parameters.My  = extract(P, "thorax_center", 2);
		parameters.Mz  = extract(P, "thorax_center", 3);
		parameters.TSx = extract(P, "TSprojection", 1);				// coordinates of TS projection in scapula frame
		parameters.TSy = extract(P, "TSprojection", 2);
		parameters.TSz = extract(P, "TSprojection", 3);
		parameters.AIx = extract(P, "AIprojection", 1);				// coordinates of AI projection in scapula frame
		parameters.AIy = extract(P, "AIprojection", 2);
		parameters.AIz = extract(P, "AIprojection", 3);
		
		// conoid force model
		parameters.epsconoid = extract(P, "conoid_eps", 1); 	// epsilon for conoid force-length model
		parameters.conK = extract(P, "conoid_stiffness", 1);	// stiffness of conoid
		parameters.conL = extract(P, "conoid_length", 1);		// length (m) of conoid ligament
		parameters.conOx = extract(P, "conoid_origin", 1);		// coordinates of conoid origin in clavicle frame
		parameters.conOy = extract(P, "conoid_origin", 2);
		parameters.conOz = extract(P, "conoid_origin", 3);
		parameters.conIx = extract(P, "conoid_insertion", 1); 	// coordinates of conoid insertion in scapula frame
		parameters.conIy = extract(P, "conoid_insertion", 2);
		parameters.conIz = extract(P, "conoid_insertion", 3);
		
						
		// Check that the nDofs in the model struct is the same as the constant NDOF we use internally here
		if (NDOF - extract(P, "nDofs", 1) != 0) {
			mexErrMsgTxt("Initialize: nDofs is not consistent with this version of the MEX function.");		
		}
		
		// Check that the Nsegments in the model struct is the same as NSEG we use internally here
		if (NSEG - extract(P, "nSegments", 1) != 0) {
			mexErrMsgTxt("Initialize: nSegments is not consistent with this version of the MEX function.");		
		}
		
        // from P, extract the joint positions and axes, and segment mass properties
				
		// from P, extract the muscles and their properties
		{
			mxArray *field;
			int i,j,k;
			NMUS = (int) extract(P, "nMus", 1);
			NMUSmx = (mwSize) extract(P, "nMus", 1);
			if (NMUS > MAXMUSCLES) {
				mexErrMsgTxt("Initialize: too many muscles.");
			}
			if ((field = mxGetField(P,0,"muscles") ) == NULL) {
				mexErrMsgTxt("Initialize: no muscles field found in model.");
			}
			if (!mxIsCell(field)) {
					mexErrMsgTxt("Initialize: muscles field is not a cell array.");
			}
			nonzeromomentarms = 0;
			for (i=0; i<NMUS; i++) {
				mxArray *cell, *fieldincell;
				// printf("Extracting fields from muscle %d\n", i);
				cell = mxGetCell(field, i);
				
				// extract muscle name
				if ((fieldincell = mxGetField(cell,0,"name") ) == NULL) {
					mexErrMsgTxt("Initialize: no name field found in muscle.");
				}
				mxGetString(fieldincell, muscles[i].name, NAMLEN);
				
				// extract muscle properties
				muscles[i].Lceopt 		= extract(cell, "lceopt" , 1);
				muscles[i].pennopt 		= extract(cell, "pennopt" , 1);
				muscles[i].Fmax 		= extract(cell, "fmax" , 1);
				muscles[i].Vmax 		= extract(cell, "vmax" , 1);
				muscles[i].Tact 		= extract(cell, "tact" , 1);
				muscles[i].Tdeact 		= extract(cell, "tdeact" , 1);
				muscles[i].SEEslack 	= extract(cell, "lslack" , 1);
				muscles[i].mass 		= extract(cell, "mass" , 1);
				muscles[i].PEEslack 	= extract(cell, "PEEslack" , 1);
				muscles[i].width 		= 0.56;		
				muscles[i].gmax 		= 1.5;
				muscles[i].umax		 	= 0.04;					// Thelen does the SEE differently, but this is not critical
				muscles[i].krel			= 1.0;					// Thelen does the PEE differently, but this is not critical
				muscles[i].HillA		= 0.25;					// Thelen (2003) also used this value, and we used it before

				// some other properties are derived:
				muscles[i].c3 = muscles[i].Vmax * muscles[i].HillA * (muscles[i].gmax - 1.0) / (muscles[i].HillA + 1.0);
				muscles[i].kSEE = 1.0/(muscles[i].umax*muscles[i].umax*muscles[i].SEEslack*muscles[i].SEEslack);
				muscles[i].kPEE = muscles[i].krel / (muscles[i].width * muscles[i].width);
				
				// polynomials
				muscles[i].nmusdof		= (int) extract(cell, "dof_count", 1);
				
				// calculate how many non-zero moment arms we have in the entire model
				nonzeromomentarms = nonzeromomentarms + muscles[i].nmusdof;

				if (muscles[i].nmusdof > MAXMUSDOF) {
					printf("Muscle number: %d\n", i+1);
					mexErrMsgTxt("Initialize: too many degrees of freedom in muscle path.");
				}
				for (j=0; j<muscles[i].nmusdof; j++) {
					muscles[i].musdof[j] = (int) extract(cell, "dof_indeces", j+1) - 1;	// our dof index starts at 0 in this C code
				}
				muscles[i].npolterms	= (int) extract(cell, "lparam_count", 1);
				if (muscles[i].npolterms > MAXPOLTERMS) {
					printf("Muscle number: %d\n", i+1);
					mexErrMsgTxt("Initialize: too many polynomial terms in muscle path.");
				}
						
				for (j=0; j<muscles[i].npolterms; j++) {
					muscles[i].polcoeff[j] = extract(cell, "lcoefs", j+1);
					for (k=0; k<muscles[i].nmusdof; k++) {
						muscles[i].expon[j][k] = (int) extract(cell, "lparams", muscles[i].npolterms*k+j+1);
					}
				}
				muscles[i].crossGH = (int) extract(cell, "crossesGH", 1);
				if (muscles[i].crossGH ==1) {
					muscles[i].nFpolterms[0]	= (int) extract(cell, "xparam_count", 1);
					if (muscles[i].nFpolterms[0] > MAXPOLTERMS) {
						printf("Muscle number: %d\n", i+1);
						mexErrMsgTxt("Initialize: too many polynomial terms in x direction of GH force vector.");
					}
					muscles[i].nFpolterms[1]	= (int) extract(cell, "yparam_count", 1);
					if (muscles[i].nFpolterms[1] > MAXPOLTERMS) {
						printf("Muscle number: %d\n", i+1);
						mexErrMsgTxt("Initialize: too many polynomial terms in y direction of GH force vector.");
					}
					muscles[i].nFpolterms[2]	= (int) extract(cell, "zparam_count", 1);
					if (muscles[i].nFpolterms[2] > MAXPOLTERMS) {
						printf("Muscle number: %d\n", i+1);
						mexErrMsgTxt("Initialize: too many polynomial terms in z direction of GH force vector.");
					}
							
					for (j=0; j<muscles[i].nFpolterms[0]; j++) {
						muscles[i].xFpolcoeff[j] = extract(cell, "xcoefs", j+1);
						for (k=0; k<muscles[i].nmusdof; k++) {
							muscles[i].xFexpon[j][k] = (int) extract(cell, "xparams", muscles[i].nFpolterms[0]*k+j+1);
						}
					}
					for (j=0; j<muscles[i].nFpolterms[1]; j++) {
						muscles[i].yFpolcoeff[j] = extract(cell, "ycoefs", j+1);
						for (k=0; k<muscles[i].nmusdof; k++) {
							muscles[i].yFexpon[j][k] = (int) extract(cell, "yparams", muscles[i].nFpolterms[1]*k+j+1);
						}
					}
					for (j=0; j<muscles[i].nFpolterms[2]; j++) {
						muscles[i].zFpolcoeff[j] = extract(cell, "zcoefs", j+1);
						for (k=0; k<muscles[i].nmusdof; k++) {
							muscles[i].zFexpon[j][k] = (int) extract(cell, "zparams", muscles[i].nFpolterms[2]*k+j+1);
						}
					}				
				}
			}
		}
				
		// from P, extract the passive joint properties
		// Same for all joints for now!
 		{
			int i;
			for (i=0; i<NDOF; i++) {
//				joints[i].K1 		= extract2(P,"dofs",i+1,"stiffness_K1",1);
//				joints[i].K2 		= extract2(P,"dofs",i+1,"stiffness_K2",1);
//				joints[i].B 		= extract2(P,"dofs",i+1,"damping_B",1);
//				joints[i].K1 			= 0.5;
//				joints[i].K2 			= 500;
//				joints[i].B 			= 0.1;
				joints[i].K1 			= 1.0;
				joints[i].K2 			= 5000;
				joints[i].B 			= 1.0;
				joints[i].Lowerlim  	= extract2(P,"dofs",i+1,"range",1);
				joints[i].Upperlim 		= extract2(P,"dofs",i+1,"range",2);
				joints[i].Midpoint 		= (joints[i].Lowerlim + joints[i].Upperlim)/2;
			}
		}
 						
		// set the initialized flag and compute some constants
		initialized = 1;
		NSTATES = 2*NDOF + 2*NMUS;
		NSTATESmx = 2*NDOF + 2*NMUSmx;
		if (NSTATES > MAXSTATES)
			mexErrMsgTxt("Initialize: too many states.");
		nf = 2*NDOF + 2*NMUS;					// number of elements in f(x,xdot,u,M)
				
		return;
	}
	else {
		// commands other than Initialize are done here
		
		// MEX function pointers to inputs from Matlab
		double *x, *xdot, *u, *Mextra, *M, *forF, *handcomF;
		double *f, *FGHtotal, *Fscap, *visdata;
		double *df_dx, *df_dxdot, *df_du, *MA;
		mwIndex *df_dx_irs, *df_dx_jcs;
		mwIndex *df_dxdot_irs,*df_dxdot_jcs;
		mwIndex *df_du_irs, *df_du_jcs;
		mwIndex *MA_irs, *MA_jcs;
		mwIndex Ndf_dx, Ndf_dxdot, Ndf_du, NMA;
		double *mforces;				// muscle forces
		double *moments;				// joint moments
		double *LCEopt, *limits;
		double *SEEslack, *MuscleMass;
		int i,j,k;
		double ang, angvel, d;
		double GHforce_vec[3];
		
		// Inout output variables for Autolev
		double *q, *qd, *qdd;
		double zero[NDOF];
		double dz_dq[NDOF][NDOF];
		double dz_dqd[NDOF][NDOF];
		double dz_dqdd[NDOF][NDOF];
		double F_GH[3];
		double F_SCAP[2][3],dist_SCAP[2],pos_SCAP[2][3];
		double vis[NVIS][12];
		double qTH[3], mTH[3], exF[2], handF[3];
		
		// Muscle variables
		double Lm[MAXMUSCLES];							// muscle+tendon length, based on skeleton kinematic state
		double dLm_dang[MAXMUSCLES][NDOF];				// derivatives of Lm with respect to joint angles
		double dLm_dangdang[MAXMUSCLES][NDOF][NDOF]; 	// second derivative of Lm w.r.t. joint angles
		double force[MAXMUSCLES];						// muscle forces
		double g[MAXMUSCLES];							// muscle force imbalance
		double dg_da[MAXMUSCLES], dg_dLce[MAXMUSCLES], dg_dLcedot[MAXMUSCLES], dg_dLm[MAXMUSCLES]; 	// derivatives of muscle imbalance
		double h[MAXMUSCLES];							// activation dynamics imbalance
		double dh_da[MAXMUSCLES], dh_dadot[MAXMUSCLES], dh_du[MAXMUSCLES];
		double dforce_dLm[MAXMUSCLES], dforce_dLce[MAXMUSCLES];
		
		// Joint variables
		double mom[NDOF];				// total actuation to be supplied to Autolev generated code
		double dmom_dang[NDOF][NDOF];
		double dmom_dangvel[NDOF];		// is a diagonal matrix (no coupling between joints), so no need for full storage
		double dmom_dLce[NDOF][MAXMUSCLES];
		
		int do_ghforce = 0;				// do we need to output GH force?
		int do_SEEslack = 0;
		int do_MuscleMass = 0;
		int do_scapulacontact = 0;
		int do_fscapula = 0;
		int do_GHF_only = 0;
		int do_externalangles = 0;
		int do_crossGH = 0;		
		int do_dynamics = 0;
		int do_LCEopt = 0;
		int do_vis = 0;
		int do_moments = 0;
		int do_muscleforces = 0;
		int do_momentarms = 0;
		int do_lengths = 0;
		int do_name = 0;
		int do_limits = 0;
		
		int do_derivatives 	= (nlhs > 1);						// will be true if user requested derivatives;
		if (strcmp(command,"LCEopt") == 0)
			do_LCEopt = 1;
		else if (strcmp(command,"Visualization") == 0)
			do_vis = 1;
		else if (strcmp(command,"SEEslack") == 0)
			do_SEEslack = 1;
		else if (strcmp(command,"MuscleMass") == 0)
			do_MuscleMass = 1;
		else if (strcmp(command,"Scapulacontact") == 0)
			do_scapulacontact = 1;
		else if (strcmp(command,"Jointmoments") == 0)
			do_moments = 1;
		else if (strcmp(command,"Muscleforces") == 0)
			do_muscleforces = 1;
		else if (strcmp(command,"Momentarms") == 0)
			do_momentarms = 1;
		else if (strcmp(command,"Musclelengths") == 0)
			do_lengths = 1;
		else if (strcmp(command,"Musclename") == 0)
			do_name = 1;
		else if (strcmp(command,"crossGH") == 0)
			do_crossGH = 1;
		else if (strcmp(command,"Limits") == 0)
			do_limits = 1;
		else if (strcmp(command,"Dynamics") == 0)
			do_dynamics = 1;
		else {
			printf("Command given was: %s\n", command);
			mexErrMsgTxt("Unknown command.");
		}
	
		// Give error message if model was not initialized
		if (!initialized) {
			printf("Command given: %s\n", command);
			mexErrMsgTxt("Model was not initialized.");
		}

	// Get the inputs we need, and some commands can be completed already
	
	if (do_LCEopt) {
		if (nrhs != 1) {
			mexErrMsgTxt("No inputs must be given for command LCEopt.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("Command LCEopt only has one output.");
		}
		// Create NMUS x 1 matrix for the LCEopt output
		plhs[0] = mxCreateDoubleMatrix(NMUSmx, 1, mxREAL);
		LCEopt = mxGetPr(plhs[0]);
		// Copy LCEopt from the muscle structs
		for (i=0; i<NMUS; i++) *LCEopt++ = muscles[i].Lceopt;
		return;									// and we are done!
	}
	else if (do_SEEslack) {
		if (nrhs != 1) {
			mexErrMsgTxt("No inputs must be given for command SEEslack.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("Command SEEslack only has one output.");
		}
		// Create NMUS x 1 matrix for the SEEslack output
		plhs[0] = mxCreateDoubleMatrix(NMUSmx, 1, mxREAL);
		SEEslack = mxGetPr(plhs[0]);
		// Copy SEEslack from the muscle structs
		for (i=0; i<NMUS; i++) *SEEslack++ = muscles[i].SEEslack;
		return;									// and we are done!
	}
	else if (do_MuscleMass) {
		if (nrhs != 1) {
			mexErrMsgTxt("No inputs must be given for command MuscleMass.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("Command MuscleMass only has one output.");
		}
		// Create NMUS x 1 matrix for the MuscleMass output
		plhs[0] = mxCreateDoubleMatrix(NMUSmx, 1, mxREAL);
		MuscleMass = mxGetPr(plhs[0]);
		// Copy MuscleMass from the muscle structs
		for (i=0; i<NMUS; i++) *MuscleMass++ = muscles[i].mass;
		return;									// and we are done!
	}
	else if (do_crossGH) {
		if (nrhs != 2) {
			mexErrMsgTxt("One input (muscle number) must be given for command 'crossGH'.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("Command 'crossGH' only has one output.");
		}
		i = (int) *mxGetPr(prhs[1]);
		if (i<1 || i>NMUS) {
			mexErrMsgTxt("Muscle number out of range in command 'crossGH'.");
		}				
		plhs[0] = mxCreateDoubleScalar(muscles[i-1].crossGH);
		return;
	}
	else if (do_limits) {		
		if (nrhs != 1) {
			mexErrMsgTxt("No inputs must be given for command Limits.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("Command Limits only has one output.");
		}
		// Create NDOF x 2 matrix for the output
		plhs[0] = mxCreateDoubleMatrix(2, NDOF, mxREAL);
		limits = mxGetPr(plhs[0]);
		// Copy limits from the joint structs
		for (i=0; i<NDOF; i++) {
			*limits++ = joints[i].Lowerlim;
			*limits++ = joints[i].Upperlim;
		}
		return;									// and we are done!
	}
	else if (do_name) {
		if (nrhs != 2) {
			mexErrMsgTxt("One input (muscle number) must be given with command 'Musclename'.");
		}
		if (nlhs > 1) {
			mexErrMsgTxt("Command 'Musclename' only has one output.");
		}
		i = (int) *mxGetPr(prhs[1]);
		if (i<1 || i>NMUS) {
			mexErrMsgTxt("Muscle number out of range in command 'Musclename'.");
		}				
		plhs[0] = mxCreateString(muscles[i-1].name);
		return;
	}
	else if (do_dynamics) {
		if (nrhs < 4) {
			mexErrMsgTxt("Three inputs (x,xdot,u) are needed for command 'Dynamics'.");
		}
		if (nrhs > 7) {
			mexErrMsgTxt("No more than six inputs (x,xdot,u,M,exF,handF) are allowed for command 'Dynamics'.");
		}
		// determine which optional outputs were requested with dynamics
		if (nlhs == 4) {
			do_derivatives = 1;
		}
		else if (nlhs == 5) {
			do_derivatives = 1;
			do_ghforce = 1;
		}
		else if (nlhs == 6) {
			do_derivatives = 1;
			do_ghforce = 1;
			do_fscapula = 1;
		}
		else if (nlhs == 7) {
			do_derivatives = 1;
			do_ghforce = 1;
			do_fscapula = 1;
			do_externalangles = 1;
		}
		else if (nlhs != 1) {
			mexErrMsgTxt("One, four, five, six, or seven outputs must be given for das3mex dynamics.");
		}
		x 		= getvec(prhs, 2, NSTATES);
		xdot 	= getvec(prhs, 3, NSTATES);
		u 		= getvec(prhs, 4, NMUS);
	}
	else {					// all other commands have one input (state x) and one output
		if (nrhs != 2) {
			printf("Command %s needs exactly one input (x).\n", command);
			mexErrMsgTxt("Incorrect input.");
		}
		if (nlhs > 1) {
			printf("Command %s has only one output.\n", command);
			mexErrMsgTxt("Incorrect output.");
		}
		x 		= getvec(prhs, 2,NSTATES);
		xdot	= zeros;
		u 		= zeros;
	}
		
	// Compute the muscle dynamics, and get muscle forces
	for(i=0; i<NMUS; i++) {	
		// Calculate imbalance of activation dynamics
		double rate = u[i]/muscles[i].Tact + (1-u[i])/muscles[i].Tdeact;		// rate constant
		h[i] = xdot[2*NDOF+NMUS+i] - rate*(u[i] - x[2*NDOF+NMUS+i]);
		if (do_derivatives) {
			dh_da[i] = rate;
			dh_dadot[i] = 1.0;
			dh_du[i] = -rate - (1/muscles[i].Tact - 1/muscles[i].Tdeact)*(u[i] - x[2*NDOF+NMUS+i]);
		}
			
		// Calculate muscle length Lm and derivatives dLm/dq from generalized coordinates in x
		Lm[i] = MusclePath(&muscles[i], x, dLm_dang[i], dLm_dangdang[i]);

		// Calculate muscle force imbalance, normalized to Fmax
		g[i] = MuscleDynamics(&muscles[i],
			x[2*NDOF+NMUS+i],		// active state of muscle i
			x[2*NDOF+i],			// Lce of muscle i
			xdot[2*NDOF+i],			// Lcedot
			Lm[i],					// muscle length
			&dg_da[i], &dg_dLce[i], &dg_dLcedot[i], &dg_dLm[i],
			&force[i], &dforce_dLm[i], &dforce_dLce[i]);				
	}

	// Does user only want muscle forces?
	if (do_muscleforces) {
		// Create NMUS x 1 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(NMUSmx, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);
		// Copy muscle forces from the force[] array that we already have
		for (i=0; i<NMUS; i++) {
			*mforces++ = force[i];
		}
		return;				// done!
	}

	// Does user only want muscle lengths?
	if (do_lengths) {
		// Create NMUS x 1 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(NMUSmx, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);
		// Copy lengths from the Lm[] array that we already have
		for (i=0; i<NMUS; i++) {
			*mforces++ = Lm[i];
		}
		return;				// done!
	}

	// Does user only want moment arms?
	if (do_momentarms) {
		// Create NMUS x NDOF sparse matrix for result
		plhs[0] = mxCreateSparse(NMUSmx, NDOF, nonzeromomentarms, 0);
		MA = mxGetPr(plhs[0]);
		MA_irs = mxGetIr(plhs[0]);
		MA_jcs = mxGetJc(plhs[0]);
		NMA = 0;
		
		// copy moment arms from dLm_dang and store them in sparse NMUS x NDOF matrix
		// must be done column-wise
		for (j=0; j<NDOF; j++) {			
			MA_jcs[j] = NMA;					// store element number where this column starts
			for (i=0; i<NMUS; i++) {
				if (dLm_dang[i][j] != 0) {
					MA_irs[NMA] = i;				// store row number of this matrix element
					MA[NMA] = -dLm_dang[i][j];		// store the value of this matrix element
					NMA++;		
				}
			}
		}
		MA_jcs[NDOF] = NMA;		// store final element number for the matrix
		return;					// done!
	}
	
	// Compute the joint moments
	for (i=0; i<NDOF; i++) {
		// initialize derivatives to zero
		if (do_derivatives) {
			for (j=0; j<NDOF; j++) dmom_dang[i][j] = 0.0;
			for (j=0; j<NMUS; j++) dmom_dLce[i][j] = 0.0;
			dmom_dangvel[i] = 0.0;
		}
		
		// start with passive joint moment
		ang = x[i];									// joint angle i is state variable i
		angvel = x[NDOF+i];							// and the corresponding angular velocity
		mom[i] = -joints[i].K1*(ang-joints[i].Midpoint) - joints[i].B*angvel;		// a small linear elastic moment at all angles, plus damping
		if (do_derivatives) {
			dmom_dang[i][i] = -joints[i].K1;				// and its derivatives
			dmom_dangvel[i] = -joints[i].B;
		}
		d = joints[i].Lowerlim - ang;				// are we below min angle?
		if (d > 0) {								// yes, add quadratic term
			mom[i] = mom[i] + joints[i].K2*d*d;	
			if (do_derivatives) {						// and its derivative
				dmom_dang[i][i] = dmom_dang[i][i] - 2*joints[i].K2*d;
			}
		}
		d = ang - joints[i].Upperlim;
		if (d > 0) {								// are we above max angle?
			mom[i] = mom[i] - joints[i].K2*d*d;			// yes, add quadratic term
			if (do_derivatives) {
				dmom_dang[i][i] = dmom_dang[i][i] - 2*joints[i].K2*d;
			}
		}
		
		// add the muscle moments
		for (j=0; j<NMUS; j++) {
			mom[i] = mom[i] - dLm_dang[j][i]*force[j];		// moment arm is -dLm/dang
			if (do_derivatives) {
				for (k=0; k<NDOF; k++) {
					dmom_dang[i][k] = dmom_dang[i][k] - dLm_dang[j][i]*dforce_dLm[j]*dLm_dang[j][k] - dLm_dangdang[j][i][k]*force[j];
				}
				dmom_dLce[i][j] = dmom_dLce[i][j] - dLm_dang[j][i]*dforce_dLce[j];
			}
		}
	}
	
	// Does user only want joint moments?
	if (do_moments) {
		// Create NDOF x 1 matrix for the function output
		plhs[0] = mxCreateDoubleMatrix(NDOF, 1, mxREAL);
		moments = mxGetPr(plhs[0]);
		// Copy moments data from the mom[] array that we already have
		for (i=0; i<NDOF; i++) {
			*moments++ = mom[i];
		}
		return;				// done!
	}
	
	// did the user supply the 5 extra moments as 4th input of MEX function?
	if (nrhs > 4) {
		M = getvec(prhs, 5, 5);
		// put the first three moments into mTH:
		mTH[0] = M[0];
		mTH[1] = M[1];
		mTH[2] = M[2];
		// and add the last two moments to the last two elements of mom (flexion and pronation moments):
		mom[NDOF-2] = mom[NDOF-2] + M[3];
		mom[NDOF-1] = mom[NDOF-1] + M[4];
	}
	else {
		mTH[0] = 0;
		mTH[1] = 0;
		mTH[2] = 0;	
	}

	// did the user supply the mobile arm support force as 5th input of MEX function?
	if (nrhs > 5) {
		forF = getvec(prhs, 6, 2);
		exF[0] = forF[0];
		exF[1] = forF[1];
	}
	else {
		exF[0] = 0;
		exF[1] = 0;
	}

	// did the user supply the hand force as 6th input of MEX function?
	if (nrhs > 6) {
		handcomF = getvec(prhs, 7, 3);
		handF[0] = handcomF[0];
		handF[1] = handcomF[1];
		handF[2] = handcomF[2];
	}
	else {
		handF[0] = 0;
		handF[1] = 0;
		handF[2] = 0;
	}
	
	// Call the C function that was generated by Autolev and cleaned by autolevclean.exe
	q = &x[0];
	qd = &x[NDOF];	
	qdd = &xdot[NDOF];
	das3_al(&parameters, q, qd, qdd, mTH, exF, handF, zero, dz_dq, dz_dqd, dz_dqdd, F_GH, F_SCAP, dist_SCAP, pos_SCAP, vis, qTH);
		
	// Does user only want visualization data?
	if (do_vis) {
		// Create NVIS * 12 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(NVIS, 12, mxREAL);
		visdata = mxGetPr(plhs[0]);

		// Copy visdata from the matrix that came out of Autolev
		for (i=0; i<12; i++) {
			// fill column i
			for (j=0; j<NVIS; j++) {
				*visdata++ = vis[j][i];
			}
		}
		return;			// done!
	}
	
	// Does user only want scapula contact?
	if (do_scapulacontact) {
		// Create 2 x 1 matrix for the return argument
		plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);
		for (i=0; i<2; i++) {
//			*mforces++ = dist_SCAP[i];
			*mforces++ = (pos_SCAP[i][0]/parameters.Ax)*(pos_SCAP[i][0]/parameters.Ax)+(pos_SCAP[i][1]/parameters.Ay)*(pos_SCAP[i][1]/parameters.Ay)+(pos_SCAP[i][2]/parameters.Az)*(pos_SCAP[i][2]/parameters.Az)-1;
		}
		return;				// done!
	}
	
	// Assemble the MEX function output for GH force
	if (do_ghforce) {
		// Create 3x1 matrix for the 5th return argument, the total GH reaction force
		plhs[4] = mxCreateDoubleMatrix(3, 1, mxREAL);
		FGHtotal = mxGetPr(plhs[4]);

		// first get the net force from Autolev
		for (j=0; j<3; j++) {
			FGHtotal[j] = F_GH[j];
		}
		
		// add muscle contributions
		for (i=0; i<NMUS; i++) {
			if (muscles[i].crossGH) {  
				// only include muscles that cross GH
				get_force_vector(&muscles[i], q, GHforce_vec);
				for (j=0; j<3; j++) {
					FGHtotal[j] += force[i] * GHforce_vec[j];
				}
			}
		}
	}
	
	// Does user want scapula-thorax contact forces?
	if (do_fscapula) {
		// Create 3x2 matrix for the 6th return argument, the scapula-thorax contact forces
		plhs[5] = mxCreateDoubleMatrix(3, 2, mxREAL);
		f = mxGetPr(plhs[5]);
		for (i=0; i<2; i++) {
			for (j=0; j<3; j++) {
				*f++ = F_SCAP[i][j];
			}
		}
	}
	
	// Does user want thorax-humerus angles?
	if (do_externalangles) {
		// Create 3x1 matrix for the 7th return argument, the thorax-humerus angles
		plhs[6] = mxCreateDoubleMatrix(3, 1, mxREAL);
		f = mxGetPr(plhs[6]);
		for (i=0; i<3; i++) {
			*f++ = qTH[i];
		}
	}

	// Create matrix for the 1st return argument (f)
	plhs[0] = mxCreateDoubleMatrix(NSTATESmx, 1, mxREAL);
	f = mxGetPr(plhs[0]);
		
	// the first NDOF rows of the implicit differential equation (IDE) are: qdot-dq/dt = 0
	for (i=0; i<NDOF; i++) f[i] = x[NDOF+i] - xdot[i];
	// the next NDOF rows of the IDE are the equations of motion: ZERO+MOM=0
	for (i=0; i<NDOF; i++) f[NDOF+i] = zero[i] + mom[i];
	// the next NMUS rows of the IDE are the muscle contraction dynamics
	for (i=0; i<NMUS; i++) f[2*NDOF+i] = g[i];
	// the final NMUS rows of the IDE are the muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0
	for (i=0; i<NMUS; i++) f[2*NDOF+NMUS+i] =  h[i];
		
	// Create matrix for the Jacobian results
	if (do_derivatives) {	
		// The sparse Jacobians have to be filled in column order, using Matlab sparse data structure
		//int Ndf_dx_allocated = 610+2*nonzeromomentarms;	// number of non-zero elements in df/dx
		mwSize Ndf_dx_allocated = NDOF*(1 + 2*NDOF + NMUSmx) + 3*NMUS + nonzeromomentarms; // number of non-zero elements in df/dx
		//int Ndf_dxdot_allocated = 332;					// number of non-zero elements in  df/dxdot
		mwSize Ndf_dxdot_allocated = NDOF*(NDOF+1) + 2*NMUSmx;	// number of non-zero elements in  df/dxdot
		double almostzero = 1e-10;

		// --------Second output of MEX function: Jacobian df/dx
		plhs[1] = mxCreateSparse(NSTATESmx ,NSTATESmx, Ndf_dx_allocated, 0);

		df_dx = mxGetPr(plhs[1]);
		df_dx_irs = mxGetIr(plhs[1]);
		df_dx_jcs = mxGetJc(plhs[1]);
		Ndf_dx = 0;

		// derivatives with respect to q, columns 1..NDOF of x
		for (i=0; i<NDOF; i++) {			
			df_dx_jcs[i] = Ndf_dx;				// store element number where this column starts
			
			// derivatives of ZERO+MOM with respect to q are in rows NDOF+1 to 2*NDOF
			for (j=0; j<NDOF; j++) {
				df_dx_irs[Ndf_dx] = NDOF+j;		// store row number of this matrix element
				df_dx[Ndf_dx] = dz_dq[j][i];	// store the value of this matrix element
				// add the contributions dmom/dq
				df_dx[Ndf_dx] = df_dx[Ndf_dx] + dmom_dang[j][i];
				if (fabs(df_dx[Ndf_dx]) > almostzero) Ndf_dx++;		// if zero, we can overwrite it with the next matrix element
			}
			
			// derivatives of muscle imbalance with respect to q are in rows 2*NDOF+1 to 2*NDOF+NMUS
			for (j=0; j<NMUS; j++) {
				if (dLm_dang[j][i] != 0) {						// this only exists if muscle j crosses joint i
					df_dx_irs[Ndf_dx] = 2*NDOF+j;				// store row number of this matrix element
					df_dx[Ndf_dx] = dg_dLm[j]*dLm_dang[j][i];	// store the value of this matrix element
					if (fabs(df_dx[Ndf_dx]) > almostzero) Ndf_dx++;		// if zero, we can overwrite it with the next matrix element
				}
			}
		}
		
		// derivatives with respect to qdot, columns NDOF+1 to 2*NDOF of x
		for (i=0; i<NDOF; i++) {
			df_dx_jcs[NDOF+i] = Ndf_dx;		// store element number where this column starts
			
			// derivatives of (qdot-dq/dt) with respect to qdot are diagonal, in rows 1 to NDOF
			df_dx_irs[Ndf_dx] = i;		// store row number of this matrix element
			df_dx[Ndf_dx] = 1.0;		// store the value of this matrix element
			Ndf_dx++;
			
			// derivatives of ZERO+MOM with respect to qdot are in rows NDOF+1 to 2*NDOF
			for (j=0; j<NDOF; j++) {
				df_dx_irs[Ndf_dx] = NDOF+j;			// store row number of this matrix element
				df_dx[Ndf_dx] = dz_dqd[j][i];		// store the value of this matrix element
				// add the contributions dmom/dqdot (remember dmom_dangvel is diagonal)
				if (j==i) df_dx[Ndf_dx] = df_dx[Ndf_dx] + dmom_dangvel[i];
				if (fabs(df_dx[Ndf_dx]) > almostzero) Ndf_dx++;
			}		
		}
		
		// derivatives with respect to Lce, columns 2*NDOF+1 to 2*NDOF+NMUS
		for (i=0; i<NMUS; i++) {			
			df_dx_jcs[2*NDOF+i] = Ndf_dx;		// store element number where this column starts
		
			// derivatives of ZERO+MOM with respect to Lce are in rows NDOF+1 to 2*NDOF
			for (j=0; j<NDOF; j++) {
				df_dx_irs[Ndf_dx] = NDOF+j;		// store row number of this matrix element
				df_dx[Ndf_dx] = dmom_dLce[j][i];
				if (fabs(df_dx[Ndf_dx]) > almostzero) Ndf_dx++;
			}
			
			// derivatives of muscle force balance with respect to Lce are diagonal, rows 2*NDOF+1 to 2*NDOF+NMUS
			df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
			df_dx[Ndf_dx] = dg_dLce[i];			// store the value of this matrix element
			if (fabs(df_dx[Ndf_dx]) > almostzero) Ndf_dx++;				
		}
		
		// derivatives with respect to Act, columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
		for (i=0; i<NMUS; i++) {			
			df_dx_jcs[2*NDOF+NMUS+i] = Ndf_dx;		// store element number where this column starts
		
			// derivatives of muscle force balance with respect to Act are diagonal, rows 2*NDOF+1 to 2*NDOF+NMUS
			df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
			df_dx[Ndf_dx] = dg_da[i];		// store the value of this matrix element
			if (df_dx[Ndf_dx] != 0.0) Ndf_dx++;				
		
			// derivatives of activation dynamics with respect to Act are diagonal, rows 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
			df_dx_irs[Ndf_dx] = 2*NDOF+NMUS+i;		// store row number of this matrix element
			df_dx[Ndf_dx] = dh_da[i];		// store the value of this matrix element
			if (fabs(df_dx[Ndf_dx]) > almostzero) Ndf_dx++;				
		}
		df_dx_jcs[NSTATES] = Ndf_dx;		// store final element number

		// --------Third output of MEX function: Jacobian df/dxdot
		plhs[2] = mxCreateSparse(NSTATESmx, NSTATESmx, Ndf_dxdot_allocated, 0);
		df_dxdot = mxGetPr(plhs[2]);
		df_dxdot_irs = mxGetIr(plhs[2]);
		df_dxdot_jcs = mxGetJc(plhs[2]);
		Ndf_dxdot = 0;

		// derivatives with respect to dq/dt, columns 1..NDOF of xdot
		for (i=0; i<NDOF; i++) {			
			df_dxdot_jcs[i] = Ndf_dxdot;			// store element number where this column starts

			// derivatives of (qdot-dq/dt) with respect to dq/dt are diagonal, in rows 1 to NDOF
			df_dxdot_irs[Ndf_dxdot] = i;		// store row number of this matrix element
			df_dxdot[Ndf_dxdot] = -1.0;			// store the value of this matrix element
			if (fabs(df_dxdot[Ndf_dxdot]) > almostzero) Ndf_dxdot++;
		}
		
		// derivatives with respect to dqdot/dt, columns NDOF+1 to 2*NDOF of xdot
		for (i=0; i<NDOF; i++) {
			df_dxdot_jcs[NDOF+i] = Ndf_dxdot;		// store element number where this column starts
			
			// derivatives of ZERO+MOM with respect to qdd are in rows NDOF to 2*NDOF
			for (j=0; j<NDOF; j++) {
				df_dxdot_irs[Ndf_dxdot] = NDOF+j;		// store row number of this matrix element
				df_dxdot[Ndf_dxdot] = dz_dqdd[j][i];	// store the value of this matrix element
				if (fabs(df_dxdot[Ndf_dxdot]) > almostzero) Ndf_dxdot++;
			}
		}
		
		// derivatives with respect to Lcedot, columns 2*NDOF+1 to 2*NDOF+NMUS
		for (i=0; i<NMUS; i++) {			
			df_dxdot_jcs[2*NDOF+i] = Ndf_dxdot;		// store element number where this column starts
						
			// derivatives of muscle force balance with respect to Lcedot are diagonal, rows 2*NDOF to 2*NDOF+NMUS
			df_dxdot_irs[Ndf_dxdot] = 2*NDOF+i;		// store row number of this matrix element
			df_dxdot[Ndf_dxdot] = dg_dLcedot[i];	// store the value of this matrix element
			if (fabs(df_dxdot[Ndf_dxdot]) > almostzero) Ndf_dxdot++;
		}
		
		// derivatives with respect to Actdot, columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
		for (i=0; i<NMUS; i++) {			
			df_dxdot_jcs[2*NDOF+NMUS+i] = Ndf_dxdot;		// store element number where this column starts
			
			// derivatives of activation dynamics with respect to Actdot are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
			df_dxdot_irs[Ndf_dxdot] = 2*NDOF+NMUS+i;		// store row number of this matrix element
			df_dxdot[Ndf_dxdot] = dh_dadot[i];				// store the value of this matrix element
			if (fabs(df_dxdot[Ndf_dxdot]) > almostzero) Ndf_dxdot++;				
		}
		df_dxdot_jcs[NSTATES] = Ndf_dxdot;		// store final element number
		
		// --------Fourth output of MEX function: Jacobian df/du
		plhs[3] = mxCreateSparse(NSTATESmx, NMUSmx, NMUSmx, 0);
		df_du = mxGetPr(plhs[3]);
		df_du_irs = mxGetIr(plhs[3]);
		df_du_jcs = mxGetJc(plhs[3]);
		Ndf_du = 0;
		
		// derivatives with respect to u of each muscle
		for (i=0; i<NMUS; i++) {			
			df_du_jcs[i] = Ndf_du;		// store element number where this column starts
		
			// derivatives of activation dynamics with respect to u are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
			df_du_irs[Ndf_du] = 2*NDOF+NMUS+i;		// store row number of this matrix element
			df_du[Ndf_du] = dh_du[i];				// store the value of this matrix element
			Ndf_du++;				
		}
		df_du_jcs[NMUS] = Ndf_du;		// store final element number
		
	}
	}
}
