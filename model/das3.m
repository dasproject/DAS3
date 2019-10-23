% This file documents the use of the MEX function das3.
% The actual function is the file das3.mexw32 or das3.mexw64

% The MEX function has several ways in which it can be used.

%==========================================================================
% Initialization

% load the "model" structure that contains the model parameters
load model_struct;

das3('Initialize',model);
% This needs to be done first.  During initialization, the model will be initialized from
% the "model" structure

%==========================================================================
% Dynamics
function [f, dfdx, dfdxdot, dfdu, FGH, FSCAP, qTH] = das3('Dynamics',x, xdot, u, M, exF, handF)
% This is to evaluate the model dynamics in the implicit form f(x, xdot, u) = 0.
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%	xdot	(nstates x 1) 	State derivatives
%	u		(nmus x 1) 	Muscle excitations

% Optional inputs
% 	M		(5 x 1)		Moments applied to the thorax-humerus YZY axes and the elbow flexion and supination axes
%   exF     (2 x 1)     Vertical force of amplitude exF(2) applied to the ulna at a distance of exF(1) (meters) from the elbow 
%						(to simulate a mobile arm support)
%	handF	(3 x 1) 	Force applied to the centre of mass of the hand 
%						(defined in the global frame: +X is laterally, +Y is upwards and +Z is posteriorly)

% Outputs
%	f		(nstates x 1) 			Dynamics imbalance

% Optional outputs
%	dfdx	(nstates x nstates sparse) 	Jacobian of f with respect to x
%	dfdxdot	(nstates x nstates sparse) 	Jacobian of f with respect to xdot
%	dfdu	(nstates x nmus sparse)	Jacobian of f with respect to u
%	FGH		(3 x 1)				3D GH contact force, acting on scapula, expressed in scapula reference frame
%   FSCAP 	(3 x 2)				3D Contact forces acting on TS and AI, expressed in thorax reference frame
%	qTH		(3 x 1)				angles between thorax and humerus (YZY sequence)

% Notes
% 	1. FGH is only correct when dynamics are satisfied, i.e. f is zero
%	2. The three Jacobians must always be requested together, or not at all.

%==========================================================================
% Extract mass parameters of the muscle elements
function MusMass = das3('MuscleMass')
%
% Outputs
%	MussMass	(nmus x 1)	Mass for all muscle elements (in kg)

%==========================================================================
% Extract LCEopt parameters of the muscle elements
function LCEopt = das3('LCEopt')
%
% Outputs
%	Lceopt	(nmus x 1)	Optimal fiber length parameters for all muscle elements (in meters)

%==========================================================================
% Extract SEEslack parameters of the muscle elements
function SEEslack = das3('SEEslack')
%
% Outputs
%	SEEslack (nmus x 1)	Slack length of series elastic element for all muscle elements (in meters)

%==========================================================================
% Extract dof limits
function limits = das3('Limits')
%
% Outputs
%   limits	(2 x ndof)	Lower and upper limits of the ndof dofs, in radians

%==========================================================================
% Compute stick figure coordinates
function vis = das3('Visualization', x)
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%
% Outputs
%   vis   	(nSegments x 12)	Position(3) and orientation(3x3) of nSegments segments
%
% Notes
%	1.	Output only depends on the first ndofs elements of x (the joint angles)

%==========================================================================
% Compute scapula contact
function F_contact = das3('Scapulacontact', x)
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%
% Outputs
%   F_contact (2 x 1)	Thorax surface equation solved for TS and AI

%==========================================================================
% Compute joint moments
function moments = das3('Jointmoments', x)
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%
% Outputs
%   moments (ndof x 1)			Joint moments (N m)

%==========================================================================
% Compute muscle forces
function forces = das3('Muscleforces', x)
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%
% Outputs
%   forces	(nmus x 1)			Muscle forces (N)

%==========================================================================
% Compute moment arms
function momentarms = das3('Momentarms', x)
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%
% Outputs
%   momentarms	(nmus x ndof sparse)	Moment arms of nmus muscle elements at ndof joints (meters)

%==========================================================================
% Compute muscle-tendon lengths
function lengths = das3('Musclelengths', x)
%
% Inputs
%	x		(nstates x 1) 	Model state, consists of ndof generalized coordinates (rad), ndof generalized velocities (rad/s),
%						nmus CE lengths (relative to LCEopt), and nmus muscle active states (between 0 and 1).
%
% Outputs
%   lengths	(nmus x 1)	Length of nmus muscle elements (meters)
%==========================================================================
% Find muscle name
function name = das3('Musclename', number)
%
% Inputs
%	number	(scalar) 	Number of a muscle element, must be in the range 1..nmus
%
% Outputs
%   name	(string)	Name of the muscle element, as define on the das3.bio file
%==========================================================================
% Find muscle name
function crossGH_flag = das3('crossGH', number)
%
% Inputs
%	number	(scalar) 	Number of a muscle element, must be in the range 1..nmus
%
% Outputs
%   crossGH_flag	(scalar)	1 if the muscle crosses GH, 0 if it doesn't