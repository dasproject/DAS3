% das3.al
% This is an input file for Autolev to generate dynamic equations of motion for the DAS3 model.

AUTOZ ON
OVERWRITE ALL
PAUSE 0
UnitSystem  kg,meter,sec

MotionVariables' SCy''
MotionVariables' SCz''
MotionVariables' SCx''
MotionVariables' ACy''
MotionVariables' ACz''
MotionVariables' ACx''
MotionVariables' GHy''
MotionVariables' GHz''
MotionVariables' GHyy''
MotionVariables' ELx''
MotionVariables' PSy''

%------------------------------------
% Newtonian frame
%------------------------------------
Newtonian ground

%------------------------------------
% Segment thorax
%------------------------------------
Body thorax
Constants par__thoraxMass
Mass thorax = par__thoraxMass
Constants par__thoraxIxx, par__thoraxIyy, par__thoraxIzz, par__thoraxIxy, par__thoraxIyz, par__thoraxIxz
Inertia thorax, par__thoraxIxx, par__thoraxIyy, par__thoraxIzz, par__thoraxIxy, par__thoraxIyz, par__thoraxIxz

Frames thorax1, thorax2
Constants par__baseXrot, par__baseYrot, par__baseZrot
Simprot(ground, thorax1, 1, par__baseXrot)
Simprot(thorax1, thorax2, 2, par__baseYrot)
Simprot(thorax2, thorax,  3, par__baseZrot);	% orientation of thorax

Point thoraxJoint
Constants par__baseX, par__baseY, par__baseZ
P_groundO_thoraxJoint> = Vector(ground, par__baseX, par__baseY, par__baseZ)	% location in parent
V_thoraxJoint_ground> = 0> % velocity of joint
A_thoraxJoint_ground> = 0> % acceleration of joint

Constants par__thoraxCMx, par__thoraxCMy, par__thoraxCMz
P_thoraxJoint_thoraxO> = Vector(thorax, par__thoraxCMx, par__thoraxCMy, par__thoraxCMz)  % location of CM
V_thoraxO_ground> = 0> % velocity of CM
A_thoraxO_ground> = 0> % acceleration of CM
W_thorax_ground> = 0>
ALF_thorax_ground> = 0>

%------------------------------------
% Segment clavicle1
%------------------------------------
Frame clavicle1

SCy_dot = SCy'
SCy_dotdot = SCy''

Constants par__SCyAxisX, par__SCyAxisY, par__SCyAxisZ
SCyAxis> = Vector(thorax, par__SCyAxisX, par__SCyAxisY, par__SCyAxisZ)
Simprot(thorax,clavicle1,SCyAxis>,SCy)

W_clavicle1_ground> = W_thorax_ground> + SCy_dot * SCyAxis>     % generates angular velocity of child body
ALF_clavicle1_ground> = ALF_thorax_ground> + SCy_dotdot * SCyAxis>     % generates angular acceleration of child body

Point clavicle1Joint
Constants par__sc1X, par__sc1Y, par__sc1Z
P_thoraxJoint_clavicle1Joint> 	= Vector(thorax, par__sc1X, par__sc1Y, par__sc1Z)	% location in parent
V2pts(ground, thorax, thoraxJoint, clavicle1Joint) % velocity of joint
A2pts(ground, thorax, thoraxJoint, clavicle1Joint) % acceleration of joint

%------------------------------------
% Segment clavicle2
%------------------------------------
Frame clavicle2

SCz_dot = SCz'
SCz_dotdot = SCz''

Constants par__SCzAxisX, par__SCzAxisY, par__SCzAxisZ
SCzAxis> = Vector(clavicle1, par__SCzAxisX, par__SCzAxisY, par__SCzAxisZ)
Simprot(clavicle1,clavicle2,SCzAxis>,SCz)

W_clavicle2_ground> = W_clavicle1_ground> + SCz_dot * SCzAxis>     % generates angular velocity of child body
ALF_clavicle2_ground> = ALF_clavicle1_ground> + SCz_dotdot * SCzAxis>     % generates angular acceleration of child body

Point clavicle2Joint
Constants par__sc2X, par__sc2Y, par__sc2Z
P_clavicle1Joint_clavicle2Joint> 	= Vector(clavicle1, par__sc2X, par__sc2Y, par__sc2Z)	% location in parent
V2pts(ground, clavicle1, clavicle1Joint, clavicle2Joint) % velocity of joint
A2pts(ground, clavicle1, clavicle1Joint, clavicle2Joint) % acceleration of joint

%------------------------------------
% Segment clavicler
%------------------------------------
Body clavicler
Constants par__claviclerMass
Mass clavicler = par__claviclerMass
Constants par__claviclerIxx, par__claviclerIyy, par__claviclerIzz, par__claviclerIxy, par__claviclerIyz, par__claviclerIxz
Inertia clavicler, par__claviclerIxx, par__claviclerIyy, par__claviclerIzz, par__claviclerIxy, par__claviclerIyz, par__claviclerIxz

SCx_dot = SCx'
SCx_dotdot = SCx''

Constants par__SCxAxisX, par__SCxAxisY, par__SCxAxisZ
SCxAxis> = Vector(clavicle2, par__SCxAxisX, par__SCxAxisY, par__SCxAxisZ)
Simprot(clavicle2,clavicler,SCxAxis>,SCx)

W_clavicler_ground> = W_clavicle2_ground> + SCx_dot * SCxAxis>     % generates angular velocity of child body
ALF_clavicler_ground> = ALF_clavicle2_ground> + SCx_dotdot * SCxAxis>     % generates angular acceleration of child body

Point claviclerJoint
Constants par__sc3X, par__sc3Y, par__sc3Z
P_clavicle2Joint_claviclerJoint> 	= Vector(clavicle2, par__sc3X, par__sc3Y, par__sc3Z)	% location in parent
V2pts(ground, clavicle2, clavicle2Joint, claviclerJoint) % velocity of joint
A2pts(ground, clavicle2, clavicle2Joint, claviclerJoint) % acceleration of joint

Constants par__claviclerCMx, par__claviclerCMy, par__claviclerCMz
P_claviclerJoint_claviclerO> = Vector(clavicler, par__claviclerCMx, par__claviclerCMy, par__claviclerCMz)  % location of CM
V2pts(ground, clavicler, claviclerJoint, claviclerO) % velocity of CM
A2pts(ground, clavicler, claviclerJoint, claviclerO) % acceleration of CM

%------------------------------------
% Segment scapula1
%------------------------------------
Frame scapula1

ACy_dot = ACy'
ACy_dotdot = ACy''

Constants par__ACyAxisX, par__ACyAxisY, par__ACyAxisZ
ACyAxis> = Vector(clavicler, par__ACyAxisX, par__ACyAxisY, par__ACyAxisZ)
Simprot(clavicler,scapula1,ACyAxis>,ACy)

W_scapula1_ground> = W_clavicler_ground> + ACy_dot * ACyAxis>     % generates angular velocity of child body
ALF_scapula1_ground> = ALF_clavicler_ground> + ACy_dotdot * ACyAxis>     % generates angular acceleration of child body

Point scapula1Joint
Constants par__ac1X, par__ac1Y, par__ac1Z
P_claviclerJoint_scapula1Joint> 	= Vector(clavicler, par__ac1X, par__ac1Y, par__ac1Z)	% location in parent
V2pts(ground, clavicler, claviclerJoint, scapula1Joint) % velocity of joint
A2pts(ground, clavicler, claviclerJoint, scapula1Joint) % acceleration of joint

%------------------------------------
% Segment scapula2
%------------------------------------
Frame scapula2

ACz_dot = ACz'
ACz_dotdot = ACz''

Constants par__ACzAxisX, par__ACzAxisY, par__ACzAxisZ
ACzAxis> = Vector(scapula1, par__ACzAxisX, par__ACzAxisY, par__ACzAxisZ)
Simprot(scapula1,scapula2,ACzAxis>,ACz)

W_scapula2_ground> = W_scapula1_ground> + ACz_dot * ACzAxis>     % generates angular velocity of child body
ALF_scapula2_ground> = ALF_scapula1_ground> + ACz_dotdot * ACzAxis>     % generates angular acceleration of child body

Point scapula2Joint
Constants par__ac2X, par__ac2Y, par__ac2Z
P_scapula1Joint_scapula2Joint> 	= Vector(scapula1, par__ac2X, par__ac2Y, par__ac2Z)	% location in parent
V2pts(ground, scapula1, scapula1Joint, scapula2Joint) % velocity of joint
A2pts(ground, scapula1, scapula1Joint, scapula2Joint) % acceleration of joint

%------------------------------------
% Segment scapular
%------------------------------------
Body scapular
Constants par__scapularMass
Mass scapular = par__scapularMass
Constants par__scapularIxx, par__scapularIyy, par__scapularIzz, par__scapularIxy, par__scapularIyz, par__scapularIxz
Inertia scapular, par__scapularIxx, par__scapularIyy, par__scapularIzz, par__scapularIxy, par__scapularIyz, par__scapularIxz

ACx_dot = ACx'
ACx_dotdot = ACx''

Constants par__ACxAxisX, par__ACxAxisY, par__ACxAxisZ
ACxAxis> = Vector(scapula2, par__ACxAxisX, par__ACxAxisY, par__ACxAxisZ)
Simprot(scapula2,scapular,ACxAxis>,ACx)

W_scapular_ground> = W_scapula2_ground> + ACx_dot * ACxAxis>     % generates angular velocity of child body
ALF_scapular_ground> = ALF_scapula2_ground> + ACx_dotdot * ACxAxis>     % generates angular acceleration of child body

Point scapularJoint
Constants par__ac3X, par__ac3Y, par__ac3Z
P_scapula2Joint_scapularJoint> 	= Vector(scapula2, par__ac3X, par__ac3Y, par__ac3Z)	% location in parent
V2pts(ground, scapula2, scapula2Joint, scapularJoint) % velocity of joint
A2pts(ground, scapula2, scapula2Joint, scapularJoint) % acceleration of joint

Constants par__scapularCMx, par__scapularCMy, par__scapularCMz
P_scapularJoint_scapularO> = Vector(scapular, par__scapularCMx, par__scapularCMy, par__scapularCMz)  % location of CM
V2pts(ground, scapular, scapularJoint, scapularO) % velocity of CM
A2pts(ground, scapular, scapularJoint, scapularO) % acceleration of CM

%------------------------------------
% Segment humerus1
%------------------------------------
Frame humerus1

GHy_dot = GHy'
GHy_dotdot = GHy''

Constants par__GHyAxisX, par__GHyAxisY, par__GHyAxisZ
GHyAxis> = Vector(scapular, par__GHyAxisX, par__GHyAxisY, par__GHyAxisZ)
Simprot(scapular,humerus1,GHyAxis>,GHy)

W_humerus1_ground> = W_scapular_ground> + GHy_dot * GHyAxis>     % generates angular velocity of child body
ALF_humerus1_ground> = ALF_scapular_ground> + GHy_dotdot * GHyAxis>     % generates angular acceleration of child body

Point humerus1Joint
Constants par__gh1X, par__gh1Y, par__gh1Z
P_scapularJoint_humerus1Joint> 	= Vector(scapular, par__gh1X, par__gh1Y, par__gh1Z)	% location in parent
V2pts(ground, scapular, scapularJoint, humerus1Joint) % velocity of joint
A2pts(ground, scapular, scapularJoint, humerus1Joint) % acceleration of joint

%------------------------------------
% Segment humerus2
%------------------------------------
Frame humerus2

GHz_dot = GHz'
GHz_dotdot = GHz''

Constants par__GHzAxisX, par__GHzAxisY, par__GHzAxisZ
GHzAxis> = Vector(humerus1, par__GHzAxisX, par__GHzAxisY, par__GHzAxisZ)
Simprot(humerus1,humerus2,GHzAxis>,GHz)

W_humerus2_ground> = W_humerus1_ground> + GHz_dot * GHzAxis>     % generates angular velocity of child body
ALF_humerus2_ground> = ALF_humerus1_ground> + GHz_dotdot * GHzAxis>     % generates angular acceleration of child body

Point humerus2Joint
Constants par__gh2X, par__gh2Y, par__gh2Z
P_humerus1Joint_humerus2Joint> 	= Vector(humerus1, par__gh2X, par__gh2Y, par__gh2Z)	% location in parent
V2pts(ground, humerus1, humerus1Joint, humerus2Joint) % velocity of joint
A2pts(ground, humerus1, humerus1Joint, humerus2Joint) % acceleration of joint

%------------------------------------
% Segment humerusr
%------------------------------------
Body humerusr
Constants par__humerusrMass
Mass humerusr = par__humerusrMass
Constants par__humerusrIxx, par__humerusrIyy, par__humerusrIzz, par__humerusrIxy, par__humerusrIyz, par__humerusrIxz
Inertia humerusr, par__humerusrIxx, par__humerusrIyy, par__humerusrIzz, par__humerusrIxy, par__humerusrIyz, par__humerusrIxz

GHyy_dot = GHyy'
GHyy_dotdot = GHyy''

Constants par__GHyyAxisX, par__GHyyAxisY, par__GHyyAxisZ
GHyyAxis> = Vector(humerus2, par__GHyyAxisX, par__GHyyAxisY, par__GHyyAxisZ)
Simprot(humerus2,humerusr,GHyyAxis>,GHyy)

W_humerusr_ground> = W_humerus2_ground> + GHyy_dot * GHyyAxis>     % generates angular velocity of child body
ALF_humerusr_ground> = ALF_humerus2_ground> + GHyy_dotdot * GHyyAxis>     % generates angular acceleration of child body

Point humerusrJoint
Constants par__gh3X, par__gh3Y, par__gh3Z
P_humerus2Joint_humerusrJoint> 	= Vector(humerus2, par__gh3X, par__gh3Y, par__gh3Z)	% location in parent
V2pts(ground, humerus2, humerus2Joint, humerusrJoint) % velocity of joint
A2pts(ground, humerus2, humerus2Joint, humerusrJoint) % acceleration of joint

Constants par__humerusrCMx, par__humerusrCMy, par__humerusrCMz
P_humerusrJoint_humerusrO> = Vector(humerusr, par__humerusrCMx, par__humerusrCMy, par__humerusrCMz)  % location of CM
V2pts(ground, humerusr, humerusrJoint, humerusrO) % velocity of CM
A2pts(ground, humerusr, humerusrJoint, humerusrO) % acceleration of CM

%------------------------------------
% Segment ulnar
%------------------------------------
Body ulnar
Constants par__ulnarMass
Mass ulnar = par__ulnarMass
Constants par__ulnarIxx, par__ulnarIyy, par__ulnarIzz, par__ulnarIxy, par__ulnarIyz, par__ulnarIxz
Inertia ulnar, par__ulnarIxx, par__ulnarIyy, par__ulnarIzz, par__ulnarIxy, par__ulnarIyz, par__ulnarIxz

ELx_dot = ELx'
ELx_dotdot = ELx''

Constants par__ELxAxisX, par__ELxAxisY, par__ELxAxisZ
ELxAxis> = Vector(humerusr, par__ELxAxisX, par__ELxAxisY, par__ELxAxisZ)
Simprot(humerusr,ulnar,ELxAxis>,ELx)

W_ulnar_ground> = W_humerusr_ground> + ELx_dot * ELxAxis>     % generates angular velocity of child body
ALF_ulnar_ground> = ALF_humerusr_ground> + ELx_dotdot * ELxAxis>     % generates angular acceleration of child body

Point ulnarJoint
Constants par__huX, par__huY, par__huZ
P_humerusrJoint_ulnarJoint> 	= Vector(humerusr, par__huX, par__huY, par__huZ)	% location in parent
V2pts(ground, humerusr, humerusrJoint, ulnarJoint) % velocity of joint
A2pts(ground, humerusr, humerusrJoint, ulnarJoint) % acceleration of joint

Constants par__ulnarCMx, par__ulnarCMy, par__ulnarCMz
P_ulnarJoint_ulnarO> = Vector(ulnar, par__ulnarCMx, par__ulnarCMy, par__ulnarCMz)  % location of CM
V2pts(ground, ulnar, ulnarJoint, ulnarO) % velocity of CM
A2pts(ground, ulnar, ulnarJoint, ulnarO) % acceleration of CM

%------------------------------------
% Segment radiusr
%------------------------------------
Body radiusr
Constants par__radiusrMass
Mass radiusr = par__radiusrMass
Constants par__radiusrIxx, par__radiusrIyy, par__radiusrIzz, par__radiusrIxy, par__radiusrIyz, par__radiusrIxz
Inertia radiusr, par__radiusrIxx, par__radiusrIyy, par__radiusrIzz, par__radiusrIxy, par__radiusrIyz, par__radiusrIxz

PSy_dot = PSy'
PSy_dotdot = PSy''

Constants par__PSyAxisX, par__PSyAxisY, par__PSyAxisZ
PSyAxis> = Vector(ulnar, par__PSyAxisX, par__PSyAxisY, par__PSyAxisZ)
Simprot(ulnar,radiusr,PSyAxis>,PSy)

W_radiusr_ground> = W_ulnar_ground> + PSy_dot * PSyAxis>     % generates angular velocity of child body
ALF_radiusr_ground> = ALF_ulnar_ground> + PSy_dotdot * PSyAxis>     % generates angular acceleration of child body

Point radiusrJoint
Constants par__urX, par__urY, par__urZ
P_ulnarJoint_radiusrJoint> 	= Vector(ulnar, par__urX, par__urY, par__urZ)	% location in parent
V2pts(ground, ulnar, ulnarJoint, radiusrJoint) % velocity of joint
A2pts(ground, ulnar, ulnarJoint, radiusrJoint) % acceleration of joint

Constants par__radiusrCMx, par__radiusrCMy, par__radiusrCMz
P_radiusrJoint_radiusrO> = Vector(radiusr, par__radiusrCMx, par__radiusrCMy, par__radiusrCMz)  % location of CM
V2pts(ground, radiusr, radiusrJoint, radiusrO) % velocity of CM
A2pts(ground, radiusr, radiusrJoint, radiusrO) % acceleration of CM

%------------------------------------
% Segment handr
%------------------------------------
Body handr
Constants par__handrMass
Mass handr = par__handrMass
Constants par__handrIxx, par__handrIyy, par__handrIzz, par__handrIxy, par__handrIyz, par__handrIxz
Inertia handr, par__handrIxx, par__handrIyy, par__handrIzz, par__handrIxy, par__handrIyz, par__handrIxz

Frames handr1, handr2
Constants par__rcXrot, par__rcYrot, par__rcZrot
Simprot(radiusr, handr1, 1, par__rcXrot)
Simprot(handr1, handr2, 2, par__rcYrot)
Simprot(handr2, handr,  3, par__rcZrot);	% orientation of handr

Point handrJoint
Constants par__rcX, par__rcY, par__rcZ
P_radiusrJoint_handrJoint> = Vector(radiusr, par__rcX, par__rcY, par__rcZ)	% location in parent
V2pts(ground, radiusr, radiusrJoint, handrJoint) % velocity of joint
A2pts(ground, radiusr, radiusrJoint, handrJoint) % acceleration of joint

Constants par__handrCMx, par__handrCMy, par__handrCMz
P_handrJoint_handrO> = Vector(handr, par__handrCMx, par__handrCMy, par__handrCMz)  % location of CM
V2pts(ground, handr, handrJoint, handrO) % velocity of CM
A2pts(ground, handr, handrJoint, handrO) % acceleration of CM

%------------------------------------
% Apply Gravity
%------------------------------------
Constants par__GravAxisX, par__GravAxisY, par__GravAxisZ
GravityAxis> = Vector(ground, par__GravAxisX, par__GravAxisY, par__GravAxisZ)
Gravity( GravityAxis> )

%--------------------------------------------------------------------
% Apply force on the forearm (ulna) to simulate a mobile arm support
%--------------------------------------------------------------------
Points ForF
Variables  distF, ampF
P_ulnarJoint_ForF> = Vector(ulnar,0,-distF,0)	
V2pts(ground, ulnar, ulnarJoint, ForF)
Force_ForF> = ampF*ground2> 

%--------------------------------------------------------------------
% Apply force to the CoM of the hand
%--------------------------------------------------------------------
Variables  handFx, handFy, handFz
Force_handrO> = handFx*ground1> + handFy*ground2> + handFz*ground3> 

%--------------------------------------------------------------------
% Calculate and apply conoid force
%--------------------------------------------------------------------
Points conO, conI							% Conoid origin (on the clavicle) and insertion (on the scapula)
Constants par__epsconoid, par__conK, par__conL				% parameters for the conoid ligament model
Constants par__conOx, par__conOy, par__conOz                % coordinates of conoid origin in clavicle frame
Constants par__conIx, par__conIy, par__conIz                % coordinates of conoid insertion in scapula frame
P_claviclerJoint_conO> = Vector(clavicler,par__conOx,par__conOy,par__conOz)					
P_scapularJoint_conI> = Vector(scapular,par__conIx,par__conIy,par__conIz)				
V2pts(ground, clavicler, claviclerJoint, conO)
V2pts(ground, scapular, scapularJoint, conI)		
LOI = Mag( P_conO_conI> ) 												% Distance between conoid origin and insertion
Stretch = LOI - par__conL												% Stretch of conoid "spring"
StretchPositive := 0.5*(Stretch + sqrt(Stretch^2 + par__epsconoid^2))	% Eliminate negative stretches
Uvec> = P_conO_conI> / LOI 												% Unit vector from conoid origin to insertion
Force(conO/conI, -par__conK*StretchPositive*Uvec> ) 				    % Spring force

%--------------------------------------------------------------------
% Determine the YZY rotation angles from thorax to humerus
% Apply the external torques on the axes of this joint coordinate system
%--------------------------------------------------------------------
% This is the symbolic form of an YZY rotation matrix:
%  frames a,b,c,d
%  variables q1,q2,q3
%  simprot(a,b,2,q1);
%  simprot(b,c,3,q2);
%  simprot(c,d,2,q3);
%  a_d
% Result[1,1] = COS(q1)*COS(q2)*COS(q3) - SIN(q1)*SIN(q3)
% Result[1,2] = -SIN(q2)*COS(q1)
% Result[1,3] = SIN(q1)*COS(q3) + SIN(q3)*COS(q1)*COS(q2)
% Result[2,1] = SIN(q2)*COS(q3)
% Result[2,2] = COS(q2)
% Result[2,3] = SIN(q2)*SIN(q3)
% Result[3,1] = -SIN(q3)*COS(q1) - SIN(q1)*COS(q2)*COS(q3)
% Result[3,2] = SIN(q1)*SIN(q2)
% Result[3,3] = COS(q1)*COS(q3) - SIN(q1)*SIN(q3)*COS(q2)

% so we can solve the YZY angles as follows (if we are not in gimbal lock where sin(q2) is zero)
qTH = [														  &
	sym_atan2(thorax_humerusr[3,2] , (-thorax_humerusr[1,2])) 	; &		% plane of elevation, between -PI and PI
	acos(thorax_humerusr[2,2])								; &		% elevation angle, between 0 and PI
	sym_atan2(thorax_humerusr[2,3] , thorax_humerusr[2,1])	  &		% humerus internal rotation, between -PI and PI
	]	
encode qTH

% calculate the floating axis (where we need to apply the Z moment)
floatingaxis> = cross(thorax2>,humerusr2>)
floatingaxis> = floatingaxis>/mag(floatingaxis>)		% normalize to unit vector

% apply an external moment on humerus that is equivalent to the JCS moments that we wanted
variables MTHy, MTHz, MTHyy								% the three JCS moments between thorax and humerus
variables MHx, MHy, MHz									% the XYZ component of the (equivalent) external moment
MH> = vector(thorax, MHx, MHy, MHz)

% equations that say that JCS moments are projections of the external moment:
equation[1] = dot(MH>, thorax2>) 		- MTHy
equation[2] = dot(MH>, floatingaxis>) 	- MTHz
equation[3] = dot(MH>, humerusr2>) 	- MTHyy

% now we can solve the external moment (will fail in gimbal lock!):
solve(equation, MHx, MHy, MHz)

% and we apply this external moment to the humerus
Torque(thorax/humerusr, MH>)		

%--------------------------------------------------
% Scapulothoracic contact forces on TS and AI
%--------------------------------------------------
Constants par__Mx,par__My,par__Mz                        % center of the thorax ellipsoid
Constants par__Ax,par__Ay,par__Az                        % radii of the thorax ellipsoid
Constants par__epscontact, par__kcontact             % parameters for scapula-thorax contact model

Points CenterEllipsoid
P_thoraxJoint_CenterEllipsoid> = Vector(thorax, par__Mx, par__My, par__Mz)

% compute the coordinates of the ellipsoid center relative to thorax
PxCE = dot(P_thoraxO_CenterEllipsoid>, thorax1>)
PyCE = dot(P_thoraxO_CenterEllipsoid>, thorax2>)
PzCE = dot(P_thoraxO_CenterEllipsoid>, thorax3>)

% Create the point TS and its location x,y,z on the scapula
Points TS													
Constants par__TSx, par__TSy, par__TSz
P_scapularJoint_TS> = Vector(scapular, par__TSx, par__TSy, par__TSz)
V2pts(ground, scapular, scapularJoint, TS)			% Autolev will need velocity of contact point in equations of motion

% compute its global coordinates (relative to thorax)
PxTS = dot(P_thoraxO_TS>, thorax1>)
PyTS = dot(P_thoraxO_TS>, thorax2>)
PzTS = dot(P_thoraxO_TS>, thorax3>)

% computes the "distance" to thorax surface
FTS = ((PxTS-PxCE)/par__Ax)^2 + ((PyTS-PyCE)/par__Ay)^2 + ((PzTS-PzCE)/par__Az)^2 - 1

% attenuate F when it is positive
FminusTS = 0.5*(FTS - sqrt(FTS^2 + par__epscontact^2))

% compute the forces						
FxTS = -par__kcontact*PxTS*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Ax^2*FminusTS
FyTS = -par__kcontact*PyTS*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Ay^2*FminusTS
FzTS = -par__kcontact*PzTS*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Az^2*FminusTS

% apply the force to the point TS
Force_TS> += FxTS*thorax1> + FyTS*thorax2> + FzTS*thorax3>

% Create the point AI and its location x,y,z on the scapula
Points AI													
Constants par__AIx, par__AIy, par__AIz
P_scapularJoint_AI> = Vector(scapular, par__AIx, par__AIy, par__AIz)
V2pts(ground, scapular, scapularJoint, AI)			% Autolev will need velocity of contact point in equations of motion

% compute its global coordinates (relative to thorax)
PxAI = dot(P_thoraxO_AI>, thorax1>)
PyAI = dot(P_thoraxO_AI>, thorax2>)
PzAI = dot(P_thoraxO_AI>, thorax3>)

% computes the "distance" to thorax surface
FAI = ((PxAI-PxCE)/par__Ax)^2 + ((PyAI-PyCE)/par__Ay)^2 + ((PzAI-PzCE)/par__Az)^2 - 1

% attenuate F when it is positive
FminusAI = 0.5*(FAI - sqrt(FAI^2 + par__epscontact^2))

% compute the forces						
FxAI = -par__kcontact*PxAI*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Ax^2*FminusAI
FyAI = -par__kcontact*PyAI*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Ay^2*FminusAI
FzAI = -par__kcontact*PzAI*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Az^2*FminusAI

% apply the force to the point AI
Force_AI> += FxAI*thorax1> + FyAI*thorax2> + FzAI*thorax3>

% Also encode those so we have them in the MEX function
F_SCAP = [FxTS , FyTS, FzTS; FxAI , FyAI, FzAI];
% Also encode the TS and AI coordinates and the "distance" to thorax surface
dist_SCAP = [FTS; FAI];
pos_SCAP = [PxTS , PyTS, PzTS; PxAI , PyAI, PzAI];
encode F_SCAP, dist_SCAP, pos_SCAP

%----------------------------------------------------------
% GH net reaction force acting on the Scapula
%----------------------------------------------------------

% compute net reaction force from weight and accelerations of arm segments
% this will be the actual reaction force if the accelerations satisfy the equations of motion (Zero = 0)
FGH> = 		par__humerusrMass * (-A_humerusrO_ground> + GravityAxis>)	&
		+	par__ulnarMass    * (-A_ulnarO_ground>    + GravityAxis>)	&
		+	par__radiusrMass  * (-A_radiusrO_ground>  + GravityAxis>)	&
		+	par__handrMass    * (-A_handrO_ground>    + GravityAxis>)	

% express the result as XYZ force components in Scapula reference frame
F_GH = [	dot(FGH>, scapular1>) ; dot(FGH>, scapular2>) ; dot(FGH>, scapular3>) ]

% we need to report F_GH back to the simulation program to help detect shoulder dislocation
% Note that this is the net intersegmental force, not accounting for muscles.  Muscle contributions
% are done in das3.c
encode F_GH		

%--------------------------------------------------------------------
% Generate equations of motion
%--------------------------------------------------------------------

Zero = Fr() + FrStar()

%---------------------------------------------------------------------------------------
% For visualization, we export position and orientation of the segments in a matrix.
% This output is only used for model testing in Matlab.  Required processor time is
% very small compared to dynamics, so we just generate this always.
%---------------------------------------------------------------------------------------

vis = [dot(P_groundO_thoraxJoint>, ground1>), dot(P_groundO_thoraxJoint>, ground2>), dot(P_groundO_thoraxJoint>, ground3>), rows(ground_thorax,1), rows(ground_thorax,2), rows(ground_thorax,3);&
		dot(P_groundO_clavicle1Joint>, ground1>), dot(P_groundO_clavicle1Joint>, ground2>), dot(P_groundO_clavicle1Joint>, ground3>), rows(ground_clavicle1,1), rows(ground_clavicle1,2), rows(ground_clavicle1,3);&
		dot(P_groundO_clavicle2Joint>, ground1>), dot(P_groundO_clavicle2Joint>, ground2>), dot(P_groundO_clavicle2Joint>, ground3>), rows(ground_clavicle2,1), rows(ground_clavicle2,2), rows(ground_clavicle2,3);&
		dot(P_groundO_claviclerJoint>, ground1>), dot(P_groundO_claviclerJoint>, ground2>), dot(P_groundO_claviclerJoint>, ground3>), rows(ground_clavicler,1), rows(ground_clavicler,2), rows(ground_clavicler,3);&
		dot(P_groundO_scapula1Joint>, ground1>), dot(P_groundO_scapula1Joint>, ground2>), dot(P_groundO_scapula1Joint>, ground3>), rows(ground_scapula1,1), rows(ground_scapula1,2), rows(ground_scapula1,3);&
		dot(P_groundO_scapula2Joint>, ground1>), dot(P_groundO_scapula2Joint>, ground2>), dot(P_groundO_scapula2Joint>, ground3>), rows(ground_scapula2,1), rows(ground_scapula2,2), rows(ground_scapula2,3);&
		dot(P_groundO_scapularJoint>, ground1>), dot(P_groundO_scapularJoint>, ground2>), dot(P_groundO_scapularJoint>, ground3>), rows(ground_scapular,1), rows(ground_scapular,2), rows(ground_scapular,3);&
		dot(P_groundO_humerus1Joint>, ground1>), dot(P_groundO_humerus1Joint>, ground2>), dot(P_groundO_humerus1Joint>, ground3>), rows(ground_humerus1,1), rows(ground_humerus1,2), rows(ground_humerus1,3);&
		dot(P_groundO_humerus2Joint>, ground1>), dot(P_groundO_humerus2Joint>, ground2>), dot(P_groundO_humerus2Joint>, ground3>), rows(ground_humerus2,1), rows(ground_humerus2,2), rows(ground_humerus2,3);&
		dot(P_groundO_humerusrJoint>, ground1>), dot(P_groundO_humerusrJoint>, ground2>), dot(P_groundO_humerusrJoint>, ground3>), rows(ground_humerusr,1), rows(ground_humerusr,2), rows(ground_humerusr,3);&
		dot(P_groundO_ulnarJoint>, ground1>), dot(P_groundO_ulnarJoint>, ground2>), dot(P_groundO_ulnarJoint>, ground3>), rows(ground_ulnar,1), rows(ground_ulnar,2), rows(ground_ulnar,3);&
		dot(P_groundO_radiusrJoint>, ground1>), dot(P_groundO_radiusrJoint>, ground2>), dot(P_groundO_radiusrJoint>, ground3>), rows(ground_radiusr,1), rows(ground_radiusr,2), rows(ground_radiusr,3);&
		dot(P_groundO_handrJoint>, ground1>), dot(P_groundO_handrJoint>, ground2>), dot(P_groundO_handrJoint>, ground3>), rows(ground_handr,1), rows(ground_handr,2), rows(ground_handr,3)];

encode vis

%-----------------------------------------------------------------------------------------------------
% For implicit dynamics: generate symbolic expressions for Zero and Jacobians d(Zero)/d(q,qd,qdd)
%-----------------------------------------------------------------------------------------------------

% q are the generalized coordinates
q = [SCy, &
SCz, &
SCx, &
ACy, &
ACz, &
ACx, &
GHy, &
GHz, &
GHyy, &
ELx, &
PSy]

% generalized velocities
qd = [SCy', &
SCz', &
SCx', &
ACy', &
ACz', &
ACx', &
GHy', &
GHz', &
GHyy', &
ELx', &
PSy']

% generalized accelerations
qdd = [SCy'', &
SCz'', &
SCx'', &
ACy'', &
ACz'', &
ACx'', &
GHy'', &
GHz'', &
GHyy'', &
ELx'', &
PSy'']

% generate expressions for implicit equation of motion and its Jacobians
dz_dq = ZEE(D(Zero,q))
dz_dqd = ZEE(D(Zero,qd))
dz_dqdd = ZEE(D(Zero,qdd))
encode Zero, dz_dq, dz_dqd, dz_dqdd

% write all Encoded expressions to C file
Code Algebraic() das3_al_raw.c

EXIT
