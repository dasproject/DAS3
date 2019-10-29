function das3_makeMEX(model,MEXtemplate,newMEX)
	% Creates MEX file and header file from Opensim model osimfile
	% MEXtemplate is the file used to fill in all the MEX content that is
	% independent of the model
	% the new MEX file is saved in newMEX.c
	%
	% Dimitra Blana, September 2014

	[pathstr,fname,~] = fileparts(newMEX);
	mexhfile = fullfile(pathstr,[fname,'.h']);

	% Open MEX template for reading and new MEX file for writing
	infile = fopen(MEXtemplate, 'rt');
	if infile==-1
		errordlg('Could not open MEX template file for reading','File Error');
		return;
	end
	outfile = fopen(newMEX, 'wt');
	if outfile==-1
		errordlg('Could not open new MEX file for writing','File Error');
		return;
	end
	mexh_file = fopen(mexhfile, 'wt');
	if mexh_file==-1
		errordlg('Could not open MEX header file for writing','File Error');
		return;
	end

	% Start filling in the header file
	fprintf(mexh_file,'// %s.h\n\n',fname);

	fprintf(mexh_file,'// This file must be included in both %s.c and Autolev %s_al.c to ensure that declarations\n',fname, fname);
	fprintf(mexh_file,'// and constants remain consistent between them.\n\n');

	fprintf(mexh_file,'#define NDOF %d			// number of kinematic degrees of freedom\n',model.nDofs);
	fprintf(mexh_file,'#define NVIS %d			// number of rigid bodies for visualization\n\n',model.nSegments-1);

	fprintf(mexh_file,'// struct with model parameters that must be shared with Autolev code\n');
	fprintf(mexh_file,'typedef struct {\n');
	fprintf(mexh_file,'	// Gravity parameters\n');
	fprintf(mexh_file,'	double GravAxisX, GravAxisY, GravAxisZ;\n\n');

	fprintf(mexh_file,'	// scapula-thorax contact force model\n');
	fprintf(mexh_file,'	double epscontact, kcontact;            // epsilon and stiffness in model\n');
	fprintf(mexh_file,'	double Ax, Ay, Az, Mx, My, Mz;          // radii and center of the thorax ellipsoid\n');
	fprintf(mexh_file,'	double TSx, TSy, TSz, AIx, AIy, AIz;    // coordinates of TS and AI projections in scapula frame\n\n');

	fprintf(mexh_file,'	// conoid force model\n');
	fprintf(mexh_file,'	double epsconoid, conK, conL;       // epsilon, stiffness and length of conoid ligament\n');
	fprintf(mexh_file,'	double  conOx, conOy, conOz;        // coordinates of conoid origin in clavicle frame\n');
	fprintf(mexh_file,'	double  conIx, conIy, conIz;        // coordinates of conoid insertion in scapula frame\n');

	% Start filling in the c MEX file
	fprintf(outfile,'/*======================================================================================\n');
	fprintf(outfile,' *\n');
	fprintf(outfile,' * %s.c\n',fname);
	fprintf(outfile,' *\n');
	fprintf(outfile,' * Implicit differential equation for 3D musculoskeletal model : f(x,dx/dt,u) = 0\n');
	fprintf(outfile,' *\n');
	fprintf(outfile,' * This is the source code for the MEX function %s.mexw32\n',fname);
	fprintf(outfile,' * The musculoskeletal model is documented in the file %s.pdf\n',fname);
	fprintf(outfile,' * The function API documentation is in %s.m\n',fname);
	fprintf(outfile,' *\n');
	fprintf(outfile,' *=====================================================================================*/\n\n');
	fprintf(outfile,'#include <stdio.h>\n');
	fprintf(outfile,'#include <stdlib.h>\n');
	fprintf(outfile,'#include <math.h>\n');
	fprintf(outfile,'#include <string.h>\n');
	fprintf(outfile,'#include "mex.h"\n');
	fprintf(outfile,'#include "%s.h"\n\n',fname);

	% Copy everything from the template MEX file up to the "#define NSEG"
	tline = fgetl(infile);
	while (feof(infile)==0)&&(~strncmp(tline, '#define NSEG',12))
		fprintf(outfile,'%s\n',tline);
		tline = fgetl(infile);
	end
	fprintf(outfile,'#define NSEG %d  // number of segments\n',model.nSegments);
	tline = fgetl(infile);

	while (feof(infile)==0)&&(~strcmp(tline, '        // from P, extract the joint positions and axes, and segment mass properties'))
		fprintf(outfile,'%s\n',tline);
		tline = fgetl(infile);
	end
	fprintf(outfile,'%s\n',tline);

	% six transforms in "custom joints"
	% in this version, we ignore translations and only look at rotations!
	%transcoords = {'r1_coordinates', 'r2_coordinates', 'r3_coordinates', 't1_coordinates', 't2_coordinates', 't3_coordinates'};
	transcoords = {'r1_coordinates', 'r2_coordinates', 'r3_coordinates'};
	fprintf(outfile,'        // Joint Properties: \n');
	fprintf(mexh_file,'\n	// Joint parameters\n');
	fprintf(outfile,'        {\n');
	for ijnt=1:model.nJoints
		current_joint = model.joints{ijnt};
		fprintf(outfile,'		parameters.%sX = extract2(P,"joints",%d,"location",1);\n',current_joint.name,ijnt);		
		fprintf(outfile,'		parameters.%sY = extract2(P,"joints",%d,"location",2);\n',current_joint.name,ijnt);		
		fprintf(outfile,'		parameters.%sZ = extract2(P,"joints",%d,"location",3);\n',current_joint.name,ijnt);		
		fprintf(mexh_file,'	double %sX, %sY, %sZ;                   // coordinates of joint in parent frame\n',current_joint.name,current_joint.name,current_joint.name);
		
		% find transforms for current joint
		fields = isfield(current_joint, transcoords);
		transforms = find(fields);

		% if there are no transforms, this is a weld joint
		if isempty(transforms)
			fprintf(outfile,'		parameters.%sXrot = extract2(P,"joints",%d,"orientation",1);\n',current_joint.name,ijnt);
			fprintf(outfile,'		parameters.%sYrot = extract2(P,"joints",%d,"orientation",2);\n',current_joint.name,ijnt);
			fprintf(outfile,'		parameters.%sZrot = extract2(P,"joints",%d,"orientation",3);\n\n',current_joint.name,ijnt);
			fprintf(mexh_file,'	double %sXrot, %sYrot, %sZrot;      // orientation of weld joint in parent frame\n\n',current_joint.name,current_joint.name,current_joint.name);
		else
			for i_intframes = 1:length(transforms)
				current_transform_coords = transcoords{transforms(i_intframes)};
				current_transform = strrep(current_transform_coords, '_coordinates', '');
				current_transform_coefs = [current_transform '_coefs'];
				coefs = getfield(current_joint, current_transform_coefs);  
				dofs = getfield(current_joint, current_transform_coords);       
				dof = dofs{1};  % only one coordinate supported

				% check if this is a constrained coordinate
				for iconst=1:model.nConstraints
					if strcmp(model.constraints{iconst}.dependent_dof,dof)
						inddof = model.constraints{iconst}.independent_dof{1};
						indcoefs = model.constraints{iconst}.coefs;
						if length(indcoefs)==2
							if indcoefs(1)~=0 || indcoefs(2)~=1
								fprintf(outfile, '		parameters.%sd0 = extract2(P,"constraints",%d,"coefs",1);\n',dof,iconst);
								fprintf(outfile, '		parameters.%sd1 = extract2(P,"constraints",%d,"coefs",2);\n',dof,iconst);
								fprintf(mexh_file,'	double %sd0, %sd1;              // intercept and slope of linear coupled coordinate function\n',dof,dof);                        
							end
						else
							fprintf(outfile, '		parameters.%sd1 = extract2(P,"constraints",%d,"coefs",1);\n',dof,iconst);
							fprintf(outfile, '		parameters.%sd2 = extract2(P,"constraints",%d,"coefs",2);\n',dof,iconst);
							fprintf(outfile, '		parameters.%sd3 = extract2(P,"constraints",%d,"coefs",3);\n',dof,iconst);
							fprintf(outfile, '		parameters.%sd4 = extract2(P,"constraints",%d,"coefs",4);\n',dof,iconst);
							fprintf(outfile, '		parameters.%sd5 = extract2(P,"constraints",%d,"coefs",5);\n',dof,iconst);
							fprintf(outfile, '		parameters.%sd6 = extract2(P,"constraints",%d,"coefs",6);\n',dof,iconst);
							fprintf(mexh_file,'	double %sd1, %sd2, %sd3, %sd4, %sd5, %sd6;  // parameters of coupled coordinate spline function\n',dof,dof,dof,dof,dof,dof); 
						end
						break;                
					end
				end
				
				if length(coefs)==2
					if coefs(1)~=0 || coefs(2)~=1
						fprintf(outfile, '		parameters.%s0 = extract2(P,"joints",%d,"%s_coefs",1);\n',dof,ijnt,current_transform);
						fprintf(outfile, '		parameters.%s1 = extract2(P,"joints",%d,"%s_coefs",2);\n',dof,ijnt,current_transform);
						fprintf(mexh_file,'	double %s0, %s1;                // intercept and slope of linear transform function\n',dof,dof); 
					end
				else
					fprintf(outfile, '		parameters.%s1 = extract2(P,"joints",%d,"%s_coefs",1);\n',dof,ijnt,current_transform);
					fprintf(outfile, '		parameters.%s2 = extract2(P,"joints",%d,"%s_coefs",2);\n',dof,ijnt,current_transform);
					fprintf(outfile, '		parameters.%s3 = extract2(P,"joints",%d,"%s_coefs",3);\n',dof,ijnt,current_transform);
					fprintf(outfile, '		parameters.%s4 = extract2(P,"joints",%d,"%s_coefs",4);\n',dof,ijnt,current_transform);
					fprintf(outfile, '		parameters.%s5 = extract2(P,"joints",%d,"%s_coefs",5);\n',dof,ijnt,current_transform);
					fprintf(outfile, '		parameters.%s6 = extract2(P,"joints",%d,"%s_coefs",6);\n',dof,ijnt,current_transform);
					fprintf(mexh_file,'	double %s1, %s2, %s3, %s4, %s5, %s6;    // parameters of spline transform function\n',dof,dof,dof,dof,dof,dof); 
				end
				fprintf(outfile, '		parameters.%sAxisX = extract2(P,"joints",%d,"%s_axis",1);\n',dof,ijnt,current_transform);
				fprintf(outfile, '		parameters.%sAxisY = extract2(P,"joints",%d,"%s_axis",2);\n',dof,ijnt,current_transform);
				fprintf(outfile, '		parameters.%sAxisZ = extract2(P,"joints",%d,"%s_axis",3);\n',dof,ijnt,current_transform);
				fprintf(outfile,'		normalize(&parameters.%sAxisX,    &parameters.%sAxisY,    &parameters.%sAxisZ);\n\n',dof,dof,dof);	
				fprintf(mexh_file,'	double %sAxisX, %sAxisY, %sAxisZ;       // rotation axis orientation in parent frame\n\n',dof,dof,dof);           
			end         
		end
	end
	fprintf(outfile,'        }\n\n');

	fprintf(outfile,'        // Segment Properties: \n');
	fprintf(mexh_file,'\n	// Segment parameters\n');
	fprintf(outfile,'        {\n');

	for iseg = 1:model.nSegments
		current_seg = model.segments{iseg};
		if current_seg.mass<=0.0001, continue; end
		fprintf(outfile,'        EXTRACTINERTIAL(%s, %d);\n',current_seg.name,iseg);
		fprintf(mexh_file,'	double %sMass; // mass\n',current_seg.name);
		fprintf(mexh_file,'	double %sCMx, %sCMy, %sCMz; // center of mass\n',current_seg.name,current_seg.name,current_seg.name);
		fprintf(mexh_file,'	double %sIxx, %sIyy, %sIzz, %sIxy, %sIyz, %sIxz; // moments of inertia\n\n',current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name);
	end
	fprintf(outfile,'        }\n');

	tline = fgetl(infile);
	while (feof(infile)==0)
		fprintf(outfile,'%s\n',tline);
		tline = fgetl(infile);
	end
	fprintf(outfile,'%s\n',tline);

	fprintf(mexh_file,'} param_struct;\n\n');

    % watch out: the prototype for the Autolev C function must be identical
    % to the function definition in autolevclean (in das3_buildmodel.m)
    % if not, hopefully the C compiler will give an error message
	fprintf(mexh_file,'// prototype for the Autolev C function\n');
	fprintf(mexh_file,'void %s_al(\n',fname);
	fprintf(mexh_file,'	param_struct* param,				// input: pointer to struct containing parameter\n');
	fprintf(mexh_file,'	double q[NDOF], 					// input: joint angles\n');
	fprintf(mexh_file,'	double qd[NDOF], 					// input: angular velocities\n');
	fprintf(mexh_file,'	double qdd[NDOF],					// input: angular accelerations\n');
	fprintf(mexh_file,'	double mTH[3],                      // input: external moments applied to thoracohumeral axes\n'); 
	fprintf(mexh_file,'	double exF[2],                      // input: forearm support force\n');
	fprintf(mexh_file,'	double handF[3],                    // input: external force applied to the CoM of the hand\n');
	fprintf(mexh_file,'	double Zero[NDOF],					// output: dynamics imbalance\n');
	fprintf(mexh_file,'	double dz_dq[NDOF][NDOF], 			// output: Jacobian of Zero with respect to q\n');
	fprintf(mexh_file,'	double dz_dqd[NDOF][NDOF], 			// output: Jacobian of Zero with respect to qdot\n');
	fprintf(mexh_file,'	double dz_dqdd[NDOF][NDOF], 		// output: Jacobian of Zero with respect to qdotdot\n');
	fprintf(mexh_file,'	double F_GH[3],                     // output: 3D GH contact force, acting on scapula, expressed in scapula reference frame\n');
	fprintf(mexh_file,'	double F_SCAP[2][3],                // output: 3D Contact forces acting on TS and AI, expressed in thorax reference frame\n');
	fprintf(mexh_file,'	double dist_SCAP[2],                // output: "distance" of TS and AI to thorax surface\n');
	fprintf(mexh_file,'	double pos_SCAP[2][3],              // output: TS and AI coordinates in thorax frame\n');
	fprintf(mexh_file,'	double vis[NVIS][12],				// output: position and orientation of bones\n');
	fprintf(mexh_file,'	double qTH[3]);                     // output: angles between thorax and humerus (YZY sequence)\n');

	fclose(mexh_file);
	fclose(infile);
	fclose(outfile);
end