function das3_makeAutolev(model, Autolevfile)
	% Creates Autolev file from Opensim model
	%
	% First version by Dimitra Blana, October 2014

	% Open Autolev file for writing
	autolev_file = fopen(Autolevfile, 'wt');
	if autolev_file==-1
		error('Could not create %s. It may be open in another application', Autolevfile);
		return;
    end
    fprintf('Writing %s...\n', Autolevfile);

	fprintf(autolev_file, '%% %s\n',Autolevfile);
	fprintf(autolev_file, '%% This is an input file for Autolev to generate dynamic equations of motion for the %s model.\n\n',model.name);
	fprintf(autolev_file, 'AUTOZ ON\n');
	fprintf(autolev_file, 'OVERWRITE ALL\n');
	fprintf(autolev_file, 'PAUSE 0\n');
	fprintf(autolev_file, 'UnitSystem  kg,meter,sec\n\n');

	% Add the motion variables
	for idof = 1:model.nDofs
		fprintf(autolev_file, 'MotionVariables'' %s''''\n',model.dofs{idof}.name);
	end
	fprintf(autolev_file, '\n');

	% Add segments

	% six transforms in "custom joints"
	% in this version, we ignore translations and only look at rotations!
	%transcoords = {'r1_coordinates', 'r2_coordinates', 'r3_coordinates', 't1_coordinates', 't2_coordinates', 't3_coordinates'};
	transcoords = {'r1_coordinates', 'r2_coordinates', 'r3_coordinates'};

	% First, create the Newtonian frame, called 'ground'
    % (Opensim 4.0 has a ground, but it is not a body, so it was not picked up
    % by readosim.m)
    newtonian = 'ground';
	fprintf(autolev_file, '%%------------------------------------\n');
    fprintf(autolev_file, '%% Newtonian frame\n');
	fprintf(autolev_file, '%%------------------------------------\n');
	fprintf(autolev_file, 'Newtonian %s\n\n',newtonian);

	% Add all other segments
    allsegs = 1:model.nSegments;
	seen_dofs = cell(1,1);
	dof_index=1;

	for iseg = allsegs
		current_seg = model.segments{iseg};
		% find the joint linking this body/frame with the parent
		for ijnt=1:model.nJoints
			if strcmp(model.joints{ijnt}.name,current_seg.parent_joint), break; end
		end
		current_joint = model.joints{ijnt};
		prev_seg_name = current_joint.parent_segment;
		
		fprintf(autolev_file, '%%------------------------------------\n');
		fprintf(autolev_file, '%% Segment %s\n',current_seg.name);
		fprintf(autolev_file, '%%------------------------------------\n');
		
		% find transforms for current joint
		fields = isfield(current_joint, transcoords);
		transforms = find(fields);
		
		% if there are no transforms, this is a weld joint
		if isempty(transforms)
			if current_seg.mass<=0.0001
				% this is ~massless, so make it a frame
				fprintf(autolev_file, 'Frame %s\n\n',current_seg.name);
			else
				fprintf(autolev_file, 'Body %s\n',current_seg.name);
				fprintf(autolev_file, 'Constants par__%sMass\n',current_seg.name);
				fprintf(autolev_file, 'Mass %s = par__%sMass\n',current_seg.name,current_seg.name);
				fprintf(autolev_file, 'Constants par__%sIxx, par__%sIyy, par__%sIzz, par__%sIxy, par__%sIyz, par__%sIxz\n',current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name);
				fprintf(autolev_file, 'Inertia %s, par__%sIxx, par__%sIyy, par__%sIzz, par__%sIxy, par__%sIyz, par__%sIxz\n\n',current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name);
			end
			
			% orientation is specified as an X-Y-Z rotation
			fprintf(autolev_file, 'Frames %s1, %s2\n',current_seg.name,current_seg.name);
			fprintf(autolev_file, 'Constants par__%sXrot, par__%sYrot, par__%sZrot\n',current_joint.name,current_joint.name,current_joint.name);
			fprintf(autolev_file, 'Simprot(%s, %s1, 1, par__%sXrot)\n',prev_seg_name,current_seg.name,current_joint.name);
			fprintf(autolev_file, 'Simprot(%s1, %s2, 2, par__%sYrot)\n',current_seg.name,current_seg.name,current_joint.name);
			fprintf(autolev_file, 'Simprot(%s2, %s,  3, par__%sZrot);	%% orientation of %s\n\n',current_seg.name,current_seg.name,current_joint.name,current_seg.name);
			
			fprintf(autolev_file, 'Point %sJoint\n',current_seg.name);
			
			if strcmp(prev_seg_name, newtonian)
				% This body is attached to the newtonian frame
				fprintf(autolev_file, 'Constants par__%sX, par__%sY, par__%sZ\n',current_joint.name,current_joint.name,current_joint.name);
				fprintf(autolev_file, 'P_%sO_%sJoint> = Vector(%s, par__%sX, par__%sY, par__%sZ)	%% location in parent\n',newtonian,current_seg.name,newtonian,current_joint.name,current_joint.name,current_joint.name);
				fprintf(autolev_file, 'V_%sJoint_%s> = 0> %% velocity of joint\n', current_seg.name, newtonian);
				fprintf(autolev_file, 'A_%sJoint_%s> = 0> %% acceleration of joint\n\n', current_seg.name, newtonian);
				fprintf(autolev_file, 'Constants par__%sCMx, par__%sCMy, par__%sCMz\n',current_seg.name,current_seg.name,current_seg.name);
				fprintf(autolev_file, 'P_%sJoint_%sO> = Vector(%s, par__%sCMx, par__%sCMy, par__%sCMz)  %% location of CM\n',current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name);
				fprintf(autolev_file, 'V_%sO_%s> = 0> %% velocity of CM\n', current_seg.name, newtonian);
				fprintf(autolev_file, 'A_%sO_%s> = 0> %% acceleration of CM\n', current_seg.name, newtonian);
				fprintf(autolev_file, 'W_%s_%s> = 0>\n', current_seg.name, newtonian);
				fprintf(autolev_file, 'ALF_%s_%s> = 0>\n\n', current_seg.name, newtonian);
			else
				fprintf(autolev_file, 'Constants par__%sX, par__%sY, par__%sZ\n',current_joint.name,current_joint.name,current_joint.name);
				fprintf(autolev_file, 'P_%sJoint_%sJoint> = Vector(%s, par__%sX, par__%sY, par__%sZ)	%% location in parent\n',prev_seg_name,current_seg.name,prev_seg_name,current_joint.name,current_joint.name,current_joint.name);
				fprintf(autolev_file, 'V2pts(%s, %s, %sJoint, %sJoint) %% velocity of joint\n', newtonian, prev_seg_name, prev_seg_name, current_seg.name);
				fprintf(autolev_file, 'A2pts(%s, %s, %sJoint, %sJoint) %% acceleration of joint\n\n', newtonian, prev_seg_name, prev_seg_name, current_seg.name);
				fprintf(autolev_file, 'Constants par__%sCMx, par__%sCMy, par__%sCMz\n',current_seg.name,current_seg.name,current_seg.name);
				fprintf(autolev_file, 'P_%sJoint_%sO> = Vector(%s, par__%sCMx, par__%sCMy, par__%sCMz)  %% location of CM\n',current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name,current_seg.name);
				fprintf(autolev_file, 'V2pts(%s, %s, %sJoint, %sO) %% velocity of CM\n', newtonian, current_seg.name, current_seg.name, current_seg.name);
				fprintf(autolev_file, 'A2pts(%s, %s, %sJoint, %sO) %% acceleration of CM\n\n', newtonian, current_seg.name, current_seg.name, current_seg.name);
			end
		else
			
			% if there are more than one transforms, we need to add
			% intermediate frames
			
			for i_intframes = 1:length(transforms)
				if i_intframes<length(transforms)
					% intermediate frame
					newframename = [current_seg.name num2str(i_intframes)];
					newframename_mass = 0;
					fprintf(autolev_file, 'Frame %s\n\n',newframename);
				else
					newframename = current_seg.name;
					newframename_mass = current_seg.mass;
					if  newframename_mass<=0.0001
						% this is ~massless, so make it a frame
						fprintf(autolev_file, 'Frame %s\n\n',newframename);
					else
						fprintf(autolev_file, 'Body %s\n',newframename);
						fprintf(autolev_file, 'Constants par__%sMass\n',newframename);
						fprintf(autolev_file, 'Mass %s = par__%sMass\n',newframename,newframename);
						fprintf(autolev_file, 'Constants par__%sIxx, par__%sIyy, par__%sIzz, par__%sIxy, par__%sIyz, par__%sIxz\n',newframename,newframename,newframename,newframename,newframename,newframename);
						fprintf(autolev_file, 'Inertia %s, par__%sIxx, par__%sIyy, par__%sIzz, par__%sIxy, par__%sIyz, par__%sIxz\n\n',newframename,newframename,newframename,newframename,newframename,newframename,newframename);
					end
				end
				
				% find the current transform
				current_transform_coords = transcoords{transforms(i_intframes)};
				current_transform = strrep(current_transform_coords, '_coordinates', '');
				current_transform_coefs = [current_transform '_coefs'];
				coefs = getfield(current_joint, current_transform_coefs);
				dofs = getfield(current_joint, current_transform_coords);
				dof = dofs{1};  % only one coordinate supported
				seen_dofs{dof_index} = dof;
				dof_index = dof_index+1;
				
				% check if this is a constrained coordinate
				constrainedflag = 0;
				for iconst=1:model.nConstraints
					if strcmp(model.constraints{iconst}.dependent_dof,dof)
						inddof = model.constraints{iconst}.independent_dof{1};
						indcoefs = model.constraints{iconst}.coefs;
						constrainedflag=1;
						if length(indcoefs)==2
							% this is a linear function
							if indcoefs(1)==0 && indcoefs(2)==1
								% the constrained coordinate equals the
								% independent coordinate
								fprintf(autolev_file, '%s = %s\n',dof,inddof);
								fprintf(autolev_file, '%s_dot = %s''\n',dof,inddof);
								fprintf(autolev_file, '%s_dotdot = %s''''\n',dof,inddof);
							else
								% it is a line
								fprintf(autolev_file, 'Constants par__%sd0,par__%sd1\n',dof,dof);
								fprintf(autolev_file, '%s = par__%sd0 + par__%sd1*%s\n',dof,dof,dof,inddof);
								fprintf(autolev_file, '%s_dot = par__%sd1*%s''\n',dof,dof,inddof);
								fprintf(autolev_file, '%s_dotdot = par__%sd1*%s''''\n\n',dof,dof,inddof);
							end
						else
							% it is a spline
							fprintf(autolev_file, 'Constants par__%sd1,par__%sd2,par__%sd3,par__%sd4,par__%sd5,par__%sd6\n',dof,dof,dof,dof,dof,dof);
							fprintf(autolev_file, '%s = spline(3,%s,par__%sd1,par__%sd2,par__%sd3,par__%sd4,par__%sd5,par__%sd6)\n',dof,inddof,dof,dof,dof,dof,dof,dof);
							fprintf(autolev_file, '%s_dot = Dt(%s)\n',dof,dof);
							fprintf(autolev_file, '%s_dotdot = Dt(%s_dot)\n\n',dof,dof);
						end
						break;
					end
				end
				% check if this dof has been encountered earlier in the model
				seen_dof_flag = 0;
				for iseendof = 1:length(seen_dofs)-1
					if strcmp(seen_dofs{iseendof},dof)
						seen_dof_flag = 1;
						break;
					end
				end
				% if it has not been encountered, write its derivatives
				if ~constrainedflag && ~seen_dof_flag
					fprintf(autolev_file, '%s_dot = %s''\n',dof,dof);
					fprintf(autolev_file, '%s_dotdot = %s''''\n\n',dof,dof);
				end
				% now we have dof_dot and dof_dotdot whether dof is independent
				% or not
				
				fprintf(autolev_file, 'Constants par__%sAxisX, par__%sAxisY, par__%sAxisZ\n',dof,dof,dof);
				fprintf(autolev_file, '%sAxis> = Vector(%s, par__%sAxisX, par__%sAxisY, par__%sAxisZ)\n',dof,prev_seg_name,dof,dof,dof);
				% now check what kind of function the transform is
				if length(coefs)==2
					% it is an equality
					if coefs(1)==0 && coefs(2)==1
						fprintf(autolev_file, 'Simprot(%s,%s,%sAxis>,%s)\n\n',prev_seg_name,newframename,dof,dof);
						fprintf(autolev_file, 'W_%s_%s> = W_%s_%s> + %s_dot * %sAxis>     %% generates angular velocity of child body\n',newframename, newtonian, prev_seg_name, newtonian, dof, dof);
						fprintf(autolev_file, 'ALF_%s_%s> = ALF_%s_%s> + %s_dotdot * %sAxis>     %% generates angular acceleration of child body\n\n',newframename, newtonian, prev_seg_name, newtonian, dof, dof);
					else
						% it is a line
						fprintf(autolev_file, 'Constants par__%s0,par__%s1\n',dof,dof);
						dof_temp = sprintf('par__%s0 + par__%s1*%s',dof,dof,dof);
						fprintf(autolev_file, 'Simprot(%s,%s,%sAxis>,%s)\n\n',prev_seg_name,newframename,dof,dof_temp);
						fprintf(autolev_file, 'W_%s_%s> = W_%s_%s> + par__%s1*%s_dot * %sAxis>     %% generates angular velocity of child body\n',newframename, newtonian, prev_seg_name, newtonian, dof, dof, dof);
						fprintf(autolev_file, 'ALF_%s_%s> = ALF_%s_%s> + par__%s1*%s_dotdot * %sAxis>     %% generates angular acceleration of child body\n\n',newframename, newtonian, prev_seg_name, newtonian, dof, dof, dof);
					end
				else
					% it is a spline
					fprintf(autolev_file, 'Constants par__%s1,par__%s2,par__%s3,par__%s4,par__%s5,par__%s6\n',dof,dof,dof,dof,dof,dof);
					dof_temp = sprintf('spline(3,%s,par__%s1,par__%s2,par__%s3,par__%s4,par__%s5,par__%s6)',dof,dof,dof,dof,dof,dof,dof);
					fprintf(autolev_file, 'Simprot(%s,%s,%sAxis>,%s)\n\n',prev_seg_name,newframename,dof,dof_temp);
					fprintf(autolev_file, '%s_dot = Dt(%s)\n',dof_temp,dof_temp);
					fprintf(autolev_file, '%s_dotdot = Dt(%s_dot)\n\n',dof_temp,dof_temp);
					fprintf(autolev_file, 'W_%s_%s> = W_%s_%s> + %s_dot * %sAxis>     %% generates angular velocity of child body\n',newframename, newtonian, prev_seg_name, newtonian, dof_temp, dof);
					fprintf(autolev_file, 'ALF_%s_%s> = ALF_%s_%s> + %s_dotdot * %sAxis>     %% generates angular acceleration of child body\n\n',newframename, newtonian, prev_seg_name, newtonian, dof_temp, dof);
				end
				
				if i_intframes==length(transforms)
					% only for one body per segment (not the intermediate ones)                 
					fprintf(autolev_file, 'Point %sJoint\n',newframename);
					if strcmp(current_joint.parent_segment, newtonian)
						% This body is attached to the newtonian frame
						fprintf(autolev_file, 'Constants par__%sX, par__%sY, par__%sZ\n',current_joint.name,current_joint.name,current_joint.name);
						fprintf(autolev_file, 'P_%sO_%sJoint> = Vector(%s, par__%sX, par__%sY, par__%sZ)	%% location in parent\n',current_joint.parent_segment,newframename,current_joint.parent_segment,current_joint.name,current_joint.name,current_joint.name);
						fprintf(autolev_file, 'V_%sJoint_%s> = 0> %% velocity of joint\n', newframename, newtonian);
						fprintf(autolev_file, 'A_%sJoint_%s> = 0> %% acceleration of joint\n\n', newframename, newtonian);
					else
						fprintf(autolev_file, 'Constants par__%sX, par__%sY, par__%sZ\n',current_joint.name,current_joint.name,current_joint.name);
						fprintf(autolev_file, 'P_%sJoint_%sJoint> 	= Vector(%s, par__%sX, par__%sY, par__%sZ)	%% location in parent\n',current_joint.parent_segment,newframename,current_joint.parent_segment,current_joint.name,current_joint.name,current_joint.name);
						fprintf(autolev_file, 'V2pts(%s, %s, %sJoint, %sJoint) %% velocity of joint\n', newtonian, current_joint.parent_segment, current_joint.parent_segment, newframename);
						fprintf(autolev_file, 'A2pts(%s, %s, %sJoint, %sJoint) %% acceleration of joint\n\n', newtonian, current_joint.parent_segment, current_joint.parent_segment, newframename);
					end
					if newframename_mass>0.0001
						fprintf(autolev_file, 'Constants par__%sCMx, par__%sCMy, par__%sCMz\n',newframename,newframename,newframename);
						fprintf(autolev_file, 'P_%sJoint_%sO> = Vector(%s, par__%sCMx, par__%sCMy, par__%sCMz)  %% location of CM\n',newframename,newframename,newframename,newframename,newframename,newframename);
						fprintf(autolev_file, 'V2pts(%s, %s, %sJoint, %sO) %% velocity of CM\n', newtonian, newframename, newframename, newframename);
						fprintf(autolev_file, 'A2pts(%s, %s, %sJoint, %sO) %% acceleration of CM\n\n', newtonian, newframename, newframename, newframename);
					end
				end
				prev_seg_name = newframename;
			end
		end
	end

	% Gravity
	fprintf(autolev_file, '%%------------------------------------\n');
	fprintf(autolev_file, '%% Apply Gravity\n');
	fprintf(autolev_file, '%%------------------------------------\n');
	fprintf(autolev_file, 'Constants par__GravAxisX, par__GravAxisY, par__GravAxisZ\n');
	fprintf(autolev_file, 'GravityAxis> = Vector(%s, par__GravAxisX, par__GravAxisY, par__GravAxisZ)\n',newtonian);
	fprintf(autolev_file, 'Gravity( GravityAxis> )\n\n');

	% Forearm support
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% Apply force on the forearm (ulna) to simulate a mobile arm support\n');
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, 'Points ForF\n');
	fprintf(autolev_file, 'Variables  distF, ampF\n');
	fprintf(autolev_file, 'P_ulnarJoint_ForF> = Vector(ulnar,0,-distF,0)	\n');							
	fprintf(autolev_file, 'V2pts(%s, ulnar, ulnarJoint, ForF)\n',newtonian);
	fprintf(autolev_file, 'Force_ForF> = ampF*%s2> \n\n',newtonian);

	% Force at the hand
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% Apply force to the CoM of the hand\n');
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, 'Variables  handFx, handFy, handFz\n');
	fprintf(autolev_file, 'Force_handrO> = handFx*%s1> + handFy*%s2> + handFz*%s3> \n\n',newtonian,newtonian,newtonian);

	% Conoid force
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% Calculate and apply conoid force\n');
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, 'Points conO, conI							%% Conoid origin (on the clavicle) and insertion (on the scapula)\n');
	fprintf(autolev_file, 'Constants par__epsconoid, par__conK, par__conL				%% parameters for the conoid ligament model\n');
	fprintf(autolev_file, 'Constants par__conOx, par__conOy, par__conOz                %% coordinates of conoid origin in clavicle frame\n');
	fprintf(autolev_file, 'Constants par__conIx, par__conIy, par__conIz                %% coordinates of conoid insertion in scapula frame\n');
	fprintf(autolev_file, 'P_claviclerJoint_conO> = Vector(clavicler,par__conOx,par__conOy,par__conOz)					\n');
	fprintf(autolev_file, 'P_scapularJoint_conI> = Vector(scapular,par__conIx,par__conIy,par__conIz)				\n');
	fprintf(autolev_file, 'V2pts(%s, clavicler, claviclerJoint, conO)\n',newtonian);
	fprintf(autolev_file, 'V2pts(%s, scapular, scapularJoint, conI)		\n',newtonian);
	fprintf(autolev_file, 'LOI = Mag( P_conO_conI> ) 												%% Distance between conoid origin and insertion\n');
	fprintf(autolev_file, 'Stretch = LOI - par__conL												%% Stretch of conoid "spring"\n');
	fprintf(autolev_file, 'StretchPositive := 0.5*(Stretch + sqrt(Stretch^2 + par__epsconoid^2))	%% Eliminate negative stretches\n');
	fprintf(autolev_file, 'Uvec> = P_conO_conI> / LOI 												%% Unit vector from conoid origin to insertion\n');
	fprintf(autolev_file, 'Force(conO/conI, -par__conK*StretchPositive*Uvec> ) 				    %% Spring force\n\n');

	% Thoraco-humeral angles
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% Determine the YZY rotation angles from thorax to humerus\n');
	fprintf(autolev_file, '%% Apply the external torques on the axes of this joint coordinate system\n');
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% This is the symbolic form of an YZY rotation matrix:\n');
	fprintf(autolev_file, '%%  frames a,b,c,d\n');
	fprintf(autolev_file, '%%  variables q1,q2,q3\n');
	fprintf(autolev_file, '%%  simprot(a,b,2,q1);\n');
	fprintf(autolev_file, '%%  simprot(b,c,3,q2);\n');
	fprintf(autolev_file, '%%  simprot(c,d,2,q3);\n');
	fprintf(autolev_file, '%%  a_d\n');
	fprintf(autolev_file, '%% Result[1,1] = COS(q1)*COS(q2)*COS(q3) - SIN(q1)*SIN(q3)\n');
	fprintf(autolev_file, '%% Result[1,2] = -SIN(q2)*COS(q1)\n');
	fprintf(autolev_file, '%% Result[1,3] = SIN(q1)*COS(q3) + SIN(q3)*COS(q1)*COS(q2)\n');
	fprintf(autolev_file, '%% Result[2,1] = SIN(q2)*COS(q3)\n');
	fprintf(autolev_file, '%% Result[2,2] = COS(q2)\n');
	fprintf(autolev_file, '%% Result[2,3] = SIN(q2)*SIN(q3)\n');
	fprintf(autolev_file, '%% Result[3,1] = -SIN(q3)*COS(q1) - SIN(q1)*COS(q2)*COS(q3)\n');
	fprintf(autolev_file, '%% Result[3,2] = SIN(q1)*SIN(q2)\n');
	fprintf(autolev_file, '%% Result[3,3] = COS(q1)*COS(q3) - SIN(q1)*SIN(q3)*COS(q2)\n\n');

	fprintf(autolev_file, '%% so we can solve the YZY angles as follows (if we are not in gimbal lock where sin(q2) is zero)\n');
	fprintf(autolev_file, 'qTH = [														  &\n');
	fprintf(autolev_file, '	sym_atan2(thorax_humerusr[3,2] , (-thorax_humerusr[1,2])) 	; &		%% plane of elevation, between -PI and PI\n');
	fprintf(autolev_file, '	acos(thorax_humerusr[2,2])								; &		%% elevation angle, between 0 and PI\n');
	fprintf(autolev_file, '	sym_atan2(thorax_humerusr[2,3] , thorax_humerusr[2,1])	  &		%% humerus internal rotation, between -PI and PI\n');
	fprintf(autolev_file, '	]	\n');
	fprintf(autolev_file, 'encode qTH\n\n');

	fprintf(autolev_file, '%% calculate the floating axis (where we need to apply the Z moment)\n');
	fprintf(autolev_file, 'floatingaxis> = cross(thorax2>,humerusr2>)\n');
	fprintf(autolev_file, 'floatingaxis> = floatingaxis>/mag(floatingaxis>)		%% normalize to unit vector\n\n');

	fprintf(autolev_file, '%% apply an external moment on humerus that is equivalent to the JCS moments that we wanted\n');
	fprintf(autolev_file, 'variables MTHy, MTHz, MTHyy								%% the three JCS moments between thorax and humerus\n');
	fprintf(autolev_file, 'variables MHx, MHy, MHz									%% the XYZ component of the (equivalent) external moment\n');
	fprintf(autolev_file, 'MH> = vector(thorax, MHx, MHy, MHz)\n\n');

	fprintf(autolev_file, '%% equations that say that JCS moments are projections of the external moment:\n');
	fprintf(autolev_file, 'equation[1] = dot(MH>, thorax2>) 		- MTHy\n');
	fprintf(autolev_file, 'equation[2] = dot(MH>, floatingaxis>) 	- MTHz\n');
	fprintf(autolev_file, 'equation[3] = dot(MH>, humerusr2>) 	- MTHyy\n\n');

	fprintf(autolev_file, '%% now we can solve the external moment (will fail in gimbal lock!):\n');
	fprintf(autolev_file, 'solve(equation, MHx, MHy, MHz)\n\n');

	fprintf(autolev_file, '%% and we apply this external moment to the humerus\n');
	fprintf(autolev_file, 'Torque(thorax/humerusr, MH>)		\n\n');

	% Scapulothoracic contact forces
	fprintf(autolev_file, '%%--------------------------------------------------\n');
	fprintf(autolev_file, '%% Scapulothoracic contact forces on TS and AI\n');
	fprintf(autolev_file, '%%--------------------------------------------------\n');

	scap_points = {'TS';'AI'};
	fprintf(autolev_file, 'Constants par__Mx,par__My,par__Mz                        %% center of the thorax ellipsoid\n');
	fprintf(autolev_file, 'Constants par__Ax,par__Ay,par__Az                        %% radii of the thorax ellipsoid\n');
	fprintf(autolev_file, 'Constants par__epscontact, par__kcontact             %% parameters for scapula-thorax contact model\n\n');
	fprintf(autolev_file, 'Points CenterEllipsoid\n');
	fprintf(autolev_file, 'P_thoraxJoint_CenterEllipsoid> = Vector(thorax, par__Mx, par__My, par__Mz)\n\n');

	fprintf(autolev_file, '%% compute the coordinates of the ellipsoid center relative to thorax\n');
	fprintf(autolev_file, 'PxCE = dot(P_thoraxO_CenterEllipsoid>, thorax1>)\n');
	fprintf(autolev_file, 'PyCE = dot(P_thoraxO_CenterEllipsoid>, thorax2>)\n');
	fprintf(autolev_file, 'PzCE = dot(P_thoraxO_CenterEllipsoid>, thorax3>)\n\n');



	for iscap=1:length(scap_points)
		fprintf(autolev_file, '%% Create the point %s and its location x,y,z on the scapula\n',scap_points{iscap});
		fprintf(autolev_file, 'Points %s													\n',scap_points{iscap});
		fprintf(autolev_file, 'Constants par__%sx, par__%sy, par__%sz\n',scap_points{iscap},scap_points{iscap},scap_points{iscap});
		fprintf(autolev_file, 'P_scapularJoint_%s> = Vector(scapular, par__%sx, par__%sy, par__%sz)\n',scap_points{iscap},scap_points{iscap},scap_points{iscap},scap_points{iscap});
		fprintf(autolev_file, 'V2pts(%s, scapular, scapularJoint, %s)			%% Autolev will need velocity of contact point in equations of motion\n\n',newtonian, scap_points{iscap});

		fprintf(autolev_file, '%% compute its global coordinates (relative to thorax)\n');
		fprintf(autolev_file, 'Px%s = dot(P_thoraxO_%s>, thorax1>)\n',scap_points{iscap},scap_points{iscap});
		fprintf(autolev_file, 'Py%s = dot(P_thoraxO_%s>, thorax2>)\n',scap_points{iscap},scap_points{iscap});
		fprintf(autolev_file, 'Pz%s = dot(P_thoraxO_%s>, thorax3>)\n\n',scap_points{iscap},scap_points{iscap});

		fprintf(autolev_file, '%% computes the "distance" to thorax surface\n');
		fprintf(autolev_file, 'F%s = ((Px%s-PxCE)/par__Ax)^2 + ((Py%s-PyCE)/par__Ay)^2 + ((Pz%s-PzCE)/par__Az)^2 - 1\n\n',scap_points{iscap},scap_points{iscap},scap_points{iscap},scap_points{iscap});

		fprintf(autolev_file, '%% attenuate F when it is positive\n');
		fprintf(autolev_file, 'Fminus%s = 0.5*(F%s - sqrt(F%s^2 + par__epscontact^2))\n\n',scap_points{iscap},scap_points{iscap},scap_points{iscap});

		fprintf(autolev_file, '%% compute the forces						\n');
		fprintf(autolev_file, 'Fx%s = -par__kcontact*Px%s*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Ax^2*Fminus%s\n',scap_points{iscap},scap_points{iscap},scap_points{iscap});
		fprintf(autolev_file, 'Fy%s = -par__kcontact*Py%s*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Ay^2*Fminus%s\n',scap_points{iscap},scap_points{iscap},scap_points{iscap});
		fprintf(autolev_file, 'Fz%s = -par__kcontact*Pz%s*(par__Ax^2 + par__Ay^2 + par__Az^2)/par__Az^2*Fminus%s\n\n',scap_points{iscap},scap_points{iscap},scap_points{iscap});

		fprintf(autolev_file, '%% apply the force to the point %s\n',scap_points{iscap});
		fprintf(autolev_file, 'Force_%s> += Fx%s*thorax1> + Fy%s*thorax2> + Fz%s*thorax3>\n\n',scap_points{iscap},scap_points{iscap},scap_points{iscap},scap_points{iscap});
	end

	fprintf(autolev_file, '%% Also encode those so we have them in the MEX function\n');
	fprintf(autolev_file, 'F_SCAP = [FxTS , FyTS, FzTS; FxAI , FyAI, FzAI];\n');
	fprintf(autolev_file, '%% Also encode the TS and AI coordinates and the "distance" to thorax surface\n');
	fprintf(autolev_file, 'dist_SCAP = [FTS; FAI];\n');				
	fprintf(autolev_file, 'pos_SCAP = [PxTS , PyTS, PzTS; PxAI , PyAI, PzAI];\n');				
	fprintf(autolev_file, 'encode F_SCAP, dist_SCAP, pos_SCAP\n\n');


	% GH net reaction force
	fprintf(autolev_file, '%%----------------------------------------------------------\n');
	fprintf(autolev_file, '%% GH net reaction force acting on the Scapula\n');
	fprintf(autolev_file, '%%----------------------------------------------------------\n\n');

	fprintf(autolev_file, '%% compute net reaction force from weight and accelerations of arm segments\n');
	fprintf(autolev_file, '%% this will be the actual reaction force if the accelerations satisfy the equations of motion (Zero = 0)\n');
	fprintf(autolev_file, 'FGH> = 		par__humerusrMass * (-A_humerusrO_%s> + GravityAxis>)	&\n',newtonian);
	fprintf(autolev_file, '		+	par__ulnarMass    * (-A_ulnarO_%s>    + GravityAxis>)	&\n',newtonian);
	fprintf(autolev_file, '		+	par__radiusrMass  * (-A_radiusrO_%s>  + GravityAxis>)	&\n',newtonian);
	fprintf(autolev_file, '		+	par__handrMass    * (-A_handrO_%s>    + GravityAxis>)	\n\n',newtonian);

	fprintf(autolev_file, '%% express the result as XYZ force components in Scapula reference frame\n');
	fprintf(autolev_file, 'F_GH = [	dot(FGH>, scapular1>) ; dot(FGH>, scapular2>) ; dot(FGH>, scapular3>) ]\n\n');

	fprintf(autolev_file, '%% we need to report F_GH back to the simulation program to help detect shoulder dislocation\n');
	fprintf(autolev_file, '%% Note that this is the net intersegmental force, not accounting for muscles.  Muscle contributions\n');
	fprintf(autolev_file, '%% are done in das3.c\n');
	fprintf(autolev_file, 'encode F_GH		\n\n');				


	% Equations of motion
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% Generate equations of motion\n');
	fprintf(autolev_file, '%%--------------------------------------------------------------------\n\n');

	fprintf(autolev_file, 'Zero = Fr() + FrStar()\n\n');

	% Visualization
	fprintf(autolev_file, '%%---------------------------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% For visualization, we export position and orientation of the segments in a matrix.\n');
	fprintf(autolev_file, '%% This output is only used for model testing in Matlab.  Required processor time is\n');
	fprintf(autolev_file, '%% very small compared to dynamics, so we just generate this always.\n');
	fprintf(autolev_file, '%%---------------------------------------------------------------------------------------\n\n');

	fprintf(autolev_file, 'vis = [dot(P_%sO_%sJoint>, %s1>), dot(P_%sO_%sJoint>, %s2>), dot(P_%sO_%sJoint>, %s3>), rows(%s_%s,1), rows(%s_%s,2), rows(%s_%s,3);&\n',newtonian,model.segments{allsegs(1)}.name,newtonian,newtonian,model.segments{allsegs(1)}.name,newtonian,newtonian,model.segments{allsegs(1)}.name,newtonian,newtonian,model.segments{allsegs(1)}.name,newtonian,model.segments{allsegs(1)}.name,newtonian,model.segments{allsegs(1)}.name);
	for iseg=2:length(allsegs)-1
		fprintf(autolev_file, '		dot(P_%sO_%sJoint>, %s1>), dot(P_%sO_%sJoint>, %s2>), dot(P_%sO_%sJoint>, %s3>), rows(%s_%s,1), rows(%s_%s,2), rows(%s_%s,3);&\n',newtonian,model.segments{allsegs(iseg)}.name,newtonian,newtonian,model.segments{allsegs(iseg)}.name,newtonian,newtonian,model.segments{allsegs(iseg)}.name,newtonian,newtonian,model.segments{allsegs(iseg)}.name,newtonian,model.segments{allsegs(iseg)}.name,newtonian,model.segments{allsegs(iseg)}.name);
	end
	fprintf(autolev_file, '		dot(P_%sO_%sJoint>, %s1>), dot(P_%sO_%sJoint>, %s2>), dot(P_%sO_%sJoint>, %s3>), rows(%s_%s,1), rows(%s_%s,2), rows(%s_%s,3)];\n\n',newtonian,model.segments{allsegs(end)}.name,newtonian,newtonian,model.segments{allsegs(end)}.name,newtonian,newtonian,model.segments{allsegs(end)}.name,newtonian,newtonian,model.segments{allsegs(end)}.name,newtonian,model.segments{allsegs(end)}.name,newtonian,model.segments{allsegs(end)}.name);
	fprintf(autolev_file, 'encode vis\n\n');

	% Implicit dynamics
	fprintf(autolev_file, '%%-----------------------------------------------------------------------------------------------------\n');
	fprintf(autolev_file, '%% For implicit dynamics: generate symbolic expressions for Zero and Jacobians d(Zero)/d(q,qd,qdd)\n');
	fprintf(autolev_file, '%%-----------------------------------------------------------------------------------------------------\n\n');

	fprintf(autolev_file, '%% q are the generalized coordinates\n');
	fprintf(autolev_file, 'q = [%s, &\n',model.dofs{1}.name);
	for idof = 2:model.nDofs-1
		fprintf(autolev_file, '%s, &\n',model.dofs{idof}.name);
	end
	fprintf(autolev_file, '%s]\n\n',model.dofs{model.nDofs}.name);

	fprintf(autolev_file, '%% generalized velocities\n');
	fprintf(autolev_file, 'qd = [%s'', &\n',model.dofs{1}.name);
	for idof = 2:model.nDofs-1
		fprintf(autolev_file, '%s'', &\n',model.dofs{idof}.name);
	end
	fprintf(autolev_file, '%s'']\n\n',model.dofs{model.nDofs}.name);

	fprintf(autolev_file, '%% generalized accelerations\n');
	fprintf(autolev_file, 'qdd = [%s'''', &\n',model.dofs{1}.name);
	for idof = 2:model.nDofs-1
		fprintf(autolev_file, '%s'''', &\n',model.dofs{idof}.name);
	end
	fprintf(autolev_file, '%s'''']\n\n',model.dofs{model.nDofs}.name);

	fprintf(autolev_file, '%% generate expressions for implicit equation of motion and its Jacobians\n');
	fprintf(autolev_file, 'dz_dq = ZEE(D(Zero,q))\n');
	fprintf(autolev_file, 'dz_dqd = ZEE(D(Zero,qd))\n');
	fprintf(autolev_file, 'dz_dqdd = ZEE(D(Zero,qdd))\n');
	fprintf(autolev_file, 'encode Zero, dz_dq, dz_dqd, dz_dqdd\n\n');

	fprintf(autolev_file, '%% write all Encoded expressions to C file\n');
	fprintf(autolev_file, 'Code Algebraic() das3_al_raw.c\n\n');

	fprintf(autolev_file, 'EXIT\n');

	fclose(autolev_file);

end