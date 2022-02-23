function model = das3_readosim(osimfile,musclepolyfile,GHpolyfile)
% das3_readosim initializes a Matlab struct from an Opensim model file
%
% osimfile: the Opensim model file from which to initialize
% musclepolyfile: a .mat file with polynomials for muscle-tendon lengths
% computed by das3_polynomials.m
% GHpolyfile: a .mat file with polynomials for GH force vectors
% computed by das3_polynomials.m
%
% If this Matlab struct is used to create the Autolev file or the MEX file,
% muscle path  and GH force vector information is not required
% so "musclepolyfile" and "GHpolyfile" are optional inputs
%
% Dimitra Blana, October 2014
% based on gait3d_readosim.m by Ton van den Bogert %
%
% Updated for OpenSim 4.0 by Derek Wolf, October 2019


% import OpenSim namespace
import org.opensim.modeling.*;

% load the model
Mod = Model(osimfile);

% initialize the system and get the initial state
state = Mod.initSystem();

% equilibrate all muscles
Mod.equilibrateMuscles(state);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basic model parameters
model.name = char(Mod.getName());
model.osimfile = char(Mod.getInputFileName());
file = dir(osimfile);
model.modified = file.date;

% gravity
tmpVec3f = Mod.getGravity();
model.gravity = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of joints
JointSet = Mod.getJointSet();
model.nJoints = JointSet.getSize();

% loop through joints
for ijoints = 1:model.nJoints
    currentJoint = JointSet.get(ijoints-1);
    
    % basics
    model.joints{ijoints}.name = fixname(char(currentJoint.getName()));
    
    % location and orientation in parent frame
    parent_frame=currentJoint.get_frames(0);
    tmpVec3f=parent_frame.get_translation(); % location
    model.joints{ijoints}.location = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    
    tmpVec3f=parent_frame.get_orientation(); % orientation
    model.joints{ijoints}.orientation = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    
    % name of parent segment this joint is attached to
    model.joints{ijoints}.parent_segment=fixname(char(parent_frame.findBaseFrame()));
    
    % name of segment this joint belongs to
    child_frame=currentJoint.get_frames(1);
    model.joints{ijoints}.segment=fixname(char(child_frame.findBaseFrame()));
    
    % check if orientation and translation of child frame is 0
    tmpVec3f=child_frame.get_translation();
    child_translation=[tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    tmpVec3f=child_frame.get_orientation(); % orientation
    child_orientation = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    
    if ~isequal(child_translation,[0 0 0]) || ~isequal(child_orientation,[0 0 0])
        error('Child frame translation or orientation is not [0,0,0]')
    end
    
    % spatial transforms
    currentJoint = CustomJoint.safeDownCast(currentJoint);
    % if this is a weld joint, there are no spatial transforms
    if isempty(currentJoint), continue; end
    transforms = currentJoint.getSpatialTransform();
    
    % rotations
    r1 = transforms.getTransformAxis(0);
    tmpVec3f = r1.getAxis();
    model.joints{ijoints}.r1_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    model.joints{ijoints}.r1_coefs = osimfunction2coefs(r1.getFunction());
    coord_names = r1.getCoordinateNamesInArray;
    for i_coords = 1:coord_names.getSize
        model.joints{ijoints}.r1_coordinates{i_coords} = fixname(char(coord_names.getitem(i_coords-1)));
    end
    
    r2 = transforms.getTransformAxis(1);
    tmpVec3f = r2.getAxis();
    model.joints{ijoints}.r2_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    model.joints{ijoints}.r2_coefs = osimfunction2coefs(r2.getFunction());
    coord_names = r2.getCoordinateNamesInArray;
    for i_coords = 1:coord_names.getSize
        model.joints{ijoints}.r2_coordinates{i_coords} = fixname(char(coord_names.getitem(i_coords-1)));
    end
    
    r3 = transforms.getTransformAxis(2);
    tmpVec3f = r3.getAxis();
    model.joints{ijoints}.r3_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    model.joints{ijoints}.r3_coefs = osimfunction2coefs(r3.getFunction());
    coord_names = r3.getCoordinateNamesInArray;
    for i_coords = 1:coord_names.getSize
        model.joints{ijoints}.r3_coordinates{i_coords} = fixname(char(coord_names.getitem(i_coords-1)));
    end
    
    % translations
    t1 = transforms.getTransformAxis(3);
    tmpVec3f = t1.getAxis();
    model.joints{ijoints}.t1_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    model.joints{ijoints}.t1_coefs = osimfunction2coefs(t1.getFunction());
    coord_names = t1.getCoordinateNamesInArray;
    for i_coords = 1:coord_names.getSize
        model.joints{ijoints}.t1_coordinates{i_coords} = fixname(char(coord_names.getitem(i_coords-1)));
    end
    
    t2 = transforms.getTransformAxis(4);
    tmpVec3f = t2.getAxis();
    model.joints{ijoints}.t2_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    model.joints{ijoints}.t2_coefs = osimfunction2coefs(t2.getFunction());
    coord_names = t2.getCoordinateNamesInArray;
    for i_coords = 1:coord_names.getSize
        model.joints{ijoints}.t2_coordinates{i_coords} = fixname(char(coord_names.getitem(i_coords-1)));
    end
    
    t3 = transforms.getTransformAxis(5);
    tmpVec3f = t3.getAxis();
    model.joints{ijoints}.t3_axis = [tmpVec3f.get(0) tmpVec3f.get(1) tmpVec3f.get(2)];
    model.joints{ijoints}.t3_coefs = osimfunction2coefs(t3.getFunction());
    coord_names = t3.getCoordinateNamesInArray;
    for i_coords = 1:coord_names.getSize
        model.joints{ijoints}.t3_coordinates{i_coords} = fixname(char(coord_names.getitem(i_coords-1)));
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEGMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the set of segments
BodySet = Mod.getBodySet();
model.nSegments = BodySet.getSize();
segment_names = cell(model.nSegments,1);

% loop through segments
for isegment = 1:model.nSegments
    currentSegment = BodySet.get(isegment-1);
    
    % name and mass
    model.segments{isegment}.name = fixname(char(currentSegment.getName()));
    model.segments{isegment}.mass = currentSegment.getMass();
    
    segment_names{isegment} = char(currentSegment.getName());
    
    % mass center
    tmpVec3=currentSegment.getMassCenter();
    model.segments{isegment}.mass_center = [tmpVec3.get(0) tmpVec3.get(1) tmpVec3.get(2)];
    
    % inertia matrix
    tmpVec6=currentSegment.get_inertia();
    model.segments{isegment}.inertia = zeros(3,3);
    model.segments{isegment}.inertia(1,1)=tmpVec6.get(0); % Ixx
    model.segments{isegment}.inertia(2,2)=tmpVec6.get(1); % Iyy
    model.segments{isegment}.inertia(3,3)=tmpVec6.get(2); % Izz
    model.segments{isegment}.inertia(1,2)=tmpVec6.get(3); % Ixy
    model.segments{isegment}.inertia(1,3)=tmpVec6.get(4); % Ixz
    model.segments{isegment}.inertia(2,3)=tmpVec6.get(5); % Iyz
    model.segments{isegment}.inertia(2,1)=tmpVec6.get(3); % Ixy
    model.segments{isegment}.inertia(3,1)=tmpVec6.get(4); % Ixz
    model.segments{isegment}.inertia(3,2)=tmpVec6.get(5); % Iyz
    
    
    % name of parent joint this segment is attached to
    model.segments{isegment}.parent_joint='';
    for ijoint=1:model.nJoints
        if strcmp(model.joints{ijoint}.segment,model.segments{isegment}.name)
            model.segments{isegment}.parent_joint=model.joints{ijoint}.name;
            break;
        elseif ijoint==model.nJoints % error check if no joint name was found
            error('No matches for segment name were found');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of DOFs
CoordSet = Mod.getCoordinateSet();
model.nDofs = CoordSet.getSize();

counter=0;

% loop through dofs
for idof = 1:model.nDofs
    
    currentDof = CoordSet.get(idof-1);
    
    % remove constrained dofs from the model
    if currentDof.isConstrained(state) && ~currentDof.getLocked(state), continue; end
    
    counter=counter+1;
    % make a note of locked dofs from the model
    if currentDof.getLocked(state)
        model.dofs{counter}.locked = 1;
    else
        model.dofs{counter}.locked = 0;
    end
    
    % save indeces of dofs
    model.dof_indeces(counter) = idof-1;
    
    % basics
    model.dofs{counter}.name = fixname(char(currentDof.getName()));
    model.dofs{counter}.osim_name = char(currentDof.getName());
    
    % range of motion in radians
    %     model.dofs{counter}.range(1) = currentDof.getRangeMin()+2/180*pi;
    %     model.dofs{counter}.range(2) = currentDof.getRangeMax()-2/180*pi;
    model.dofs{counter}.range(1) = currentDof.getRangeMin();
    model.dofs{counter}.range(2) = currentDof.getRangeMax();
end

model.nDofs = counter;

dof_names = cell(model.nDofs,1);
for idof = 1:model.nDofs
    dof_names{idof} = model.dofs{idof}.name;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of constraints
ConstraintSet = Mod.getConstraintSet();
model.nConstraints = ConstraintSet.getSize();

counter=0;

% loop through constraints
for iconstr = 1:model.nConstraints
    
    currentConstraint = ConstraintSet.get(iconstr-1);
    
    % remove disabled constraints from the model
    if currentConstraint.isDisabled(state), continue; end
    counter=counter+1;
    
    % These should be CoordinateCouplerConstraints
    currentConstraint = CoordinateCouplerConstraint.safeDownCast(currentConstraint);
    
    model.constraints{counter}.dependent_dof = fixname(char(currentConstraint.getDependentCoordinateName));
    ind_coords = currentConstraint.getIndependentCoordinateNames;
    for i_ind_coords = 1:ind_coords.getSize
        model.constraints{counter}.independent_dof{i_ind_coords} = fixname(char(ind_coords.getitem(i_ind_coords-1)));
    end
    model.constraints{counter}.coefs = osimfunction2coefs(currentConstraint.getFunction());
end

model.nConstraints = counter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUSCLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the set of muscles
MuscleSet = Mod.getMuscles();
model.nMus = MuscleSet.getSize();

% get the metabolic probe, that contains muscle masses
Probe = Umberger2010MuscleMetabolicsProbe.safeDownCast(Mod.getProbeSet().get('metabolics'));

counter=0;

% loop through muscles
for imus = 1:model.nMus
    currentMuscle = MuscleSet.get(imus-1);
    
    
    counter=counter+1;
    
    % basic muscle properties
    model.muscles{imus}.name = fixname(char(currentMuscle.getName()));
    model.muscles{imus}.osim_name = char(currentMuscle.getName());
    model.muscles{imus}.mass = Probe.getMuscleMass(model.muscles{imus}.osim_name);
    model.muscles{imus}.fmax = currentMuscle.getMaxIsometricForce();
    model.muscles{imus}.lceopt = currentMuscle.getOptimalFiberLength();
    model.muscles{imus}.lslack = currentMuscle.getTendonSlackLength();
    model.muscles{imus}.pennopt = currentMuscle.getPennationAngleAtOptimalFiberLength();
    model.muscles{imus}.vmax = currentMuscle.getMaxContractionVelocity();
    
    act1 = str2double(char(currentMuscle.getPropertyByName('activation1')));
    act2 = str2double(char(currentMuscle.getPropertyByName('activation2')));
    
    % From Lisa Schutte's dissertation, appendices 2 and 3.
    % "Using musculoskeletal models to explore strategies for improving
    % performance in electrical stimulation-induced leg cycle ergometry"
    % LM Schutte - 1992 - Stanford University
    % k1 = activation1
    % k2 = activation2
    % t_act = 1/(k1+k2)
    % t_deact = 1/k2
    model.muscles{imus}.tact = 1/(act1+act2);
    model.muscles{imus}.tdeact = 1/act2;
    
    currentprop = currentMuscle.getPropertyByName('active_force_length_curve');
    active_fl_o = currentprop.getValueAsObject;
    active_fl_spline = SimmSpline.safeDownCast(active_fl_o);
    num_values = active_fl_spline.getSize;
    model.muscles{imus}.Xval_active_fl = zeros(num_values,1);
    model.muscles{imus}.Yval_active_fl = zeros(num_values,1);
    for i=1:num_values
        model.muscles{imus}.Xval_active_fl(i) = active_fl_spline.getX(i-1);
        model.muscles{imus}.Yval_active_fl(i) = active_fl_spline.getY(i-1);
    end
    
    currentprop = currentMuscle.getPropertyByName('passive_force_length_curve');
    passive_fl_o = currentprop.getValueAsObject;
    passive_fl_spline = SimmSpline.safeDownCast(passive_fl_o);
    num_values = passive_fl_spline.getSize;
    model.muscles{imus}.Xval_passive_fl = zeros(num_values,1);
    model.muscles{imus}.Yval_passive_fl = zeros(num_values,1);
    for i=1:num_values
        model.muscles{imus}.Xval_passive_fl(i) = passive_fl_spline.getX(i-1);
        model.muscles{imus}.Yval_passive_fl(i) = passive_fl_spline.getY(i-1);
    end
    model.muscles{imus}.PEEslack = model.muscles{imus}.Xval_passive_fl(4);
    
    % get muscle path and DOFs this muscle crosses
    muspath = currentMuscle.getGeometryPath();
    PtSet = muspath.getPathPointSet();
    nPts = PtSet.getSize();
    
    origin_segment = PtSet.get(0).getBody();
    origin_segment_index = find(strcmp(origin_segment,segment_names));
    insertion_segment = PtSet.get(nPts-1).getBody();
    insertion_segment_index = find(strcmp(insertion_segment,segment_names));
    
    % if the origin is more distal than the insertion, flip the
    % segments
    if origin_segment_index>insertion_segment_index
        helpseg = origin_segment_index;
        origin_segment_index = insertion_segment_index;
        insertion_segment_index = helpseg;
    end
    
    current_segment_index = insertion_segment_index;
    
    dof_count = 0;
    dof_list = cell(model.nDofs,1);
    dof_indeces = zeros(model.nDofs,1);
    
    while (current_segment_index ~= origin_segment_index)
        
        % find current joint
        current_joint_name=model.segments{current_segment_index}.parent_joint;
        for ijoint=1:model.nJoints
            if strcmp(model.joints{ijoint}.name,current_joint_name)
                current_joint = JointSet.get(ijoint-1);
                break;
            elseif ijoint==model.nJoints % error check if no joint name was found
                error('No matches for segment name were found');
            end
        end
        
        ndofs = current_joint.numCoordinates();
        
        for idofs=1:ndofs
            dof_name = fixname(char(current_joint.get_coordinates(idofs-1)));
            dof_index = find(strcmp(dof_name,dof_names));
            if ~isempty(dof_index)
                dof_count=dof_count+1;
                dof_list{dof_count} = dof_name;
                dof_indeces(dof_count) = dof_index;
            end
        end
                
        parent_frame=current_joint.get_frames(0);
        parentBody=parent_frame.findBaseFrame();
        current_segment_index=find(strcmp(parentBody,segment_names));
    end
    
    [sorted_dofs, sorted_indeces] = sort(dof_indeces(1:dof_count));
    dof_list = dof_list(1:dof_count);
    model.muscles{imus}.dof_count = dof_count;
    model.muscles{imus}.dof_indeces = sorted_dofs;
    model.muscles{imus}.dof_names = dof_list(sorted_indeces);
end

model.nMus = counter;

% add polynomial muscle path information
% this is not required for the Autolev and MEX files, so it is optional
if nargin>1
    poly_mus = load(musclepolyfile);
    allmusi = zeros(length(poly_mus.mus_model));
    index=1;
    for imus=1:length(poly_mus.mus_model)
        current_mus = poly_mus.mus_model{imus};
        for iimus=1:model.nMus
            if strcmp(current_mus.name,model.muscles{iimus}.osim_name) || strcmp(current_mus.name,model.muscles{iimus}.name)
                model.muscles{iimus}.lparam_count = current_mus.num_lparams;
                model.muscles{iimus}.lparams = current_mus.lparams;
                model.muscles{iimus}.lcoefs = current_mus.lcoef;
                allmusi(index) = iimus;
                index=index+1;
                break;
            end
        end
    end
    
    if index<model.nMus
        miss_mus = setdiff(1:model.nMus,allmusi);
        miss_names = model.muscles{miss_mus(1)}.osim_name;
        for imus=2:length(miss_mus)
            miss_names = [miss_names ', ' model.muscles{miss_mus(imus)}.osim_name];
        end
        error(['No polynomial muscle path found for ', miss_names '. Please run das3_polynomials.m']);
    end
end

% add polynomial GH force vector information to the muscles that cross GH
% this is not required for the Autolev and MEX files, so it is optional
if nargin>2
    poly_GH = load(GHpolyfile);
    for imus=1:model.nMus
        current_mus = poly_GH.GH_model{imus};
        if isempty(current_mus)
            model.muscles{imus}.crossesGH = 0;
            continue;
        end
        model.muscles{imus}.crossesGH = 1;
        model.muscles{imus}.xparam_count = current_mus.xparam_count;
        model.muscles{imus}.xparams = current_mus.xparams;
        model.muscles{imus}.xcoefs = current_mus.xcoef;
        
        model.muscles{imus}.yparam_count = current_mus.yparam_count;
        model.muscles{imus}.yparams = current_mus.yparams;
        model.muscles{imus}.ycoefs = current_mus.ycoef;
        
        model.muscles{imus}.zparam_count = current_mus.zparam_count;
        model.muscles{imus}.zparams = current_mus.zparams;
        model.muscles{imus}.zcoefs = current_mus.zcoef;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARKERS - Not tested for shoulder model!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does not work with osim4 and shoulder model

%{
% get the set of markers
MarkerSet = Mod.getMarkerSet();
model.nMarkers = MarkerSet.getSize();

% loop through all markers
for imarkers = 1:model.nMarkers
    currentMarker = MarkerSet.get(imarkers-1);

    
    % basic properties
    model.markers{imarkers}.name = fixname(char(currentMarker.getName()));
    model.markers{imarkers}.segment = fixname(char(currentMarker.getBodyName()));
    model.markers{imarkers}.fixed = currentMarker.getFixed();
    
    tmpVec3 = currentMarker.getOffset();
    model.markers{imarkers}.position = [tmpVec3.get(0) tmpVec3.get(1) tmpVec3.get(2)];
    
    % find index
    for ibodys = 1:model.nSegments
        if strcmp(model.segments{ibodys}.name, model.markers{imarkers}.segment)
            model.markers{imarkers}.segmentindex = ibodys;
        end
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONOID LIGAMENT FORCE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The stiffness and length of the conoid ligament have been chosen so that
% the force-displacement curve matches human cadaver measurements described in Harris et al.
% Am J Sports Med, vol.28, 2000, pp.103-8.

model.conoid_eps	= 0.001;	% epsilon for conoid force-length model
model.conoid_stiffness = 80000;	% stiffness of conoid
model.conoid_length = 0.0174;	% length (m) of conoid ligament
model.conoid_origin =  [0.1365, 0.0206, 0.0136];	% coordinates of conoid origin in clavicle frame
model.conoid_insertion = [-0.0536, -0.0009, -0.0266];	% coordinates of conoid insertion in scapula frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCAPULA-THORAX CONTACT MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.scap_thorax_eps 	= 0.01;		% epsilon for scapula-thorax contact model ~ 0.01 is recommended
model.scap_thorax_k   	= 2e4;		% stiffness k (N/m) in scapula-thorax contact model (2e4 is recommended)
% these are the projections of TS and AI onto the thorax ellipsoid, not the actual TS and AI
% model.TSprojection = [-0.1305, -0.0181, -0.0046];
% model.AIprojection = [-0.1267, -0.1246, 0.0052];
model.TSprojection = [-0.1274   -0.0228   -0.0346]; % coordinates of TS projection in scapula frame
model.AIprojection = [-0.1257   -0.1223   -0.0275]; % coordinates of AI projection in scapula frame

% get thorax ellipsoid information from thorax wrapping object
for isegment = 1:model.nSegments
    currentSegment = BodySet.get(isegment-1);
    
    if strcmp(currentSegment.getName(),'thorax')
        ellipsoid =  WrapEllipsoid.safeDownCast(currentSegment.getWrapObjectSet.get(0));
        tmpVec3 = ellipsoid.getRadii();
        model.thorax_radii = [tmpVec3.get(0) tmpVec3.get(1) tmpVec3.get(2)]; % radii of the thorax ellipsoid
        trans = char(ellipsoid.getPropertyByName('translation'));
        model.thorax_center = sscanf(trans,'(%f %f %f)');
    end
end


end

function name = fixname(orig_name)
% no underscores in names
name = strrep(orig_name, '_', '');
while ~isletter(name(1))
    if isstrprop(name(1),'digit')
        % move numbers to the end of the name
        name = [name(2:end) name(1)];
    else
        % remove underscores and other non-alphanumeric characters
        name = name(2:end);
    end
    if isempty(name)
        error(['Name ',orig_name ' does not include any alphabetic characters. Please rename.']);
    end
end

max_name_length = 15;
rem_lets = ['-';'u';'a';'e';'i';'o';'s';'t';'l';'r';'n';'h';'d'];
index=1;
while length(name)>max_name_length && index<=length(rem_lets)
    name = strrep(name, rem_lets(index), '');
    index = index+1;
end
if length(name)>max_name_length
    error(['Name ',orig_name ' is too long. Please rename.']);
end
end

function coefs = osimfunction2coefs(f)
% reads in a function from the osim file and returns two possible outputs:
% a vector of two elements, the coefficients of a line (for linear functions)
% OR a vector of six elements, describing a spline in the correct format for Autolev

import org.opensim.modeling.*;
coefs=[];

%%%
% We want to take the scaling into account https://simtk.org/api_docs/opensim/api_docs/classOpenSim_1_1MultiplierFunction.html

% first we define a default scaling factor of 1 that will be applied to all the values :
scalingFactor = 1;

% Now check if it's scaled (if we have a multiplier function)
fx2=MultiplierFunction.safeDownCast(f);

if ~isempty(fx2) % if not empty

	% This is a multiplier function
    scalingFactor = fx2.getScale; % update the scalingFactor
    fx3 = fx2.getFunction; % extract the function inside the multiplierfunction
    f = fx3; % replace (f) to continue the function normally

end
%%%

fx=Constant.safeDownCast(f);
if ~isempty(fx)
    coefs = [fx.getValue() 0];
    return;
end

fx=LinearFunction.safeDownCast(f);
if ~isempty(fx)
    coefs = [fx.getIntercept() fx.getSlope()];
    return;
end

fx=SimmSpline.safeDownCast(f);
if ~isempty(fx)
    s = fx.getSize();
    
    knots=zeros(2,s);
    for i=1:s
        knots(1,i) = fx.getX(i-1);
        knots(2,i) = fx.getY(i-1);
    end
    
    if s==2
        % this is a line
        p = polyfit(knots(1,:),knots(2,:),1);
        coefs = [p(2) p(1)];
    elseif s==3
        % this is a 3-point spline
        % points where we want to sample the spline
        xx = knots(1,1):0.001:knots(1,end);
        % construct the natural cubic spline and differentiate
        yy = spline(knots(1,:),knots(2,:),xx);
        dydx = diff(yy)./diff(xx);
        % In Autolev, a cubic spline is done by specifying x,y
        % at the endpoints, and derivative dy/dx at the endpoints
        coefs = [knots(1,1),knots(1,end),knots(2,1),knots(2,end),dydx(1),dydx(end)];
    else
        % this is a spline with more than 3 points
        error('Splines with more than 3 points are not supported yet');
        
    end
    return;
end

if isempty(coefs)
    error('Function Type unknown');
end
end
