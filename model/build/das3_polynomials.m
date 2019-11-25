function das3_polynomials(osimfile,mydir,musclepolyfile,GHpolyfile)
% builds .mat files that contain polynomial approximations of all moment
% arms and GH lines of action for the das3 model
%
% osimfile: the Opensim model file from which to initialize
% musclepolyfile: the name of the .mat file with polynomials for muscle-tendon lengths
% GHpolyfile: the name of the .mat file with polynomials for GH lines of action
%
% Dimitra Blana, October 2014
% 
% Updated for OpenSim 4.0 by Derek Wolf, November 2019

% creates structure with information on the model
addpath('..');
model = das3_readosim(osimfile);
rmpath('..');
muscles = model.muscles;

import org.opensim.modeling.*
Mod = Model(osimfile);

% this is needed to get the GH lines of action in the scapular frame
SimEn = SimbodyEngine();
SimEn.connectSimbodyEngineToModel(Mod);
groundbody = Mod.getBodySet().get('thorax');
scapulabody = Mod.getBodySet().get('scapula_r');

default_dof = 0; % default value for all dofs
allvectxts = '';
jnt_values = cell(1,model.nDofs);
jnt_vals = cell(1,model.nDofs);

for ijnt=1:model.nDofs
    jnt_values{ijnt} = default_dof;
    jnt_vals{ijnt} = default_dof;
    if isempty(allvectxts)
        allvectxts = ['jnt_values{' num2str(ijnt) '}'];
    else
        allvectxts = [allvectxts ',jnt_values{' num2str(ijnt) '}'];
    end
end

shoulder_angles = DAS3_workspace('DAS2data');

dof_names = cell(model.nDofs,1);
for idof = 1:model.nDofs
    dof_names{idof} = model.dofs{idof}.name;
end

%% for each muscle...
for imus = 1:length(muscles)
%for imus = 1:0
    mus = muscles{imus};
    
    % give all dofs the default value
    for ijnt=1:model.nDofs
        jnt_values{ijnt} = default_dof;
    end
    
    % if the muscle only crosses the elbow, go through the range
    if mus.dof_indeces(1)>12
        for ijnt = 1:length(mus.dof_indeces)
            onedof = mus.dof_indeces(ijnt);
            lims = model.dofs{onedof}.range;
    %        values = lims(1):10:lims(2);  % every 10 degrees
            jnt_values{onedof} = lims(1):(lims(2)-lims(1))/4:lims(2);  % 5 steps
        end

        % alljnts contain all the combinations of joint angles for the range of
        % motion of this muscle
        % rows: # combinations of angles
        % columns: # angles
        try
        eval(['alljnts = combvec(' allvectxts ')'';']);
        catch err
            disp(err);
            disp(['Skipping muscle ',mus.name, ' ...sorry!']);
            continue;     
        end
    elseif mus.dof_indeces(end)<13
        % if the muscle only crosses the shoulder, use "shoulder_angles"
        alljnts = shoulder_angles;
    else
        % if the muscle crosses both shoulder and elbow, use
        % shoulder_angles + range of motion
        clear elbow
        lims = model.dofs{13}.range;
        elbow{1} = lims(1):(lims(2)-lims(1))/5:lims(2);  % 6 steps
        lims = model.dofs{14}.range;
        elbow{2} = lims(1):(lims(2)-lims(1))/5:lims(2);  % 6 steps

        % alljnts contain all the combinations of joint angles for the range of
        % motion of this muscle
        % rows: # combinations of angles
        % columns: # angles
        try
        alljnts = combvec(shoulder_angles(:,1:12)',elbow{1},elbow{2})';
        catch err
            disp(err);
            disp(['Skipping muscle ',mus.name, ' ...sorry!']);
            continue;     
        end
    end
    
    try
        % instead of moment arms directly from Opensim, calculate them using -dL/dq
        % if the muscle crosses GH, also get the lines of action
        if any(strcmp('GHy',mus.dof_names))
            [alllengths, allmomarms, allGHfvecs] = opensim_get_polyvalues(Mod, alljnts, mus.name, imus, mus.dof_indeces, SimEn, groundbody, scapulabody);
            save([mydir '\path_' mus.name],'alljnts','alllengths','allmomarms');        
            save([mydir '\GHfvec_' mus.name],'alljnts','allGHfvecs');        
            disp([mus.name, ' lengths, moment arms and GH force vectors saved.']);
        else
            [alllengths, allmomarms] = opensim_get_polyvalues(Mod, alljnts, mus.name, imus, mus.dof_indeces);
            save([mydir '\path_' mus.name],'alljnts','alllengths','allmomarms');        
            disp([mus.name, ' lengths and moment arms saved.']);
        end
        
        make_mot_file([mydir '\angles_' mus.name '.mot'],alljnts,dof_names);
        disp(['Opensim motion file for muscle ', mus.name, ' created.']);
    catch err
        disp(err);
        return;
    end

    clear mus alljnts alllengths allmomarms jnt_values all_dof_names
end

%% generate polynomial approximations of lengths
if nargin>2, musclepath_poly(muscles,mydir,musclepolyfile); end

%% generate polynomial approximations of GH lines of action
if nargin>3, GH_poly(muscles,mydir,GHpolyfile); end


%=============================================================================================
function [lengths, minusdLdq, GHfvecs] = opensim_get_polyvalues(Mod, angles, Mus, iMus, Dofs, SimEn, groundbody,scapulabody)
% This function calculates the length and -dL/dq of muscle "Mus" about dof set 
% "Dofs" at a given angle matrix "angles" of opensim model "Mod"
%
% "Angles" can be a vector (one hand position) or a matrix (one hand
% position per row)
% "Dofs" is a vector of indeces
%
% Adapted from opensim_get_momarm.m 
% Dimitra Blana, March 2012
%
% 28/3/2012: Use setValue (with active constraints) instead of set state
% 1/10/2014: Simplified how the muscles and dofs are accessed
%
% 11/25/19: Update by Derek Wolf: iMus is used to access the muscle and a
% check is used to determine if the name (with no underscores) matches the
% name in Mus

import org.opensim.modeling.*
% initialize the system to get the initial state
state = Mod.initSystem;

Mod.equilibrateMuscles(state);

% If we only want one moment arm, put it in a cell
CoordSet = Mod.getCoordinateSet();
num_request_dofs = length(Dofs);

% get the muscle
% currentMuscle = Mod.getMuscles().get(Mus);
currentMuscle = Mod.getMuscles().get(iMus-1);
if ~strcmp(fixname(char(currentMuscle.getName())),Mus)
    error('Current muscle name is incorrect')
end



% angles matrix: one position per row
[nrows,ncols] = size(angles);
nDofs = CoordSet.getSize;

if ncols~=nDofs
    if nrows~=nDofs
        errordlg('Angle matrix not the right size','Input Error');
        minusdLdq = [];
        return;
    else
        angles = angles';
    end
end

% initialise matrices for output
lengths = zeros(size(angles,1),1);
minusdLdq = zeros(size(angles,1),num_request_dofs);
if nargout>2
    GHfvecs = zeros(size(angles,1),3); 
end

for istep = 1:size(angles,1)
    if ~mod(istep,50)
        disp(['Muscle ',char(currentMuscle.getName()), ' - step ',...
            num2str(istep),' of ', num2str(size(angles,1))]);
    end
    
    % set dof values for this step
    for idof = 1:nDofs
        currentDof = CoordSet.get(idof-1);    
        currentDof.setValue(state,angles(istep,idof),1);
    end

    % get GH force line of action (if the muscle crosses GH)
    if nargout>2
        muspath = currentMuscle.getGeometryPath();
        fdarray = ArrayPointForceDirection();
        fvec2 = Vec3();
        muspath.getPointForceDirections(state,fdarray);
        scap_pt_index = -1;
        % find the "effective" muscle attachment on the scapula
        for ipt=1:fdarray.getSize-1
            body1 = char(fdarray.get(ipt-1).frame); %4.0 uses frames not bodies
            body2 = char(fdarray.get(ipt).frame);
            if strcmp(body1,'scapula_r')&&strcmp(body2,'humerus_r')
                scap_pt_index=ipt;
                break;
            elseif strcmp(body2,'scapula_r')&&strcmp(body1,'humerus_r')
                scap_pt_index=ipt-1;
                break;
            end
        end
        
        if scap_pt_index==-1
            % find the "effective" muscle attachment on the clavicle
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'clavicle_r')&&strcmp(body2,'humerus_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'clavicle_r')&&strcmp(body1,'humerus_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
        end

         if scap_pt_index==-1
            % find the "effective" muscle attachment on the thorax
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'thorax')&&strcmp(body2,'humerus_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'thorax')&&strcmp(body1,'humerus_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
         end

         if scap_pt_index==-1
            % find the "effective" muscle attachment on the ulna
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'scapula_r')&&strcmp(body2,'ulna_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'scapula_r')&&strcmp(body1,'ulna_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
         end
         
         if scap_pt_index==-1
            % find the "effective" muscle attachment on the radius
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'scapula_r')&&strcmp(body2,'radius_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'scapula_r')&&strcmp(body1,'radius_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
        end

        % calculate muscle force direction at that point (in global frame)
        scap_pt = fdarray.get(scap_pt_index);
        fvec = scap_pt.direction();
        
        % transform to the scapular coordinate frame
        SimEn.transform(state,groundbody,fvec,scapulabody,fvec2);
        GHfvecs(istep,1)=fvec2.get(0);
        GHfvecs(istep,2)=fvec2.get(1);
        GHfvecs(istep,3)=fvec2.get(2);
    end    
    
    % get length
    lengths(istep) = currentMuscle.getLength(state);

    % get moment arms
    for idof = 1:num_request_dofs
        dofindex  = Dofs(idof)-1;
        currentDof = CoordSet.get(Dofs(idof)-1);    
        
        % change the value of the dof by +-0.0001 and find length:
        currentDof.setValue(state,angles(istep,Dofs(idof))-0.0001,1);
        L1 = currentMuscle.getLength(state);
                        
        currentDof.setValue(state,angles(istep,Dofs(idof))+0.0001,1);                
        L2 = currentMuscle.getLength(state);
        
        minusdLdq(istep,idof) = -(L2-L1)/0.0002;
        
        % set dof to original value
        currentDof.setValue(state,angles(istep,Dofs(idof)),1);                
    end
end


%=============================================================================================
function musclepath_poly(muscles,mydir,musclepolyfile)

% Based on pathpoly.m by Ton van den Bogert
%
% Computes best fitting polynomial for muscle-tendon length L as a function
% of kinematic degrees of freedom q.
%
% This is done from moment arms generated by Opensim, using the relationship
% <momentarm of muscle m with respect to dof q> = -dLm/dq
%
% Stopping criteria: either error < 10% of maximum moment arm or 2mm, 
% whichever is greater, or change in error < 5% of error 
%
% input: structure "muscles", created by "das3_readosim.m"
%
% The moment arms and lengths are read from .mat files named: 
% "path_<muscle name>.mat", e.g. "path_bic_l.mat"
% 
% The polynomials are written in files: "<musclepolyfile>.mat" and
% "<musclepolyfile>.txt"
%
% Dimitra Blana, December 2011
%
% UPDATE June 18 2012: both lengths and moment arms are fitted (not just
% moment arms), and at least one term that includes every degree of freedom
% the muscle crosses is required.
%
% UPDATE August 2nd 2012: instead of including at least one term for each
% dof, the equations for each dof are *scaled* to the maximum moment arm/length 
% Also, the constant term of the polynomial is estimated as well.

%% open files and read in moment arms
format long;			% so we get full precision output for polynomial coefficients

% log file (output)
logfilename = [mydir '\pathpoly.log'];
logfile = fopen(logfilename,'w');
if (logfile < 0)
    errordlg(['Could not open log file ',logfilename,' for writing'],'File Error');
    fclose('all');
    return;
end

% results file (output)
polyfilename = [mydir '\' musclepolyfile '.txt'];
polyfile = fopen(polyfilename,'w');
if (polyfile < 0)
    errordlg(['Could not open results file ',polyfilename,' for writing'],'File Error');
    fclose('all');
    return;
end

% mat file, if it exists, append, don't overwrite
matfilename = [mydir '\' musclepolyfile '.mat'];
if exist(matfilename,'file')
    polys = load(matfilename);
    mus_model = polys.mus_model;
end

%% main loop for each muscle element
for imus = 1:length(muscles)
    
    % get moment arms and lengths out of the .mat files
    musfilename = [mydir,'\path_',muscles{imus}.name,'.mat'];
    ma = load(musfilename);
        
    mus = muscles{1,imus};
    ndofs = length(mus.dof_indeces); % number of dofs spanned by this muscle
    order = 4; % polynomial order
        
    num_data = size(ma.alljnts,1);
    
    % count how many parameters the polynomial model for muscle length has
    npar = prod(1:ndofs+order)/prod(1:ndofs)/prod(1:order);
    fprintf(1,'Muscle name:      %s\n',mus.name);
    fprintf(1,'Number of DOFs:   %d\n',ndofs);
    fprintf(1,'Polynomial order: %d\n',order);
    fprintf(1,'Potential number of polynomial terms: %d\n',npar);
    tot_data = num_data*(ndofs+1);	% total number of data points
    A = zeros(tot_data, npar);      % allocate memory space for A
    b = zeros(tot_data, 1);         % allocate memory space for b

    % get angle values for the dofs the muscle crosses
    musdof_indeces = zeros(ndofs,1);
    for idof = 1:ndofs
        imusdof = mus.dof_indeces(idof);
        musdof_indeces(idof) = imusdof;
    end
    ang = (ma.alljnts(:,musdof_indeces) + 1e-6);	% protect against angle = 0.0
   
    maxmomdof = zeros(1,ndofs);
    for idof = 1:ndofs
        maxmomdof(idof) = max(abs(ma.allmomarms(:,idof)))*1000;
    end    
    maxall = max(maxmomdof);
    % this normalises all moment arms
    ml_weight = (maxmomdof)/maxall;
    % Stopping criterion: error less than 10% of maximum moment arm (in mm) 
    % for the muscle or 2mm, whichever is greater
    momarm_error = max(0.1*maxall,2); 
    %momarm_error = 0.01*maxall; 

    for idof = 1:ndofs
        % read moment arm from allmomarms matrix
        % and angles from alljnts matrix
        b((idof-1)*num_data+1:idof*num_data) = -ma.allmomarms(:,idof)/ml_weight(idof);
        
        % generate the npar polynomial terms, and for each term, add a column to A
        polylist = zeros(npar,ndofs);
        expon = zeros(num_data,ndofs);	% start with all exponents zero
        for ii=1:npar
            polylist(ii,:) = expon(1,:);
            A((idof-1)*num_data+1:idof*num_data,ii) = (expon(:,idof).*prod(ang.^expon,2)./ang(:,idof))/ml_weight(idof); % contribution of this term to moment arm idof
            % generate the next set of exponents, but remain within model order
            k = 1;
            while (1)
                expon(:,k) = expon(:,k)+1;
                if (sum(expon(1,:)) > order && ii<npar)
                    expon(:,k)=0;
                    k = k+1;
                else
                    break;
                end
            end
        end     % done generating model terms
    end		% done reading all data for this muscle
    
    % <num_data> more rows for muscle length
    % read length from alllengths vector
    % and angles from alljnts matrix
    b(ndofs*num_data+1:(ndofs+1)*num_data) = ma.alllengths;

    % generate the npar polynomial terms, and for each term, add a column to A    
    for ipar=1:npar
        help = repmat(polylist(ipar,:),size(ang,1),1);
        A(ndofs*num_data+1:(ndofs+1)*num_data,ipar) = (prod(ang.^help,2));
    end  
        
    fprintf('Total number of data points: %d\n',tot_data);
   
    % now we have all data for this muscle stored in A and b
    
    % solve the full model with npar-1 terms
    p = A\b;		% compute coefficients of the best fitting model
    bpred = A*p;	% these are the moment arms predicted by this model
    res = bpred-b;	% residuals
    RMSfull = (sqrt(sum(res.^2)/tot_data)) * 1000;		% RMS of residuals, in mm
        
    fprintf('RMS fit error of the full model is: %f mm\n',RMSfull);
    fprintf('maximum moment arm: %f mm\n',maxall);
    fprintf('maximum error allowed: %f mm\n',momarm_error);
    fprintf(logfile,'%s\n',mus.name);
    fprintf(logfile,'  RMS fit error of the full model is: %f mm\n',RMSfull);
    fprintf(logfile,'  maximum moment arm: %f mm\n',maxall);
    fprintf(logfile,'  maximum error allowed: %f mm\n',momarm_error);

    % now do stepwise regression to select polynomial terms for a smaller model
    Aselected = [];
    polylist_selected = [];
    npar_selected = 0;
    % outer loop: successively add columns to Aselected
    for i = 1:npar-1
        [~, ncolumns] = size(A);
        % inner loop: find the column of A that causes most reduction in RMS when added to Aselected
        RMSnew = zeros(1,ncolumns); % this will store the RMS errors of each expanded model
        for j = 1:ncolumns
            % add column j from A to Anew
            Anew = [Aselected A(:,j)];
            % solve new p's
            pnew = Anew\b;
            % compute new RMS fit error
            RMSnew(j) = (norm(Anew*pnew - b)/(sqrt(tot_data))) * 1000; 	% convert to mm
        end
        % now determine which expanded model had the lowest RMS
        [RMSmin, col] = min(RMSnew);
        % if the change in error is less than 5%, stop without adding this term
%        if ((i>1)&&((RMS - RMSmin)/RMS<0.01))
        if ((i>1)&&((RMS - RMSmin)/RMS<0.05))
            fprintf('Change in error: %3f. No more terms added.\n ',(RMS - RMSmin)/RMS);
            fprintf(logfile,'Change in error: %3f. No more terms added.\n ',(RMS - RMSmin)/RMS);
            break;
        end
        % otherwise add this column to Aselected
        Aselected = [Aselected A(:,col)];
        npar_selected = npar_selected + 1;
        p = Aselected\b;		% solve the coefficients (again)
        % compute some error measures
        SSE = sum((Aselected*p-b).^2); 			% summed squared errors (SSE)
        RMS = (sqrt(SSE/tot_data)) * 1000; 			% RMS error, is the same as what we had before
        GCV = SSE/((num_data - npar_selected)^2);    	% Generalized Cross Validation
        AIC = log(SSE) + 2*npar_selected/tot_data;   	% Akaike's Information Criterion
        % print what we just did, on screen and on output file
        fprintf('Model addition step %3i: Added term ',i);
        fprintf('%i ',polylist(col,:));
        fprintf('-- RMS=%6.2f, GCV=%8.3e, AIC = %6.2f\n',RMS,GCV,AIC);
        fprintf(logfile,'  Model addition step %3i: Added term ',i);
        fprintf(logfile,'%i ',polylist(col,:));
        fprintf(logfile,'-- RMS=%f, GCV=%e, AIC = %f\n',RMS,GCV,AIC);

        % remember the exponents of this polynomial term
        polylist_selected = [polylist_selected ; polylist(col,:)];
        % remove this column from A and from polylist so it is not used again
        A = [A(:,1:(col-1)) A(:,(col+1):ncolumns)];
        polylist = [polylist(1:(col-1),:);polylist((col+1):ncolumns,:)];
        % stop adding terms if RMS error in fvectors is less than
        % momarm_error
        if ((RMS<=momarm_error))
            break;
        end
    end 		% and go find the next term

    % save muscle model in structure mus_model
    mus_model{imus}.name = mus.name;
    mus_model{imus}.num_lparams = npar_selected;
    mus_model{imus}.lparams = polylist_selected; % size of lparams: num_lparams x num_dofs
    mus_model{imus}.lcoef = p;
    
    % save muscle models in .mat file
    save(matfilename,'mus_model');
    
    % write this muscle's model in the output text file
    fprintf(polyfile,mus.name);
    fprintf(polyfile,'\nparameters %d\n',npar_selected);
    fprintf(polyfile,'# exponents ... coefficient\n');
    for i = 1:npar_selected
        fprintf(polyfile,'\t');
        fprintf(polyfile,'%3d ',polylist_selected(i,:));
        fprintf(polyfile,'   %10.5e \n', p(i));
    end
    fprintf(logfile,'  %d polynomial terms were written for %s\n',npar_selected+1, mus.name);

    % plot muscle length from Opensim and polynomial
    if ~mod(imus,20)
        examine_momarms(mus_model{imus}, mus.dof_names, ma.allmomarms, ang);	
    end
    
    clear ma A b ang;
    
end 		% go process next muscle

fclose(polyfile);
fclose(logfile);


%=============================================================================================
function examine_momarms(musmodel, dof_names, moment_arms, angles)
% plot momentarm-angle data

% choose a subset of "angles" that contains only 100 points
a = 1;
b = size(angles,1);
r = a + (b-a).*rand(100,1);
indeces = ceil(sort(r));
sangles = angles(indeces,:);

% ...or use all angles
%indeces = 1:size(angles,1);
%sangles = angles;

% calculate moment arms from polynomial
pmoment_arms = zeros(length(indeces),length(dof_names));
for iframe = 1:length(indeces)
    for i=1:musmodel.num_lparams

        % add this term's contribution to the muscle length 
        term = musmodel.lcoef(i);

        for j=1:length(dof_names)
            for k=1:musmodel.lparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end

        % first derivatives of length with respect to all q's
        for  k=1:length(dof_names)
            % derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
            if ((musmodel.lparams(i,k) > 0) && (sangles(iframe,k)))	
                dterm = musmodel.lparams(i,k)*term/sangles(iframe,k);
                pmoment_arms(iframe,k) = pmoment_arms(iframe,k) + dterm;
            end
        end
    end
end

figure;
for idof=1:length(dof_names)
    subplot(length(dof_names),1,idof);
    plot(moment_arms(indeces,idof),'bx-'); hold on; plot(-pmoment_arms(:,idof),'ro-'); 	
    title([dof_names{idof}, ' momentarms for ',musmodel.name],'Interpreter', 'none'); 
end
legend('osim','poly');


%=============================================================================================
function make_mot_file(motfilename,angles,dofnames)
% builds .mot file for osim using angles
%
% Dimitra Blana, November 2011

motfile = fopen(motfilename, 'wt');
if motfile==-1
    errordlg('Could not open .mot file for writing','File Error');
    return;
end

dofnames=strjoin(['time '; dofnames]');

[num_rows,num_columns] = size(angles);
time = 0:1/(num_rows-1):1;

fprintf(motfile, 'name %s\n',motfilename);
fprintf(motfile, 'datacolumns %d\n',num_columns+1);
fprintf(motfile, 'datarows %d\n',num_rows);
fprintf(motfile, 'range 0 1\n');
fprintf(motfile, 'endheader\n');
fprintf(motfile, '%s\n',dofnames);
for irow = 1:num_rows
    fprintf(motfile, '%f\t',time(irow));
    for idof = 1:num_columns
        fprintf(motfile, '%f\t',angles(irow,idof)*180/pi);
    end
    fprintf(motfile, '\n');
end

fclose(motfile);

%=============================================================================================
function GH_poly(muscles,mydir,musclepolyfile)

% Based on pathpoly.m by Ton van den Bogert
%
% Computes best fitting polynomial for GH force vector as a function
% of kinematic degrees of freedom q.
%
% This is done from force vectors generated by Opensim,
%
% Stopping criteria: either error <10% of range, or change in error <5% of error 
%
% input: structure "muscles", created by "das3_readosim.m"
%
% The force vectors are read from .mat files named: "GHfvec_<muscle name>.mat", 
% e.g. "GHfvec_bic_l.mat"
% 
% The polynomials are written in files: "<musclepolyfile>.mat" and
% "<musclepolyfile>.txt"
%
% Dimitra Blana, October 2014

%% open files and read in force vectors
format long;			% so we get full precision output for polynomial coefficients

% log file (output)
logfilename = [mydir '\GHpoly.log'];
logfile = fopen(logfilename,'w');
if (logfile < 0)
    errordlg(['Could not open log file ',logfilename,' for writing'],'File Error');
    fclose('all');
    return;
end

% results file (output)
polyfilename = [mydir '\' musclepolyfile '.txt'];
polyfile = fopen(polyfilename,'w');
if (polyfile < 0)
    errordlg(['Could not open results file ',polyfilename,' for writing'],'File Error');
    fclose('all');
    return;
end

% mat file, if it exists, append, don't overwrite
matfilename = [mydir '\' musclepolyfile '.mat'];
% if exist(matfilename,'file')
%     polys = load(matfilename);
%     GH_model = polys.GH_model;
% end

forcedir = ['x','y','z'];
GH_model = cell(1,length(muscles));

%% main loop for each muscle element
for imus = 1:length(muscles)

    % check it the muscle crosses GH
    crossGH=0;
    for idof = 1:muscles{imus}.dof_count
        if strcmp(muscles{imus}.dof_names{idof},'GHy')
            crossGH=1; break;
        end
    end
    if ~crossGH, continue; end
    
    % get force vectors out of the .mat files
    musfilename = [mydir,'\GHfvec_',muscles{imus}.name,'.mat'];
    ma = load(musfilename);
        
    mus = muscles{1,imus};
    ndofs = length(mus.dof_indeces); % number of dofs spanned by this muscle
    order = 4; % polynomial order
        
    num_data = size(ma.alljnts,1);
    
    % count how many parameters the polynomial model for GH force vector has
    npar = prod(1:ndofs+order)/prod(1:ndofs)/prod(1:order);
    fprintf(1,'Muscle name:      %s\n',mus.name);
    fprintf(1,'Number of DOFs:   %d\n',ndofs);
    fprintf(1,'Polynomial order: %d\n',order);
    fprintf(1,'Potential number of polynomial terms: %d\n',npar);

    % get angle values for the dofs the muscle crosses
    musdof_indeces = zeros(ndofs,1);
    for idof = 1:ndofs
        imusdof = mus.dof_indeces(idof);
        musdof_indeces(idof) = imusdof;
    end
   
    % GH fvector data for this muscle come from ndofs files, with names
    % listed in the config file
    for i = 1:3  % three components of the fvector (x-y-z)
        A = zeros(num_data, npar);	% allocate memory space for A

        % read data from GHfvector matrix
        % take fvector and angles from the right columns
        b = ma.allGHfvecs(:,i);
        ang = (ma.alljnts(:,musdof_indeces) + 1e-6);	% protect against angle = 0.0

        % generate the npar polynomial terms, and for each term, add a column to A
        polylist = zeros(npar,ndofs);
        expon = zeros(num_data,ndofs);	% start with all exponents zero
        for ii=1:npar
            polylist(ii,:) = expon(1,:);
            A(:,ii) = prod(ang.^expon,2); % contribution of this term to force component i
            % generate the next set of exponents, but remain within model order
            k = 1;
            while (1)
                expon(:,k) = expon(:,k)+1;
                if (sum(expon(1,:)) > order && ii<npar)
                    expon(:,k)=0;
                    k = k+1;
                else
                    break;
                end
            end
        end     % done generating model terms
        fprintf('Total number of data points: %d\n',num_data);
        
        vecrange = abs(max(b)-min(b));
        
        % now we have all data for this muscle stored in A and b
        % solve the full model with npar terms
        p = A\b;		% compute coefficients of the best fitting model

        bpred = A*p;	% these are the force vectors predicted by this model
        res = bpred-b;	% residuals
        RMSfull = (sqrt(sum(res.^2)/num_data));		% RMS of (normalized) residuals
        fprintf('RMS fit error of the full model is: %f (normalized)\n',RMSfull);
        fprintf('min amplitude: %f\n',min(b));
        fprintf('max amplitude: %f\n',max(b));
        fprintf('range of amplitudes: %f\n',vecrange);
        fprintf(logfile,'%s\n',mus.name);
        fprintf(logfile,'  RMS fit error of the full model is: %f (normalized)\n',RMSfull);
        fprintf(logfile,'range of amplitudes: %f\n',vecrange);

        % now do stepwise regression to select polynomial terms for a smaller model
        Aselected = [];
        polylist_selected = [];
        npar_selected = 0;
        % outer loop: successively add columns to Aselected
        for ii = 1:npar-1
            ncolumns = size(A,2);
            % inner loop: find the column of A that causes most reduction in RMS when added to Aselected
            RMSnew = zeros(1,ncolumns); % this will store the RMS errors of each expanded model
            for j = 1:ncolumns
                % add column j from A to Anew
                Anew = [Aselected A(:,j)];
                % solve new p's
                pnew = Anew\b;
                % compute new RMS fit error
                RMSnew(j) = (norm(Anew*pnew - b)/(sqrt(num_data))); 
            end
            % now determine which expanded model had the lowest RMS
            [RMSmin, col] = min(RMSnew);
            % if the change in error is less than 5%, stop without adding this term 
            if ((ii>1)&&((RMS - RMSmin)/RMS<0.05))
                fprintf('Change in error: %3f. No more terms added.\n ',(RMS - RMSmin)/RMS);
                fprintf(logfile,'Change in error: %3f. No more terms added.\n ',(RMS - RMSmin)/RMS);
                break;
            end
            % otherwise add this column to Aselected
            Aselected = [Aselected A(:,col)];npar_selected = npar_selected + 1;
            p = Aselected\b;		% solve the coefficients (again)
            % compute some error measures
            SSE = sum((Aselected*p-b).^2); 			% summed squared errors (SSE)
            RMS = (sqrt(SSE/num_data)); 			% RMS error, is the same as what we had before
            
            GCV = SSE/((num_data - npar_selected)^2);    	% Generalized Cross Validation
            AIC = log(SSE) + 2*npar_selected/num_data;   	% Akaike's Information Criterion
            % print what we just did, on screen and on output file
            fprintf('Model addition step %3i: Added term ',ii);
            fprintf('%i ',polylist(col,:));
            fprintf('-- RMS=%6.2f, GCV=%8.3e, AIC = %6.2f\n',RMS,GCV,AIC);
            fprintf(logfile,'  Model addition step %3i: Added term ',ii);
            fprintf(logfile,'%i ',polylist(col,:));
            fprintf(logfile,'-- RMS=%f, GCV=%e, AIC = %f\n',RMS,GCV,AIC);

            % remember the exponents of this polynomial term
            polylist_selected = [polylist_selected ; polylist(col,:)];
            % remove this column from A and from polylist so it is not used again
            A = [A(:,1:(col-1)) A(:,(col+1):ncolumns)];
            polylist = [polylist(1:(col-1),:);polylist((col+1):ncolumns,:)];           
            % stop adding terms if RMS error in fvectors is less than 10% of range, 
            if (RMS<=0.1*vecrange)
                break;
            end
        end 		% and go find the next term

        % write this muscle's model on the output file
        fprintf(polyfile,'parameters %d\n',npar_selected);
        fprintf(polyfile,'# exponents ... coefficient for %c force\n',forcedir(i));
        for ii = 1:npar_selected
            fprintf(polyfile,'\t');
            fprintf(polyfile,'%3d ',polylist_selected(ii,:));
            fprintf(polyfile,'   %10.5e \n', p(ii));
        end
        fprintf(logfile,'  %d polynomial terms were written for %s %c force\n',npar_selected,mus.name,forcedir(i));
        fprintf('  %d polynomial terms were written for %s %c force\n',npar_selected,mus.name,forcedir(i));

        switch i
            case 1
                % save muscle model in structure mus_model
                GH_model{imus}.name = mus.name;
                GH_model{imus}.xparam_count = npar_selected;
                GH_model{imus}.xparams = polylist_selected; % size of lparams: num_lparams x num_dofs
                GH_model{imus}.xcoef = p;
            case 2
                % save muscle model in structure mus_model
                GH_model{imus}.yparam_count = npar_selected;
                GH_model{imus}.yparams = polylist_selected; % size of lparams: num_lparams x num_dofs
                GH_model{imus}.ycoef = p;
            case 3
                % save muscle model in structure mus_model
                GH_model{imus}.zparam_count = npar_selected;
                GH_model{imus}.zparams = polylist_selected; % size of lparams: num_lparams x num_dofs
                GH_model{imus}.zcoef = p;
        end
        clear A  B;
        
    end		% done reading this force direction
    
    % save muscle models in .mat file
    save(matfilename,'GH_model');   
    
    % plot GH fvector from Opensim and polynomial
    if ~mod(imus,20)
        examine_fvectors(GH_model{imus},ma.allGHfvecs,ang);	
    end
    
    clear ma ang;
    
end 		% go process next muscle

fclose(polyfile);
fclose(logfile);


%=============================================================================================
function examine_fvectors(musmodel,fvectors, angles)
% plot momentarm-angle data

% choose a subset of "angles" that contains only 100 points
a = 1;
b = size(angles,1);
r = a + (b-a).*rand(100,1);
indeces = ceil(sort(r));
sangles = angles(indeces,:);

% ...or use all angles
%indeces = 1:size(angles,1);
%sangles = angles;

% calculate fvectors from polynomial
pfvectors = zeros(length(indeces),3);
for iframe = 1:length(indeces)
    for i=1:musmodel.xparam_count
        % add this term's contribution to the x vector
        term = musmodel.xcoef(i);
        for j=1:size(musmodel.xparams,2)
            for k=1:musmodel.xparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end
        pfvectors(iframe,1) = pfvectors(iframe,1) + term;
    end
    for i=1:musmodel.yparam_count
        % add this term's contribution to the y vector
        term = musmodel.ycoef(i);
        for j=1:size(musmodel.yparams,2)
            for k=1:musmodel.yparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end
        pfvectors(iframe,2) = pfvectors(iframe,2) + term;
    end
    for i=1:musmodel.zparam_count
        % add this term's contribution to the z vector
        term = musmodel.zcoef(i);
        for j=1:size(musmodel.zparams,2)
            for k=1:musmodel.zparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end
        pfvectors(iframe,3) = pfvectors(iframe,3) + term;
    end
end

figure;
dir_names = {'x','y','z'};
for idir=1:3
    subplot(3,1,idir);
    plot(fvectors(indeces,idir),'bx-'); hold on; plot(pfvectors(:,idir),'ro-'); 	
    title([dir_names{idir}, ' GH force vectors for ',musmodel.name],'Interpreter', 'none'); 
end
legend('osim','poly');


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