function das3_buildmodel
	% builds the das3 MEX function
    
	% load the Opensim model
    addpath('..');
	model = das3_readosim('../das3.osim');
    rmpath('..');
    
    % make the Autolev source file das3.al
    das3_makeAutolev(model, 'das3.al');  
	
    % Use Autolev to generate the kinematics/dynamics C code
    minutes = 4;
    runautolev('das3.al', minutes);

    % clean the C code generated by Autolev to make it a callable function
    autolevclean(model, 'das3_al_raw.c', 'das3_al.c');
	
	% make das3.h and das3.c
	das3_makeMEX(model, 'das3_mex_template.c', 'das3.c');
      
    % Compile the C code to make the MEX binary
    % On Windows, use the MinGW compiler (https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler)
    % Remove the -g after successful testing.  This turns on compiler
    % optimization, for a faster MEX function, but compilation is a lot
    % slower
    MEXOPTIONS = '-g -largeArrayDims';
    fprintf('Compiling the C code...\n');
    if findstr(MEXOPTIONS,'-g')
        hours = 2; 
        fprintf('Estimated completion at %s...\n', datestr(addtodate(now, hours, 'hour')));
    else
        minutes = 2; 
        fprintf('Estimated completion at %s...\n', datestr(addtodate(now, minutes, 'minute')));
    end
    mex(MEXOPTIONS,'das3.c','das3_al.c');  % creates das3.mexw64 on windows
    
    % create the .mat files with polynomials for muscle length and muscle
    % contributions to GH reaction force
    das3_polynomials('../das3.osim', './tmp', 'musclepoly.mat', 'GHpoly.mat')
	
end
%==============================================================================================
function runautolev(filename, minutes)
	% Autolev location depends on computer
	[~,computername] = system('hostname');
	if strfind(computername,'LRI-102855')     
		autolev = 'C:\Program Files\Autolev\al.exe';
	else
		autolev = [];
    end
    
    % name of files generated by autolev
    cfile = strrep(filename, '.al', '_al_raw.c');
    infile = strrep(filename, '.al', '_al_raw.in');
   		
	% Run Autolev
	if ~exist(autolev)
		disp(['Error: Autolev needs to generate ' filename '_raw.c, but you do not have Autolev installed.']);
		disp('Ask Ton to run this part of the build process.');
		disp(['Hit ENTER to continue with existing ' cfile ', or CTRL-C to quit']);
		pause
    else
        % we need to delete the .c and .in files so Autolev won't ask for overwrite permission
		warning('off', 'MATLAB:DELETE:FileNotFound');	% don't warn me if the file does not exist
		delete(cfile);		
		delete(infile);
		warning('on', 'MATLAB:DELETE:FileNotFound');
		
		fprintf('Autolev is generating %s\n', cfile);
		fprintf('Estimated completion at %s...\n', datestr(addtodate(now, minutes, 'minute')));
		
		if ispc
			system(['"' autolev '" ' filename ' > nul']);		% double quotes needed because autolev has spaces in its path
		else
			system(['wine "' autolev '" ' filename ' > nul']);	% Linux/Mac can use wine (assuming autolev installed in the same place)
		end
		% check if Autolev was successful
		if ~exist(cfile)
			warning(['Autolev error in ' filename]);
		end
		delete(infile);
		fprintf('Done.\n');
	end
end
%===================================================================================================
function autolevclean(model, rawcfile, cfile)
	disp('Cleaning C code generated by Autolev...');

	% open the raw C file that came from Autolev
	fid1 = fopen(rawcfile,'r');
	if (fid1 == -1)
		error('Could not open %s\n', rawcfile);
	end
	
	% write the clean C file
	fid2 = fopen(cfile,'w');
	if (fid2 == -1)
		error('Could not write %s\n', cfile);
	end
	
	% write the function header
	fprintf(fid2,'// This file contains C code generated by Autolev\n\n');
	fprintf(fid2,'#include "das3.h"\n');
	fprintf(fid2,'#include <math.h>\n\n');
	fprintf(fid2,'void das3_al(param_struct* par, double q[NDOF], double qd[NDOF], double qdd[NDOF],\n');
	fprintf(fid2,'   double mTH[3], double exF[2], double handF[3], double Zero[NDOF], double dz_dq[NDOF][NDOF],\n');
	fprintf(fid2,'   double dz_dqd[NDOF][NDOF], double dz_dqdd[NDOF][NDOF],\n');
	fprintf(fid2,'   double F_GH[3], double F_SCAP[2][3], double dist_SCAP[2], double pos_SCAP[2][3],\n');
	fprintf(fid2,'   double vis[NVIS][12], double qTH[3]) {\n');

	% make a macro for the 'sign' function
	fprintf(fid2,'#define sign(x) ((x > 0) - (x < 0))\n');

	% generate C code to copy q, qd, and qdd into scalar variables
	for i=1:model.nDofs
		fprintf(fid2,'	double %s = q[%1d];\n', model.dofs{i}.name,i-1);
		fprintf(fid2,'	double %sp = qd[%1d];\n', model.dofs{i}.name,i-1);
		fprintf(fid2,'	double %spp = qdd[%1d];\n', model.dofs{i}.name,i-1);
    end		

	% declare some internal variables
	fprintf(fid2,'	double MHx, MHy, MHz;\n');

	% generate C code to copy the 3 thorax-humerus into scalar variables
	fprintf(fid2,'	double MTHy  = mTH[0];\n');
	fprintf(fid2,'	double MTHz  = mTH[1];\n');
	fprintf(fid2,'	double MTHyy = mTH[2];\n');

	% generate C code to copy forearm force description (arm support) into scalar variables
	fprintf(fid2,'	double distF = exF[0];\n');
	fprintf(fid2,'	double ampF  = exF[1];\n');

	% generate C code to copy hand force description into scalar variables
	fprintf(fid2,'	double handFx = handF[0];\n');
	fprintf(fid2,'	double handFy = handF[1];\n');
	fprintf(fid2,'	double handFz = handF[2];\n');

	% generate C code to declare the temporary variables used in contact force model and conoid force model
	fprintf(fid2,'	double PxCE,PyCE,PzCE;\n');		
	fprintf(fid2,'	double PxTS,PyTS,PzTS,FxTS,FyTS,FzTS,FTS,FminusTS;\n');		
	fprintf(fid2,'	double PxAI,PyAI,PzAI,FxAI,FyAI,FzAI,FAI,FminusAI;\n');	
	fprintf(fid2,'	double LOI, Stretch, StretchPositive;\n');		

	% copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while ~feof(fid1)
		line = fgetl(fid1);
		if strncmp(line, 'double   Pi,DEGtoRAD,RADtoDEG,z[', 32)
			zlength = line(33:min(findstr(line,']'))-1);
			fprintf(fid2,'\tstatic double z[%d];\n',str2num(zlength));		% make sure there is enough room for all Zs
		end
		if strcmp(line, '/* Evaluate constants */') 				% Z[] code starts here
			copying = 1;
        elseif strcmp(line, '/* Evaluate output quantities */') 	% Z[] code ends here
			copying = 0;
        elseif strcmp(line, '/* Write output to screen and to output file(s) */')   % encoded variables code starts here
			copying = 1;
        elseif strcmp(line, '  Encode[0] = 0.0;') 								   % and stops here
			copying = 0;
        elseif copying
			line = strrep(line, 'par__', 'par->');			% change par__ into par->
			fprintf(fid2,'%s\n',line);
		end
	end
	
	% close the input file
	fclose(fid1);
	
	% close the output file
	fprintf(fid2,'}\n');
	fclose(fid2);

end
%===================================================================================================
