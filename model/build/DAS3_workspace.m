function angles = DAS3_workspace(input1,input2,input3)
% find workspace for DAS3

ndof = 11;
nmus = 138;
nstates = 2*ndof + 2*nmus;
dofnames = { 'SC_y' 'SC_z' 'SC_x' 'AC_y' 'AC_z' 'AC_x' 'GH_y' 'GH_z', 'GH_yy' 'EL_x', 'PS_y'};

x = zeros(ndof,1);
x(10) = pi/18;
x(11) = pi/2;

angles=[];

if nargin==1
    das2angles = load(input1);
    allx = das2angles.GHdata;
    for index=1:size(allx,1)
        x(1:9) = estimate_shoulder_angles(allx(index,1)*pi/180,allx(index,2)*pi/180,allx(index,3)*pi/180);
        angles=[angles x];  
    end
    angles = angles';

elseif nargin==3
    GHy = input1;
    GHz = input2;
    GHyy = input3;
    for hum_thory=GHy
        for hum_thorz=GHz
            for hum_thoryy=GHyy
                x(1:9) = estimate_shoulder_angles(hum_thory*pi/180,hum_thorz*pi/180,hum_thoryy*pi/180);
                angles=[angles x];               
            end
        end
    end
    angles = angles';
end
