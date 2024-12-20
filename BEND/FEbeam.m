% Basis beam finite element program
% Computes and plots input receptance vs. angular frequency
% Should be modified according to exercise problems
% Created 21/9-2014 by JSJ
function FEbeam(ne,eta)

% input parameters
% ne                    % Number of elements
% eta					% Loss factor

% Material parameters (to be inserted)
L = ?;				% Beam length
mp = ?;				% Mass per length
B = ?;				% Bending stiffness

% Element length
le = L/ne;

% Number of nodes and degrees of freedom computed
nn = ne+1;				% Number of nodes
neqn = 2*nn;			% Number of equations
disp(['Number of DOF ' sprintf('%d',neqn) ...
     ' Number of elements ' sprintf('%d',ne)]);

 % Global variables initialised
S = sparse(neqn,neqn);                    	% Stiffness matrix
M = sparse(neqn,neqn);                    	% Mass matrix
F = sparse(neqn,1);                          % Force vector

 % Build force vector (location to be inserted via correct DOF)
fdof = ?;                              % DOF for force
F(fdof,1) = 1;                              % unit force

% Element values of mp and B (change if inhomogeneous)
rhoA = mp*ones(ne,1);
EI = B*(1+1i*eta)*ones(ne,1);

% Matrices assembled
for e = 1:ne

    % global dofnumbers for local element
    edof = [2*e-1 2*e 2*(e+1)-1 2*(e+1)];

    % create local stiffness matrix
    se = (2*EI(e,1)/le^3)*[6    3*le   -6    3*le
                           3*le 2*le^2 -3*le le^2
                           -6   -3*le  6     -3*le
                           3*le le^2   -3*le 2*le^2];
                       
    % put into global stiffness matrix
    S(edof,edof) = S(edof,edof) + se;

    % create local mass matrix
    me = (rhoA(e,1)*le/840)*[312    44*le   108    -26*le
                             44*le  8*le^2  26*le  -6*le^2
                             108    26*le   312    -44*le
                             -26*le -6*le^2 -44*le 8*le^2];
    
    % put into global stiffness matrix
    M(edof,edof) = M(edof,edof) + me;

end

% loop over angular frequencies
nom = 5000;					% number of ang. frq's
om1 = 1;                    % start ang. frq
om2 = 50000;				% end ang. frq
for iom = 1:nom
	om(iom,1) = om1 + (iom-1)*(om2-om1)/(nom-1);
	X = (S-om(iom,1)^2*M)\F;
	% save input receptance
	Xinp(iom,1) = X(fdof,1);
end

% plot results
figure(1);
plot(om,log(abs(Xinp)));