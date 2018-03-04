function test_wkm (test)
%TEST_WKM Test the output of the MATLAB function WKM which computes the 
%   wake-kernel matrix functions COSH(SQRT(A)) and SINHC(SQRT(A)).
%   The Symbolic Math Toolbox (VPA) is assumed to be present for the tests.
%   Tested using MATLAB R2017b.
%   See also WKM, VPA
%
%   USAGE: test_wkm (test);
%     where, test = {'mct_testmats', 'expm_testmats', 'logm_testmats'}
%   'mct_testmats' : suite of 52 matrices (in matrix.m) 
%                    from Higham's <a href="matlab: web('http://www.maths.manchester.ac.uk/~higham/mctoolbox')">Matrix Computation Toolbox</a>
%   'expm_testmats': suite of 32 matrices; 
%                    Al-Mohy and Higham used them to test expm
%   'logm_testmats': suite of  8 matrices; 
%                    Higham used them to test logm
%   Note: These matrices are used as is. No scaling is done to avoid overflows
%
%   The output is organized as a 2D array. Each row corrresponds to the output
%   of its corresponding matrix. The columns are organized as follows.
%
%   id, tWKM, tFUNM, tVPA, errwkmcA, errwkmsA, errfunmcA, errfunmsA, cestfcA, cestfsA
%
%   where the acronyms in the previous sentence mean
%     id        : matrix identifier text
%     tWKM      : time for the Nadukandi--Higham algorithm (WKM)
%     tFUNM     : time for the Davies--Higham algorithm (FUNM)
%     tVPA      : time to compute wave kernels using the Davis trick in vpa
%     errWKMcA  : relative forward error 1st wave kernel using WKM 
%     errWKMsA  : relative forward error 2nd wave kernel using WKM 
%     errFUNMcA : relative forward error 1st wave kernel using FUNM 
%     errFUNMsA : relative forward error 2nd wave kernel using FUNM 
%     cestfcA   : 1st wave kernel condition number estimate cond(fc,A)*eps 
%     cestfsA   : 2nd wave kernel condition number estimate cond(fs,A)*eps 
%
%   Note: cestfcA and cesrfsA are computed in VPA. So this part will be slow.
%   The function FUNM_CONDEST1 (source code included in TEST_WKM) is used to 
%   compute the condition number estimates of the matrix functions and it is a 
%   part of Higham's 
%   <a href="matlab: web('http://www.maths.manchester.ac.uk/~higham/mftoolbox')">Matrix Function Toolbox</a>
%
%   TODO: use WKM for condition number estimation. write a wrapper to WKM that
%         whose output is only the 2nd wave kernel.
%
%   References:
%   P. Nadukandi and N. J. Higham, Computing the Wave-Kernel Matrix Functions,
%      Manchester Institute for the Mathematical Sciences, <a href="matlab: web('http://eprints.maths.manchester.ac.uk/2621/')">MIMS EPrint 2018.4</a>, 
%      2018.
%
%   Funding:
%   Nadukandi was supported by an individual fellowship from the 
%     European Union's Horizon 2020 research and innovation programme under the 
%     Marie Sklodowska-Curie grant <a href="matlab: web('https://cordis.europa.eu/project/rcn/200435_en.html')">702138</a>.
%   Higham was supported by Engineering and Physical Sciences Research Council 
%     grant <a href="matlab: web('http://gow.epsrc.ac.uk/NGBOViewGrant.aspx?GrantRef=EP/P020720/1')">EP/P020720/1</a>.

  d_vpa = 250; % set digits to use in variable precision arithmetic (vpa)
  hf_vpa =@(s,d) vpa(sym(s,'f'),d); %copied from Nick's myvpa.m; cf. help sym
  % IMP : use hf_vpa if you want vpa on the hardware float representations
  % E.g. vpa(0.1, 50) = 0.1, but
  %   hf_vpa(0.1, 50) = 0.10000000000000000555111512312578270211815834045410
  % as sym(0.1,'f') = '1.999999999999a'*2^(-4) (mantissa: 4-bit -> hexadecimal)

  % Use a multiplier to ensure that some scaling will occur for the test
  % matrices in 'mct_testmats', 'expm_testmats' and 'logm_testmats'.
  multiplier = false;
  % multiplier = true;
  multiplier_size = 60;

  switch lower(test) % IMP: MATLAB switch does not "fall through"
    case 'mct_testmats'
      n = 15;
      randmult = rand(1,matrix(0))*multiplier_size;
      for k = 1:matrix(0)
        A = matrix(k,n);
        if multiplier, A = multiplier_size*A; end % to enforce constant scaling
        % if multiplier, A = randmult(k)*A; end % to enforce random scaling
        compute_and_display(['matrix(', num2str(k, '% 2d'),', 15)'], A)
      end

    case 'expm_testmats'
      n = 15;
      [A,numexp] = expm_testmats(1);
      randmult = rand(1,numexp)*multiplier_size;
      for k = 1:numexp
        A = expm_testmats(k,n);
        if multiplier, A = multiplier_size*A; end % to enforce constant scaling
        % if multiplier, A = randmult(k)*A; end % to enforce random scaling
        compute_and_display(['expm_testmats(', num2str(k, '% 2d'),', 15)'], A)
      end

    case 'logm_testmats'
      n = 15;
      [A,numlog] = logm_testmats(1);
      randmult = rand(1,numlog)*multiplier_size;
      for k = 1:numlog
        A = logm_testmats(k,n);
        if multiplier, A = multiplier_size*A; end % to enforce constant scaling
        % if multiplier, A = randmult(k)*A; end % to enforce random scaling
        compute_and_display(['logm_testmats(', num2str(k, '% 2d'),', 15)'], A)
      end

    otherwise
      error('Unknown test.');
  end

end

function fcA = fcA_vpa (A)
  fc =@(x) cosh(sqrt(x));

  d_vpa = 250; % set digits to use in variable precision arithmetic (vpa)
  hf_vpa =@(s,d) vpa(sym(s,'f'),d); %copied from Nick's myvpa.m; cf. help sym
  d_old = digits(); digits(d_vpa);
  
  % E.B.Davies trick: ensures eigenvalues are distinct; diagonalisation possible
  del_A = randn(length(A)); % random perturbation
  del_A = 10^(-d_vpa/2)*del_A/norm(del_A,1); % make it of norm half of vpa
  [V, D] = eig(hf_vpa(A, d_vpa) + del_A); % E.B.Davies's trick
  
  fcA = double(V*diag(fc(diag(D)))/V);
  digits(d_old);
end

function fsA = fsA_vpa (A)
  sinhc =@(x) sinh(x)./x;
  fs =@(x) sinhc(sqrt(x));

  d_vpa = 250; % set digits to use in variable precision arithmetic (vpa)
  hf_vpa =@(s,d) vpa(sym(s,'f'),d); %copied from Nick's myvpa.m; cf. help sym
  d_old = digits(); digits(d_vpa);
  
  % E.B.Davies trick: ensures eigenvalues are distinct; diagonalisation possible
  del_A = randn(length(A)); % random perturbation
  del_A = 10^(-d_vpa/2)*del_A/norm(del_A,1); % make it of norm half of vpa 
  [V, D] = eig(hf_vpa(A, d_vpa) + del_A); % E.B.Davies's trick
  
  fsA = double(V*diag(fs(diag(D)))/V);
  digits(d_old);
end

function compute_and_display (tag, A)

  % Compute the fc(A) and fs(A) in vpa to compute approximation errors.
  tic;
  fcA = fcA_vpa(A);
  fsA = fsA_vpa(A);
  t_vpa = toc;
  
  % Compute the fc(A) and fs(A) using Nadukandi--Higham algorithm.
  tic;
  [rcA, rsA] = wkm(A);
  t_apx = toc;
  err_cA = norm(rcA-fcA,1)/norm(fcA,1);
  err_sA = norm(rsA-fsA,1)/norm(fsA,1);
  
  % Compute the fc(A) and fs(A) using funm.
  % get the k'th derivative of cosh(sqrt) and sinhc(sqrt) at x
  % to be used as argument to the Davies--Higham MATLAB function FUNM
  % syntax: funm(A, dfc_k); to compute cosh( sqrt(A)) using funm
  % syntax: funm(A, dfs_k); to compute sinhc(sqrt(A)) using funm

  % Functions to compute k'th derivative of cosh(sqrt) and sinhc(sqrt) at x 
  % Reference: see page 3 of MIMS EPrint 2018.4
  dfc_k =@(x,k) hypergeom([],[(1/2)+k],(x/4))/prod(k  +(1:k));
  dfs_k =@(x,k) hypergeom([],[(3/2)+k],(x/4))/prod(k+1+(0:k));

  tic;
  % funm_cA = funm(A, @fun_coshsqrt); % NJH's suggestion. Symbolic derivative
  funm_cA = funm(A, dfc_k);
  funm_sA = funm(A, dfs_k);
  t_funm = toc;
  err_funm_cA = norm(funm_cA-fcA,1)/norm(fcA,1);
  err_funm_sA = norm(funm_sA-fsA,1)/norm(fsA,1);

  % Estimate the conditioning of fc(A) and fs(A)
  try
    cnd_fcA = funm_condest1(A,@fcA_vpa)*eps;
    cnd_fsA = funm_condest1(A,@fsA_vpa)*eps;
    cnd_fcA_str = num2str(cnd_fcA, '%1.2e');
    cnd_fsA_str = num2str(cnd_fsA, '%1.2e');
  catch % error: NaN or Inf found.
    cnd_fcA_str = 'NaN or Inf';
    cnd_fsA_str = 'NaN or Inf';
  end

  % Display results
  disp([tag, ', ',...
        't_apx  = ',num2str(t_apx,  '%1.2e'),', '...
        't_funm = ',num2str(t_funm, '%1.2e'),', ',...
        't_vpa  = ',num2str(t_vpa,  '%1.2e'),', ',...
        'err_cA = ',num2str(err_cA, '%1.2e'),', ',...
        'err_sA = ',num2str(err_sA, '%1.2e'),', ',...
        'err_funm_cA = ',num2str(err_funm_cA, '%1.2e'),', ',...
        'err_funm_sA = ',num2str(err_funm_sA, '%1.2e'),', ',...
        'cond(fc,A)*eps = ',cnd_fcA_str, ', ',...
        'cond(fs,A)*eps = ',cnd_fsA_str]);

end

%=============================================================================
% Below are some matlab script files taken from Nick's MCToolBox and MFToolBox
% matrix.m, expm_testmats.m, logm_testmats.m and fun_condest1.m
%-----------------------------------------------------------------------------

function A = matrix(k, n)
%MATRIX  Test matrices accessed by number.
%        MATRIX(K, N) is the N-by-N instance of matrix number K in
%        a set of test matrices comprising those in MATLAB plus those
%        in the Matrix Computation Toolbox,
%        with all other parameters set to their default.
%        N.B. - Only those matrices which are full and take an arbitrary
%               dimension N are included.
%             - Some of these matrices are random.
%        MATRIX(K) is a string containing the name of the K'th matrix.
%        MATRIX(0) is the number of matrices, i.e. the upper limit for K.
%        Thus to set A to each N-by-N test matrix in turn use a loop like
%             for k=1:matrix(0)
%                 A = matrix(k, N);
%                 Aname = matrix(k);   % The name of the matrix
%             end
%        MATRIX(-1) returns the version number and date of the
%        Matrix Computation Toolbox.
%        MATRIX with no arguments lists the names and numbers of the M-files in the
%        collection.

%         References:
%         N. J. Higham, Accuracy and Stability of Numerical Algorithms,
%            Second edition, Society for Industrial and Applied Mathematics,
%            Philadelphia, PA, 2002; sec. 20.5.

  % Matrices from gallery:
  matrices = [%
  'cauchy  '; 'chebspec'; 'chebvand'; 'chow    ';
  'circul  '; 'clement '; 'condex  ';
  'cycol   '; 'dramadah'; 'fiedler ';
  'forsythe'; 'frank   '; 'gearmat '; 'grcar   ';
  'invhess '; 'invol   '; 'ipjfact '; 'jordbloc';
  'kahan   '; 'kms     '; 'krylov  '; 'lehmer  ';
  'lesp    '; 'lotkin  '; 'minij   '; 'moler   ';
  'orthog  '; 'parter  '; 'pei     '; 'prolate ';
  'randcolu'; 'randcorr'; 'rando   '; 'randsvd ';
  'redheff '; 'riemann '; 'ris     '; 'smoke   ';
  'toeppd  '; 'triw    ';];
  n_gall = length(matrices);
  
  % Other MATLAB matrices:
  matrices = [matrices;
  'hilb    '; 'invhilb '; 'magic   '; 'pascal  ';
  'rand    '; 'randn   ';];
  n_MATLAB = length(matrices);
  
  % Matrices from Matrix Computation Toolbox:
  matrices = [matrices;
  'augment '; 'gfpp    '; 'magic   '; 'makejcf ';
  'rschur  '; 'vand    '];
  n_mats = length(matrices);

  if nargin == 0
  
     rows = ceil(n_mats/5);
     temp = zeros(rows,5);
     temp(1:n_mats) = 1:n_mats;
  
     for i = 1:rows
        for j = 1:5
          if temp(i,j) == 0, continue, end
          fprintf(['%2.0f: ' sprintf('%s',matrices(temp(i,j),:)) '  '], ...
                  temp(i,j))
        end
        fprintf('\n')
     end
     fprintf('Matrices 1 to %1.0f are from MATLAB\.', n_MATLAB)
  
  elseif nargin == 1
     if k == 0
        A = length(matrices);
     elseif k > 0
        A = deblank(matrices(k,:));
     else
        % Version number and date of collection.
        A = 'Version 1.2, September 5 2002';
     end
  else
     if k <= n_gall
        A = eval( ['gallery(''' deblank(matrices(k,:)) ''',n)'] );
     else
        A = eval( [deblank(matrices(k,:)) '(n)'] );
     end
  end

end

function [A, nmats] = expm_testmats(k,n)
%EXPM_TESTMATS  Test matrices for matrix exponential.
%   [A, NMATS] = EXPM_TESTMATS(K,N) selects the K'th test matrix.
%    NMATS is the number of test matrices available.
%    N sets the dimension of those matrices for which the dimension
%    is variable.

  if nargout > 1, nmats = 32; A = 0; end
  if nargin < 1, return; end
  if nargin < 2, n = 8; end

  switch k
  
     case 1
     % \cite[Test 1]{ward77}.
     A = [4 2 0; 1 4 1; 1 1 4];
  
     case 2
     % \cite[Test 2]{ward77}.
     A = [29.87942128909879     .7815750847907159 -2.289519314033932
            .7815750847907159 25.72656945571064    8.680737820540137
          -2.289519314033932   8.680737820540137  34.39400925519054];
  
     case 3
     % \cite[Test 3]{ward77}.
     A = [-131 19 18;
          -390 56 54;
          -387 57 52];
  
     case 4
     % \cite[Test 4]{ward77}.
     A = gallery('forsythe',10,1e-10,0);
  
     case 5
     % \cite[p. 370]{naha95}.
     T = [1 10 100; 1 9 100; 1 11 99];
     A = T*[-0.001 0 0; 0 -1 0; 0 0 -100]/T;
  
     case 6
     % \cite[Ex.~2]{kela98}.
     A = [0.1 1e6; 0 0.1];
  
     case 7
     % \cite[p.~655]{kela98}.
     A = [0  3.8e3 0    0   0
          0 -3.8e3 1    0   0
          0 0     -1  5.5e6 0
          0 0      0 -5.5e6 2.7e7
          0 0      0   0   -2.7e7];
  
     case 8
     % \cite[Ex.~3.10]{dipa00}
     w = 1.3; x = 1e6; n = 8;
     A = (1/n) * [w*ones(n/2) x*ones(n/2)
                  zeros(n/2)  -w*ones(n/2)];
  
     case 9
     A = rosser;
     A = 2.05*A/norm(A,1);  % Bad case for expm re. cost.
  
     case 10
     A = [0 1e4;
          -1e4 0];  % exp = [cos(x) sin(x); - sin(x) cos(x)], x = 100;
  
     case 11
     A = 1e2*triu(randn(n),1);  % Nilpotent.
  
     case 12 % log of Cholesky factor of Pascal matrix. See \cite{edst03}.
     A = zeros(n); A(n+1:n+1:n^2) = 1:n-1;
  
     case 13 % \cite[p.~206]{kela89}
     A = [48 -49 50 49; 0 -2 100 0; 0 -1 -2 1; -50 50 50 -52];
  
     case 14 % \cite[p.~7, Ex I]{pang85}
     A = [0    30 1   1  1  1
         -100   0 1   1  1  1
          0     0 0  -6  1  1
          0     0 500 0  1  1
          0     0 0   0  0  200
          0     0 0   0 -15 0];
  
     case 15 % \cite[p.~9, Ex II]{pang85}
     % My interpretation of their matrix for arbitrary n.
     % N = 31 corresponds to the matrix in above ref.
     A = gallery('triw',n,1);  m = (n-1)/2;
     A = A - diag(diag(A)) + diag(-m:m)*sqrt(-1);
     for i = 1:n-1, A(i,i+1) = -2*(n-1)-2 + 4*i; end
  
     case 16 % \cite[p.~10, Ex III]{pang85}
     A = gallery('triw',n,1,1);
     A = A - diag(diag(A)) + diag(-(n-1)/2:(n-1)/2);
  
     case 17
     % \cite[Ex.~5]{kela89}.
     A = [0 1e6; 0 0];   % Same as case 6 but with ei'val 0.1 -> 0.
  
     case 18
     % \cite[(52)]{jemc05}.
     g = [0.6 0.6 4.0]; b = [2.0 0.75];
     A = [-g(1)       0    g(1)*b(1)
            0        -g(2) g(2)*b(2)
          -g(1)*g(3)  g(3) -g(3)*(1-g(1)*b(1))];
  
     case 19
     % \cite[(55)]{jemc05}.
     g = [1.5 0.5 3.0 2.0 0.4 0.03]; b = [0.6 7.0];
     A1 = [-g(5)     0      0
            0      -g(1)    0
          g(4)     g(4)   -g(3)];
     A2 = [-g(6)       0    g(6)*b(2)
            0        -g(2)  g(2)*b(1)
            0         g(4) -g(4)];
     A = [zeros(3) eye(3); A2 A1];
  
     case 20
     % \cite[Ex.~3]{kela98}.
     A = [-1 1e7; 0 -1e7];
  
     case 21
     % \cite[(21)]{mopa03}.
     Thalf = [3.8235*60*24 3.10 26.8 19.9]/60;  % Half lives in seconds/
     a = log(2)./Thalf;  % decay constant
     A = diag(-a) + diag(a(1:end-1),-1);
  
     case 22
     % \cite[(26)]{mopa03}.
     a1 = 0.01145;
     a2 = 0.2270;
     A = [-a1              0  0
           0.3594*a1     -a2  0
           0.6406*a1     a2  0];
  
     case 23
     % \cite[Table 1]{kase99}.
     a = [4.916e-18
          3.329e-7
          8.983e-14
          2.852e-13
          1.373e-11
          2.098e-6
          9.850e-10
          1.601e-6
          5.796e-8
          0.000];
     A = diag(-a) + diag(a(1:end-1),-1);
  
     case 24
         % Jitse Niesen sent me this example.
         lambda = 1e6 * 1i;
         mu = 1/2*(-1+sqrt(1+4*lambda));
         A = [ 0, 1; lambda, -1 ] - mu*eye(2);
  
     case 25 % Awad
  
         A = [1 1e17;0 1];
  
     case 26 % Awad
  
         b = 1e3; x = 1e10;
         A = [ 1-b/2   b/2 ; -b/2   1+b/2 ];
         A = [A          x*ones(2);
              zeros(2)       -A    ];
  
     case 27 % Awad
         b = 1e5;
         A = [ 1-b/2   b/2 ; -b/2   1+b/2 ];
  
     case 28 % Awad
         b = 1e3;
         A = [ 1-b/2   b/2 ; -b^4/2   1+b/2 ];
  
     case 29 % Edited by Prashanth: 29/12/2017. inserted the godunov_demo matrix 
     % p.10, S. K. Godunov, "Modern Aspects of Linear Algebra", AMS Vol 175, 1998.
     A = [  289,  2064,   336,  128,   80,    32,   16;
           1152,    30,  1312,  512,  288,   128,   32;
            -29, -2000,   756,  384, 1008,   224,   48;
            512,   128,   640,    0,  640,   512,  128;
           1053,  2256,  -504, -384, -756,   800,  208;
           -287,   -16,  1712, -128, 1968,   -30, 2032;
          -2176,  -287, -1565, -512, -541, -1152, -289];
  
     case 30
     % \cite[(14.17), p. 141]{trem05}.
     A = 10*[0 1 2; -0.01 0 3; 0 0 0];
  
     case 31
     A = triu(schur(gallery('invol',13),'complex'),1);
  
     case 32
     % \cite{kuda10}
     alpha = 1; beta = 1;  % No values are given in the paper, unfortunately.
     A = -eye(n) + alpha/2*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
     A(1,2) = beta; A(n,n-1) = beta;
  
  %    Easy example.
  %    case 13 % 7-by-7.  Idea from \cite[Matrix F, p. 81 ]{kags77a}
  %    e = -1; f = -1.1;
  %    A = blkdiag(compan([1, -4*e, 6*e^2, -4*e^3 e^4]), ...
  %                compan([1, -3*f, 3*f^2, -f^3]));
  
  %    Removed because too ill conditioned.
  %    case 8
  %    % \cite[Ex.~3]{dahi03}
  %    A = gallery('triw',4,2^(60)) - diag([17 17 2 2]);
  %    % A = A/1e4;  % Make less ill conditioned.
  
  %    case 8
  %    % My concoction.
  %    A = gallery('randsvd', 8, 1e14);
  %
  %    case 9
  %    % My concoction.
  %    A = randjorth(4,4,-1e8)/1e3;
  
  end

end

function [A, nmats] = logm_testmats(k,n)
%LOGM_TESTMATS  Test matrices for matrix logarithm.
%               [A, NMATS] = LOGM_TESTMATS(K,N) selects K'th test matrix.
%               NMATS is the number of test matrices available.
%               N sets the dimension of those matrix for which the dimension
%               is variable.

  if nargout > 1, nmats = 8; A = []; end
  if nargin < 1, return; end
  if nargin < 2, n = 8; end

  switch k
  
     case 1
     % \cite[Test 1]{casi01}.
     a = 0.1; b = 1e5; % b = 1e6; % Causes -ve ei'val error in logm_iss.
     A = (exp(a)/2)*[2+b b; -b 2-b];
  
     case 2
     % \cite[Test 2d]{casi01}.
     A = gallery(3);
  
     case 3
     % \cite[Ex. 4.4]{dipa00}.
     b = 1e6; c = 0.1; tol = 1e-8;
     % tol added so logm_x works.
     A = [exp(c) b*exp(c)
            0      exp(c) + tol];
  
     case 4
     % \cite[Ex. 6.3]{dmp06}.
     tol1 = 1e-1; tol2 = tol1*100;
     % tol added so logm_x works.
     A = [1+1e-7 1e5 1e4;
           0      1+tol1  1e5;
           0      0  1+tol2];
  
     % e=1e-8 in cases 5-8 might be harder but causes logm_iss to fail.
  
     case 5
     % Liable to give completely wrong result if log branches taken wrongly.
     e = 1e-7;
     temp = exp(sqrt(-1)*(pi-e));
     temp2 = exp(sqrt(-1)*(pi+e));
     A = [temp 1; 0 temp2];
  
     case 6
     % As case 5 with bigger (1,2) entry.
     e = 1e-7;
     temp = exp(sqrt(-1)*(pi-e));
     temp2 = exp(sqrt(-1)*(pi+e));
     A = [temp 1e3; 0 temp2];
  
     case 7
     % Gives rel_error 1e-9 for my log1p formula.
     e = 1e-7;
     temp = exp(sqrt(-1)*(pi-e));
     temp2 = temp*(1+e*i);
     A = [temp 1; 0 temp2];
  
     case 8
     % As case 7 with bigger (1,2) entry.
     e = 1e-7;
     temp = exp(sqrt(-1)*(pi-e));
     temp2 = temp*(1+e*i);
     A = [temp 1e3; 0 temp2];
  
  end

end

function [c,est] = funm_condest1(A,fun,fun_frechet,flag1,varargin)
%FUNM_CONDEST1  Estimate of 1-norm condition number of matrix function.
%    C = FUNM_CONDEST1(A,FUN,FUN_FRECHET,FLAG) produces an estimate of
%    the 1-norm relative condition number of function FUN at the matrix A.
%    FUN and FUN_FRECHET are function handles:
%      - FUN(A) evaluates the function at the matrix A.
%      - If FLAG == 0 (default)
%           FUN_FRECHET(B,E) evaluates the Frechet derivative at B
%              in the direction E;
%        if FLAG == 1
%           - FUN_FRECHET('notransp',E) evaluates the
%                    Frechet derivative at A in the direction E.
%           - FUN_FRECHET('transp',E) evaluates the
%                    Frechet derivative at A' in the direction E.
%    If FUN_FRECHET == @CS then the Frechet derivative is approximated
%    by the complex-step approximation (assuming A is real),
%    while if FUN_FRECHET == @FD then finite differences are used, the latter
%    being the default if FUN_FRECHET is empty.
%    More reliable results are obtained when FUN_FRECHET is supplied.
%    MATLAB'S NORMEST1 (block 1-norm power method) is used, with a random
%    starting matrix, so the approximation can be different each time.
%    C = FUNM_CONDEST1(A,FUN,FUN_FRECHET,FLAG,P1,P2,...) passes extra inputs
%    P1,P2,... to FUN and FUN_FRECHET.
%    [C,EST] = FUNM_CONDEST1(A,...) also returns an estimate EST of the
%    1-norm of the Frechet derivative.
%    Note: this function makes an assumption on the adjoint of the
%    Frechet derivative that, for f having a power series expansion,
%    is equivalent to the series having real coefficients.

fte_diff = 0; cstep = 0;
if nargin < 3 || isempty(fun_frechet)
   fte_diff = 1;
elseif isequal(fun_frechet,@CS)
   cstep = 1;
elseif isequal(fun_frechet,@FD)
   fte_diff = 1;
end

if nargin < 4 || isempty(flag1), flag1 = 0; end

n = length(A);
funA = feval(fun,A,varargin{:});
if fte_diff, tol = eps; end
if cstep, tol = eps^2; end

factor = norm(A,1)/norm(funA,1);

[est,v,w,iter] = normest1(@afun);
c = est*factor;

       %%%%%%%%%%%%%%%%%%%%%%%%% Nested function.
       function Z = afun(flag,X)
       %AFUN  Function to evaluate matrix products needed by NORMEST1.

       if isequal(flag,'dim')
          Z = n^2;
       elseif isequal(flag,'real')
          Z = isreal(A);
       else

          [p,q] = size(X); i = sqrt(-1);
          if p ~= n^2, error('Dimension mismatch'), end
          E = zeros(n);
          Z = zeros(n^2,q);
          for j=1:q

              E(:) = X(:,j);

              if isequal(flag,'notransp')

                 if fte_diff
                    t = sqrt(tol*norm(funA,1))/norm(E,1);
                    Y = (feval(fun,A+t*E,varargin{:}) - funA)/t;
                 elseif cstep
                    t = tol*norm(A,1)/norm(E,1);
                    Y = imag(feval(fun,A+t*i*E,varargin{:}))/t;
                 else
                    if flag1
                       Y = feval(fun_frechet,'notransp',E,varargin{:});
                    else
                       Y = feval(fun_frechet,A,E,varargin{:});
                    end
                 end

              elseif isequal(flag,'transp')

                 if fte_diff
                    t = sqrt(tol*norm(funA,1))/norm(E,1);
                    Y = (feval(fun,A'+t*E,varargin{:}) - funA')/t;
                 elseif cstep
                    t = tol*norm(A,1)/norm(E,1);
                    Y = imag(feval(fun,A'+t*i*E,varargin{:}))/t;
                 else
                    if flag1
                       Y = feval(fun_frechet,'transp',E,varargin{:});
                    else
                       Y = feval(fun_frechet,A',E,varargin{:});
                    end
                 end

              end

              Z(:,j) = Y(:);
          end

       end

       end

end
