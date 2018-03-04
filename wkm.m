function [C, S] = wkm (A)
%WKM  Evaluates the wave-kernel matrix functions.
%   [C, S] = WKM(A) evaluates the matrix functions C = COSH(SQRT(A)) and 
%   S = SINHC(SQRT(A)) at the square matrix A. For any nonzero scalar z, 
%   SINHC(z) = SINH(z)/z and SINHC(0) = 1.
%
%   C = WKM(A) evaluates only the matrix function C = COSH(SQRT(A)) and omit
%   the S = SINHC(SQRT(A)) computations.
%
%   See also EXPM, SQRTM, LOGM, FUNCTION_HANDLE.
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

  [m, s, As] = get_PadeOrder_Scaling(A); % Alg 5.1 of MIMS EPrint 2018.4 
  for ii = 1:length(As) % apply scaling
    As{ii} = As{ii}/(4^(ii*s));
  end

  [P, Q, dPxQ_PxdQ] = get_WaveKernels_diag_Pade_coeffs(m);
  P = [1;P];  Q = [1;Q];  dPxQ_PxdQ = [1;dPxQ_PxdQ];

  pA = EvalMatPolyPS(P, As);
  qA = EvalMatPolyPS(Q, As);
  I  = speye(length(A));

  % The scaling up process is intrinsically unstable.
  % For an ill conditioned problem there is little you can do to control
  % the amplification of errors in this stage.
  % 'm' and 's' were chosen such that 's' is minimal
  if nargout < 2 % the user wants just COSH(SQRT(A))
    C = qA\pA; % COSH(SQRT(A/4^s))
    for i = 1:s % scaling up using the double angle formula
      C = 2*C^2 - I;
    end
  else % compute both COSH(SQRT(A)) and SINHC(SQRT(A))
    dpxq_pxdqA  = EvalMatPolyPS(dPxQ_PxdQ ,As); % p'(A)*q(A) - p(A)*q'(A)

    [L,U] = lu(qA); % LU decomposition of qA 
    C = U\(L\pA); % COSH(SQRT(A/4^s))
    S = U\(L\(U\(L\dpxq_pxdqA))); % SINHC(SQRT(A/4^s))

    for i = 1:s % scaling up using the double angle formula
      S = S*C;
      C = 2*C^2 - I;
    end
  end

end

function [m, s, As] = get_PadeOrder_Scaling (A)
%GET_PADEORDER_SCALING returns the diagonal Pade approximation order 'm' and 
%   the required scaling 's' such that COSH(SQRT(A/4^s)) and SINHC(SQRT(A/4^s))
%   are approximated with small (~ 1/2^53) backward error and mixed 
%   forward-backward error, respectively.
%   The cell-array As stores the optimal powers of A to use in the 
%   Paterson--Stockmeyer (PS) algorithm.
%
%   This function is specialised for COSH(SQRT(A)) and SINHC(SQRT(A)).
%
%   We only consider m = {1,2,3,4,5,6,7,8,10,12,14,16,18,20} which are chosen 
%   based on the matrix multiplication cost of the PS algorithm to compute 
%   p_m(A), q_m(A) and p'_m(A)*q_m(A) - p_m(A)*q'_m(A).
%   The Higham--Tisseur (HT) 1-norm estimate (NORMEST1) is used when required.

  s = 0; % required scaling (initialisation)
  % Radii of discs (for considered 'm') within which the pade approximant has
  % an error bounded by the unit roundoff (1/2^53).
  % Ref. Table 4.2, Nadukandi--Higham paper.
  r01 = 9.42e-8;  r06 = 3;  r14 = 3;
  r02 = 2.31e-3;  r07 = 3;  r16 = 3;
  r03 = 9.14e-2;  r08 = 3;  r18 = 3;
  r04 = 6.66e-1;  r10 = 3;  r20 = 3;
  r05 = 2.36   ;  r12 = 3;
  
  As{1} = A;

  % Case m = 1. Powers of A stored: {1}
  % Estimate n02, n03 using HT::normest1
  n02 = normest1(@appAtoX,[],[],A,02)^(1/02);
  n03 = normest1(@appAtoX,[],[],A,03)^(1/03);
  a01 = max([n02;...
             n03]   );
  if a01 <= r01, m = 1; return; end
  % Available n{2,3}.

  % Case m = 2. Powers of A stored: {1--2}
  As{2} = A^2; n02 = norm(As{2},1)^(1/2);
  % Estimate n05 using HT::normest1
  n05 = normest1(@appAtoX,[],[],As{1},01,As{2},02)^(1/05);
  a02 = max([n02;...
             n05]   );
  if a02 <= r02, m = 2; return; end
  % Available n{2,3,5}.

  % Case m = 3. Powers of A stored: {1--3}
  As{3} = A*As{2}; n03 = norm(As{3},1)^(1/3); 
  % Estimate n04 and n07 using HT::normest1
  n04 = normest1(@appAtoX,[],[],As{1},01,As{3},01)^(1/04);
  n07 = normest1(@appAtoX,[],[],As{1},01,As{3},02)^(1/07);
  a03 = min(max([n02, n03;...
                 n07, n04]   ));
  if a03 <= r03, m = 3; return; end
  % Available n{2,3,4,5,7}.

  % Case m = 4. Powers of A stored: {1--4}
  As{4} = A*As{3}; n04 = norm(As{4},1)^(1/4); 
  % Estimate n09 using HT::normest1
  n09 = normest1(@appAtoX,[],[],As{1},01,As{4},02)^(1/09);
  a04 = min([a03, max([n02, n03;...
                       n09, n05]   )]);
  if a04 <= r04, m = 4; return; end
  % Available n{2,3,4,5,7,9}.

  % Case m = 5. Powers of A stored: {1--4}
  % Estimate n11 using HT::normest1
  n11 = normest1(@appAtoX,[],[],As{3},01,As{4},02)^(1/11);
  a05 = min([a04, max([n02;...
                       n11]   )]);
  if a05 <= r05, m = 5; As{5} = A*As{4}; return; end % Note: A^{1--5} returned 
  % Available n{2,3,4,5,7,9,11}.

  % Case m = 6. Powers of A stored: {1--4}
  % Estimate n13 using HT::normest1
  n13 = normest1(@appAtoX,[],[],As{1},01,As{4},03)^(1/13);
  a06 = min([a05, max([n02, n03, n04;...
                       n13, n07, n05]   )]);
  if a06 <= r06, m = 6; As{5} = A*As{4}; As{6} = A*As{5}; return; end
                                                     % Note: A^{1--6} returned
  % Available n{2--5,7,9,11,13}.

  % Case m = 7. Powers of A stored: {1--4}
  % Estimate n08, n15 using HT::normest1
  n08 = normest1(@appAtoX,[],[],         As{4},02)^(1/08);
  n15 = normest1(@appAtoX,[],[],As{3},01,As{4},03)^(1/15);
  a07 = min([a06, max([n02, n03;...
                       n15, n08]   )]);
  if a07 <= r07, m = 7; return; end
  % Available n{2--5,7--9,11,13,15}.

  % Case m = 8. Powers of A stored: {1--4}
  % Estimate n17 using HT::normest1
  n17 = normest1(@appAtoX,[],[],As{1},01,As{4},04)^(1/17);
  a08 = min([a07, max([n02;...
                       n17]   )]);
  if a08 <= r08, m = 8; return; end
  % Available n{2--5,7--9,11,13,15,17}.

  % Case m = 10. Powers of A stored: {1--5}
  As{5} = A*As{4}; n05 = norm(As{5},1)^(1/5); 
  % Estimate n10, n21 using HT::normest1
  n06 = normest1(@appAtoX,[],[],As{1},01,As{5},01)^(1/06);
  n10 = normest1(@appAtoX,[],[],         As{5},02)^(1/10);
  n21 = normest1(@appAtoX,[],[],As{1},01,As{5},04)^(1/21);
  a10 = min([a08, max([n02, n03, n03, n04, n05;...
                       n21, n10, n11, n07, n06]   )]);
  if a10 <= r10, m = 10; return; end
  % Available n{2--11,13,15,17,21}.

  % Case m = 12. Powers of A stored: {1--6}
  As{6} = A*As{5}; n06 = norm(As{6},1)^(1/6); 
  % Estimate n25 using HT::normest1
  n25 = normest1(@appAtoX,[],[],As{1},01,As{6},04)^(1/25);
  a12 = min([a10, max([n02, n03, n04, n05;...
                       n25, n13, n09, n07]   )]);
  if a12 <= r12, m = 12; return; end
  % Available n{2--11,13,15,17,21,25}.

  % Case m = 14. Powers of A stored: {1--7}
  As{7} = A*As{6}; n07 = norm(As{7},1)^(1/7); 
  % Estimate n14, n29 using HT::normest1
  n14 = normest1(@appAtoX,[],[],         As{7},02)^(1/14);
  n29 = normest1(@appAtoX,[],[],As{1},01,As{7},04)^(1/29);
  a14 = min([a12, max([n02, n03, n05;...
                       n29, n14, n08]   )]);
  if a14 <= r14, m = 14; return; end
  % Available n{2--11,13--15,17,21,25,29}.

  % Case m = 16. Powers of A stored: {1--8}
  As{8} = A*As{7}; n08 = norm(As{8},1)^(1/8); 
  % Estimate n16, n33 using HT::normest1
  n16 = normest1(@appAtoX,[],[],         As{8},02)^(1/16);
  n33 = normest1(@appAtoX,[],[],As{1},01,As{8},04)^(1/33);
  a16 = min([a14, max([n02, n03, n03, n04, n05, n06;...
                       n33, n16, n17, n11, n09, n07]   )]);
  if a16 <= r16, m = 16; return; end
  % Available n{2--11,13--17,21,25,29,33}.

  % Case m = 18. Powers of A stored: {1--9}
  As{9} = A*As{8}; n09 = norm(As{9},1)^(1/9); 
  % Estimate n19, n37 using HT::normest1
  n19 = normest1(@appAtoX,[],[],As{1},01,As{9},02)^(1/19);
  n37 = normest1(@appAtoX,[],[],As{1},01,As{9},04)^(1/37);
  a18 = min([a16, max([n02, n03, n04;...
                       n37, n19, n13]   )]);
  if a18 <= r18, m = 18; return; end
  % Available n{2--11,13--17,19,21,25,29,33,37}.

  % Case m = 20. Powers of A stored: {1--10}
  As{10} = A*As{9}; n10 = norm(As{10},1)^(1/10); 
  % Estimate n20, n41 using HT::normest1
  n20 = normest1(@appAtoX,[],[],         As{10},02)^(1/20);
  n41 = normest1(@appAtoX,[],[],As{1},01,As{10},04)^(1/41);
  a20 = min([a18, max([n02, n03, n05;...
                       n41, n20, n11]   )]);
  if a20 <= r20, m = 20; return; end
  % Available n{2--11,13--17,19--21,25,29,33,37,41}.

  % At this point A needs to be scaled.
  % The scaling up transformation is an unstable process (errors will always
  % amplify!). So we choose the smallest order 'm' which has the same scaling
  % as the highest order (m = 20) considered it in our algorithm.
  % The multiplication count mu_* jumps from m to m+1 only for {2,6,7,20}
  % m = 2 is excluded as we will have to scale down by 4^6 to enter this space
  % starting from the space corresponding to m = 6.
  s20 = ceil(0.5*log2(a20/r20));
  s07 = ceil(0.5*log2(a07/r07));
  s06 = ceil(0.5*log2(a06/r06));
  ss = [s06,s07,s20];
  mm = [ 06, 07, 20];
  kk = find(ss==ss(end),1,'first');
  m = mm(kk); s = ss(kk);
end

function [P, Q, dPxQ_PxdQ] = get_WaveKernels_diag_Pade_coeffs(m)
% P has the coeffs of x.^[1,2,...,m] in the numer (p) of the diag Pade approx.
% Q has the coeffs of x.^[1,2,...,m] in the denom (q) of the diag Pade approx.
% dPxQ_PxdQ has the coeffs of x.^[1,2,...,2*(m-1)] of (p´*q-p*q´)
%
% We only consider m = {1,2,4,6,9,12,16,20} which are chosen based the
% multiplication counts in the Paterson--Stockmeyer algorithm

  if ~ismember(m,[1,2,3,4,5,6,7,8,10,12,14,16,18,20])
    errMsg = ['We only consider m = {1,2,3,4,5,6,7,8,10,12,14,16,18,20} ',...
              'which are chosen based the multiplication counts in the',...
              'Paterson--Stockmeyer algorithm.'];
    error(errMsg);
  end

  if m == 01
    PQ = [
     4.16666666666666666667e-01 -8.33333333333333333333e-02
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [];
    return;
  end

  if m == 02
    PQ = [
     4.56349206349206349206e-01 -4.36507936507936507937e-02;
     2.07010582010582010582e-02  8.59788359788359788360e-04
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     7.93650793650793650794e-02;
    -2.59196271101033005795e-03
    ];
    return;
  end

  if m == 03
    PQ = [
     4.70595788392398561890e-01 -2.94042116076014381099e-02;
     2.73882896764252696456e-02  4.23728813559322033898e-04;
     3.72342268528709206675e-04 -3.23554348978077791637e-06
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [ 
     1.07858243451463790447e-01;
     2.43994751516101180513e-04;
    -3.77031908592849806190e-05;
     4.92776300080888626500e-07
    ];
    return;
  end

  if m == 04
    PQ = [
     4.77862206485004435256e-01 -2.21377935149955647437e-02; 
     3.08424350543706588484e-02  2.44665145201774553641e-04;
     5.87146544419631577915e-04 -1.66685394532939013794e-06;
     3.42184348618219314005e-06  6.23754467948680855643e-09
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.22391079636675537179e-01; 
     1.93348102045094054102e-03;
    -2.14815623047849010061e-05;
    -8.22682502097489228051e-08;
     2.57929906691482233200e-09;
    -1.87321322389324023517e-11
    ];
    return;
  end

  if m == 05
    PQ = [
     4.82257639391715124592e-01 -1.77423606082848754082e-02;
     3.29535721869666064782e-02  1.58085824442377515588e-04;
     7.27759945435805751819e-04 -9.06830329068752856247e-07;
     6.29636983550361826453e-06  3.45585819173468945201e-09;
     1.90341764928050555827e-08 -7.28729462656933528197e-12
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.31181945450096915850e-01;
     3.05017613941088208255e-03;
     4.43897681045943405166e-07;
    -1.99996171713876096116e-07;
     8.52360818498428687663e-10;
     3.04550663998886190882e-12;
    -4.78294695708655551606e-14;
     2.23325833649485958954e-16
    ];
    return;
  end

  if m == 06
    PQ = [
     4.85199972661551760971e-01 -1.48000273384482390290e-02; 
     3.43768869609038307826e-02  1.10233963461283630431e-04;
     8.26800030322268821585e-04 -5.38034528585256309436e-07;
     8.57191078523543125398e-06  1.85246267649409376551e-09;
     4.01155018648636115759e-08 -4.37328468117306375806e-12;
     7.08970168259802379018e-11  5.70462482331598949567e-15
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.37066611989770188609e-01; 
     3.83949962332530285053e-03;
     2.06580317259160000751e-05;
    -1.46106428324699763105e-07;
    -3.57077906032989126016e-10;
     4.62692701812690531637e-12;
    -1.09175575544553107764e-14;
    -3.35700098598898658852e-17;
     3.29738169984792220149e-19;
    -1.07779345072801629740e-21
    ];
    return;
  end

  if m == 07
    PQ = [
     4.87306595576691735341e-01 -1.26934044233082646588e-02;
     3.54011084410287894873e-02  8.11439860162551500248e-05;
     9.00226169980103686248e-04 -3.42860945735083535761e-07;
     1.03824770060651518485e-05  1.04912570735058194448e-09;
     5.96929205775457206094e-08 -2.36439612447397369460e-12;
     1.68552697695253142695e-10  3.70793813894712783488e-15;
     1.89098530697112730748e-13 -3.19910136684177567982e-18
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.41279857820050137349e-01;
     4.42560901345008389934e-03;
     3.80119972517386602820e-05;
    -2.64773083090988509726e-08;
    -8.08504268859985495261e-10;
     1.79042953703987851395e-12;
     6.22971906171734805756e-15;
    -3.78467872335860195495e-17;
     5.77048525027890798716e-20;
     1.49890352877395969478e-22;
    -1.02456051725428806112e-24;
     2.48076563914488288997e-27
    ];
    return;
  end

  if m == 08
    PQ = [
     4.88888897248517529962e-01 -1.11111027514824700378e-02;
     3.61732981204982696122e-02  6.21828295728379644459e-05;
     9.56786592467251645756e-04 -2.31096562953307114915e-07;
     1.18455335768711933388e-05  6.30479395727398407037e-10;
     7.70501557277950375985e-08 -1.30891912647731040717e-12;
     2.72636635527408396704e-10  2.04824876067076200586e-15;
     5.00875072274192971615e-13 -2.24601611287627073642e-18;
     3.78873545761552646992e-16  1.34799319033022574798e-21
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.44444461163701726591e-01;
     4.87745467967218723263e-03;
     5.26873303945067212149e-05;
     1.14674040864215575264e-07;
    -7.17034097025900028196e-10;
    -9.40792829834550356617e-13;
     8.11096373069659000554e-15;
    -8.29831212467049273376e-18;
    -3.61676563431636771388e-20;
     1.36564487184606825542e-22;
    -1.46962158513450332464e-25;
    -3.22141450008286618700e-28;
     1.63405976972618892159e-30;
    -3.05226455030956828021e-33
    ];
    return;
  end

  if m == 10
    PQ = [
     4.91106986599988200226e-01 -8.89301340001179977448e-03;
     3.72599794222139249868e-02  3.98194555531582074214e-05;
     1.03813747429589849624e-03 -1.18917369077839091929e-07;
     1.40501291242156680691e-05  2.64025802082627561944e-10;
     1.05493807230919319341e-07 -4.58696545637242442440e-13;
     4.70194661582235409366e-10  6.36659342248439102219e-16;
     1.27676037568902313278e-12 -7.03523926883094778239e-19;
     2.08230645663301758256e-15  5.98544834764512756002e-22;
     1.88771816469882461284e-18 -3.58304358729216245760e-25;
     7.38143374434111147734e-22  1.16212110636046590178e-28
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.48880639866643067118e-01;
     5.52772013176850592112e-03;
     7.57038435105662255081e-05;
     3.96014605818256075809e-07;
     3.37391306019663000899e-10;
    -2.52596281803587317521e-12;
    -1.61196139215783919459e-15;
     1.25284255435482489481e-17;
    -7.10181013998550610510e-21;
    -3.18457651897483060758e-23;
     6.61554691692825574314e-26;
    -1.86048016314971592245e-29;
    -1.28504974690431045346e-31;
     2.55100983459136823788e-34;
    -1.59915863901418206311e-37;
    -2.54092387333051452575e-40;
     7.99290703067178143955e-43;
    -9.67711401264977146157e-46
    ];
    return;
  end

  if m == 12
    PQ = [
     4.92587279209017505726e-01 -7.41272079098249427447e-03; 
     3.79879567242225230075e-02  2.76504530471034781203e-05;
     1.09378184121228565585e-03 -6.89079092177106580279e-08;
     1.56239183592874923644e-05  1.28345044309899027197e-10;
     1.27322177921172901113e-07 -1.89126792011344516549e-13;
     6.40251637006388709609e-10  2.27282688272306393470e-16;
     2.07657676592609551284e-12 -2.25743106635825220572e-19;
     4.43118746341460039135e-15  1.85097767312833773626e-22;
     6.20087059388646402719e-18 -1.23064989642292511177e-25;
     5.50124589872884771927e-21  6.34820683059209243752e-29;
     2.82319391476246355098e-24 -2.30187037525999586048e-32;
     6.43213195108112369476e-28  4.48281974192792169819e-36
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.51841225084701678118e-01; 
     5.97267573862510633365e-03;
     9.26944951714019036828e-05;
     6.43672345146636200104e-07;
     1.84186358492899498765e-09;
     3.43108710303279504200e-13;
    -6.52633029857525654463e-15;
    -1.58699275991808639391e-18;
     1.65115593683121504481e-20;
    -6.69906382811347323786e-24;
    -2.45590047702693055690e-26;
     3.31666178549789755269e-29;
     2.08295507042780117947e-33;
    -4.57141076024292340727e-35;
     4.81492814830992946436e-38;
    -3.04557803229741832314e-42;
    -4.84090269579285938620e-44;
     6.09941284017962581630e-47;
    -2.52564132422809337396e-50;
    -3.00189847533699367111e-53;
     6.46856410684045648395e-56;
    -5.49236068086897227197e-59
    ];
    return;
  end

  if m == 14
    PQ = [
     4.93645267774812775188e-01 -6.35473222518722481216e-03;
     3.85096098831950346493e-02  2.03093291219803886870e-05;
     1.13421964902531719458e-03 -4.33950417608548166695e-08;
     1.68001643240393784682e-05  6.94760100225552580593e-11;
     1.44399670533707991717e-07 -8.84682415789032264478e-14;
     7.82765253952532367326e-10  9.27389721699486294757e-17;
     2.81742087010673018904e-12 -8.16266853952786903327e-20;
     6.94610544982982292472e-15  6.08829813290988628424e-23;
     1.19144634041512783656e-17 -3.84762578658913729160e-26;
     1.42324042034667201668e-20  2.03861301340361450639e-29;
     1.16458267371902757737e-23 -8.83154117613070127533e-33;
     6.24845934600617278208e-27  2.97589145858394444659e-36;
     1.99029054533904151751e-30 -7.03082310109626404992e-40;
     2.87004069849253305441e-34  8.85663883191590192804e-44
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.53957202216292217042e-01;
     6.29609053816871881771e-03;
     1.05655797195868484427e-04;
     8.52641269843079369796e-07;
     3.40665968288961152249e-09;
     5.70151377463056838667e-12;
    -1.22191523504936569790e-15;
    -1.26810449853629352814e-17;
     2.95714044481292369239e-24;
     1.83476753093125298103e-23;
    -6.52515015868671271394e-27;
    -1.67014448664093328094e-29;
     1.64729728070888260012e-32;
     3.29035033598608860289e-36;
    -1.55138683232450958050e-38;
     9.69016385397350671549e-42;
     2.71121532238622561250e-45;
    -8.46148971047378296994e-48;
     5.56766382549588144564e-51;
     8.48925117608400796109e-56;
    -3.20374588882016687985e-57;
     2.82412156622783258321e-60;
    -8.34315116057423528488e-64;
    -7.64562710984586980469e-67;
     1.20275793283344634471e-69;
    -7.56120659493851822564e-73
    ];
    return;
  end

  if m == 16
    PQ = [
     4.94439056907362683033e-01 -5.56094309263731696716e-03; 
     3.89017396763001676381e-02  1.55445559521594550289e-05;
     1.16492614629915736318e-03 -2.90583725897129278474e-08;
     1.77112566301241556388e-05  4.07565991178991501926e-11;
     1.58052204873296607643e-07 -4.55950124884293839644e-14;
     9.02075105901299341770e-10  4.22002345457221655751e-17;
     3.47899664449148260335e-12 -3.30588153184829998294e-20;
     9.40013681951646454826e-15  2.22220333025354531635e-23;
     1.82090860867144657255e-17 -1.29036176772402837261e-26;
     2.56073156346485354075e-20  6.47448002576229292634e-30;
     2.62135872079977404686e-23 -2.78858586073525614238e-33;
     1.93687750015188625112e-26  1.01502639517306059317e-36;
     1.00874089006444433090e-29 -3.03257238945512428351e-40;
     3.52449764220010685183e-33  7.04949681389823647827e-44;
     7.44452146203099002063e-37 -1.14561068673957313905e-47;
     7.21986797815512209353e-41  9.86637610469039945661e-52
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.55544780481392032732e-01; 
     6.54169883577143060824e-03;
     1.15834845338774752779e-04;
     1.02793121525557786983e-06;
     4.88863904856669522967e-09;
     1.21693344988986581691e-11;
     1.22515788005983096393e-14;
    -6.33346391355583751016e-18;
    -1.89389701511029043138e-20;
     3.10116932537540334559e-24;
     1.71270746162635696291e-26;
    -6.01194104659642451643e-30;
    -9.99708675584284810834e-33;
     7.91154088963596196034e-36;
     1.82866963687023356641e-39;
    -5.04124930065998079457e-42;
     1.99494602072715395490e-45;
     1.02956863278156841526e-48;
    -1.45769801952815888835e-51;
     5.25455981938977591290e-55;
     1.97056747379498959611e-58;
    -3.15377594842986100222e-61;
     1.44021833238071240593e-64;
     6.42125975151928503752e-69;
    -5.08715068811374706243e-71;
     3.33387360120885932278e-74;
    -7.40344169715627081816e-78;
    -5.37003185669582064602e-81;
     6.44896679629185924723e-84;
    -3.12324055580141650840e-87
    ];
    return;
  end

  if m == 18
    PQ = [
     4.95056602570866656567e-01 -4.94339742913334343337e-03;
     3.92072463405490675211e-02  1.22783884490725710962e-05;
     1.18903279400442783447e-03 -2.03962284413635812648e-08;
     1.84371843833165692984e-05  2.54399237915857711644e-11;
     1.69185258384480581619e-07 -2.53508363631168139213e-14;
     1.00265065908197122777e-09  2.09602555634318419442e-17;
     4.06255818367578215894e-12 -1.47350442537453447709e-20;
     1.16995654112910557499e-14  8.94977808373350845898e-24;
     2.45953919004753807041e-17 -4.74360665940493169516e-27;
     3.84191842798907352157e-20  2.20515737376521504270e-30;
     4.50530058632609126369e-23 -8.99551343471547239650e-34;
     3.97940113396698267315e-26  3.20608781830183809052e-37;
     2.63615762284087612109e-29 -9.88109647537497665020e-41;
     1.29167644645363605865e-32  2.58475254816981319657e-44;
     4.55151595708924619224e-36 -5.55905166409282769139e-48;
     1.09381485315185351576e-39  9.29597967337718844155e-52;
     1.61048049731009480505e-43 -1.08357182808314752678e-55;
     1.10037082007243738943e-47  6.65974480267427152187e-60
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.56779871808399979800e-01 ;
     6.73452814532939281362e-03 ;
     1.24026214069225085105e-04 ;
     1.17572146580683277765e-06 ;
     6.24263421711736750190e-09 ;
     1.89680058047756825775e-11 ;
     3.09975648066534690099e-14 ;
     1.85326524891222202823e-17 ;
    -1.54186569114995847509e-20 ;
    -2.20976097800047915805e-23 ;
     6.37939641124043626460e-27 ;
     1.34251359857872497156e-29 ;
    -4.94207282829615675833e-33 ;
    -5.24870509134322242806e-36 ;
     3.60868083072011684291e-39 ;
     7.60572064024325295271e-43 ;
    -1.56105127033746666743e-45 ;
     4.15311293564344029032e-49 ;
     2.81585464502800545762e-52 ;
    -2.48832450988386215242e-55 ;
     4.66970002245417705775e-59 ;
     3.88188761369491242064e-62 ;
    -3.08701875727182076206e-65 ;
     7.15895679250176662946e-69 ;
     2.99744432879618085696e-72 ;
    -3.10834779222801689230e-75 ;
     1.05111861669139701514e-78 ;
     5.70316501353720130471e-83 ;
    -2.40328796088401084998e-85 ;
     1.22116207579655356904e-88 ;
    -2.11697441445246733162e-92 ;
    -1.24274313935058356306e-95 ;
     1.17779879728191501621e-98 ;
    -4.52973946650432237260e-102
    ];
    return;
  end

  if m == 20
    PQ = [
     4.95550725310104691213e-01 -4.44927468989530878710e-03;
     3.94519720000664486061e-02  9.94267834743633299251e-06;
     1.20845892338763775861e-03 -1.48592626647639751251e-08;
     1.90289044357286798906e-05  1.66814073975093532181e-11;
     1.78422582854783415129e-07 -1.49767782734111536863e-14;
     1.08821879100262654473e-09  1.11764772275502320155e-17;
     4.57605543189248276115e-12 -7.11136566507815556850e-21;
     1.38139352173452567039e-14  3.92553884103855895433e-24;
     3.08093360957904387877e-17 -1.90219607453177365090e-27;
     5.18172263059345333911e-20  8.15247725770290850404e-31;
     6.66488029053424278249e-23 -3.10254521228588909907e-34;
     6.61254900257269223436e-26  1.04904950663539491282e-37;
     5.07815161543870234510e-29 -3.14222417373734597199e-41;
     3.01262393773894531792e-32  8.27820118585584293928e-45;
     1.36900277252276284691e-35 -1.89407509644464570492e-48;
     4.68345205371052735329e-39  3.68680532543878828024e-52;
     1.16944461421110566722e-42 -5.90333503651328747265e-56;
     2.01626920797831917579e-46  7.34110894427222281046e-60;
     2.15169575011094980063e-50 -6.34612738635055971529e-64;
     1.07445780284307360596e-54  2.88023636336076489485e-68
    ];
    P = PQ(:,1);  Q = PQ(:,2); clear PQ;
    dPxQ_PxdQ = [
     1.57768117286876049092e-01 ; 
     6.88992317199591273312e-03 ;
     1.30753493301707644806e-04 ;
     1.30139053329231036129e-06 ;
     7.46203672269213523302e-09 ;
     2.56885842405466670066e-11 ;
     5.27239053526996045392e-14 ;
     5.85062600230451311717e-17 ;
     1.91013216178235626001e-20 ;
    -2.54846043836202057541e-23 ;
    -2.03272884179782235741e-26 ;
     8.19837199181104493886e-30 ;
     8.83565991555852870213e-33 ;
    -3.53457793024795070420e-36 ;
    -2.40627323362387115804e-39 ;
     1.54038880798727718958e-42 ;
     2.60479203402165919109e-46 ;
    -4.57907323346905112286e-49 ;
     8.76866548251856889023e-53 ;
     6.71725701185369928184e-56 ;
    -4.19644134586301468120e-59 ;
     3.38337138753595899292e-63 ;
     6.09218348221504120444e-66 ;
    -3.01959978186878209103e-69 ;
     2.63909366125546003766e-73 ;
     3.53333334275696263483e-76 ;
    -1.88383982128654656121e-79 ;
     3.00942635161578663386e-83 ;
     1.31295892815380441600e-86 ;
    -9.79270036480435750958e-90 ;
     2.56344908381699057693e-93 ;
     1.42548297802988123455e-97 ;
    -3.95732371608395342949e-100;
     1.60814921972033499170e-103;
    -2.23970798999931577650e-107;
    -1.08429946239780730830e-110;
     8.32151958220241778168e-114;
    -2.60320768609177097679e-117
    ];
    return;
  end
end

function P = EvalMatPolyPS(c, As)
% Evalaute a polynomial at a matrix using the Paterson--Stockmeyer algorithm
% Usage: EvalMatPolyPS([c0,c1,...,cm], {A,A^2,...,A^s})
%
% c: coefficients of polynomial in ascending order.
% As: cellarray containing all stored powers of the matrix in ascending order.
% P: polynomial evaluated at the matrix argument using PS method
%
% IMP: MATLAB uses a system commonly called "copy-on-write" to avoid making a
% copy of the input argument inside the function workspace until or unless you
% modify the input argument.

  if ~iscell(As)
    errMsg = [' ',inputname(2),' must be a cellarray. ',...
              ' ','Consider passing {',inputname(2),'}',char(10),...
              ' ','Syntax: EvalMatPolyPS(c, As)',char(10),...
              ' ','c: coefficients of polynomial in ascending order.',...
              char(10),...
              ' ','As: cellarray containing all stored powers of the',...
              ' matrix in ascending order.'];
    error(errMsg);
  end

  m = length(c)-1; % degree of the polynomial
  s = length(As); % As := {A,A^2,...,A^s}; also indication of chunk size
  r = floor(m/s); % => r chunks of size s and a possible remainder chunk
  n = length(As{1});

  % The Paterson--Stockmeyer algorithm looks like this:
  % P = c(s*0+1)*I 
  %   +           c(s*0+2)*A^(1) + c(s*0+3)*A^(2) + ... + c(s*1+1)*A^(s)
  %   +   A^(s)*{ c(s*1+2)*A^(1) + c(s*1+3)*A^(2) + ... + c(s*2+1)*A^(s)
  %   +   A^(s)*[ c(s*2+2)*A^(1) + c(s*2+3)*A^(2) + ... + c(s*3+1)*A^(s)
  %   +   ...
  %   +   A^(s)*( c(s*q+2)*A^(1) + c(s*q+3)*A^(2) + ... + c(s*r+1)*A^(s)
  %   +   A^(s)*{ c(s*r+2)*A^(1) + ... + c(m+1)*A^(m-r*s) }             )...]};
  % The strategy is similar to static scheduling of loops in parallel computing

  P = zeros(n);
  for jj = (m-r*s):-1:1
    % Enter the loop => the (r+1)th chuck exists and is the innermost
    P = P + c(jj+r*s+1)*As{jj};
  end
  for ii = r:-1:1 % Horner's method applied at the chunk level 
    P = As{s}*P;
    for jj = s:-1:1 % c-weighted accumulation of the (ii)th chunk
      P = P + c(ii*s-s+jj+1)*As{jj};
    end
  end
  P = P + c(1)*speye(n);

  % =======================================================================
  % The following code suits better for MATLAB versions >= R2015a
  % Here @mtimes is defined to work with matrices and cellarrays.
  % So inner loops get vectorised and have a succinct representation.
  % -----------------------------------------------------------------------
  %% idx = (r*s+2):(m+1); % coeff indices of the innermost (r+1)th chunk
  %% P = c(r*s+1)*speye(n) + c(idx)*As(1:length(idx))'; % innermost chunk
  %% for ii = r:-1:1 % Horner's method applied at the chunk level 
  %%   idx = (ii-1)*s + [1:s]; % indices of the (ii)th chunk
  %%   % IMP: extra parenthesis coerse evaluation order.
  %%   P = c(idx(1))*speye(n) + (c(idx(2:end))*As(1:end-1)' + (As(end)*P));
  %% end
  % -----------------------------------------------------------------------
end

function Y = appAtoX(flag,X,varargin)
% Usage:
% normest1(@appAtoX,[],[],      A4,05); % estimate norm(A^20,1) where A4 = A^4.
% normest1(@appAtoX,[],[],A ,01,A4,10); % estimate norm(A^41,1) where A4 = A^4.
%
% This functions has a nonstandard signature.
% It is designed for use by Higham--Tisseur's normest1.m
% varargin = {A1,m1,A2,m2,...,As,ms}
% You supply A1,A2,...,Am and apply them m1,m2,...,ms times to X respectively.
%
% flag can be: 'dim', 'real', 'notransp' or 'transp'
%
% if flag = 'dim'       then Y = max(size(A))
% if flag = 'real'      then Y = isreal(A)
% if flag = 'notransp'  then Y = (A^m)*X 
% if flag = 'transp'    then Y = (X'*(A^m))'

  if mod(length(varargin),2) ~= 0
    errmsg = ['Syntax error.\n',...
              '  normest1(@appAtoX,[],[],A1,m1,A2,m2,...,As,ms)'];
    error(errmsg);
  end 
  ss = length(varargin)/2;
  A = cell(ss,1);
  m = ones(ss,1);
  mIsInteger = true; AIsReal = true; mIsPositive = true;
  for i = 1:ss
    m(i) = varargin{2*i};
    A{i} = varargin{2*i-1};
    mIsInteger = mIsInteger & (m(i) == floor(m(i)));
    AIsReal = AIsReal & isreal(A{i});
    mIsPositive = mIsPositive & (m(i) >= 0);
  end

  if isequal(flag,'dim'), Y = max(size(A{1})); return; end
  if isequal(flag,'real'), Y = AIsReal; return; end
  if ~mIsInteger
    error('Noninteger application of a matrix is not contemplated.');
  end
  if ~mIsPositive
    error('This function is not meant to be used with negative powers!');
  end

  if isequal(flag,'transp')
    Y = X';
    for i = 1:ss
      for j = 1:m(i), Y = Y*A{i}; end
    end
    Y = Y';
  else
    Y = X;
    for i = 1:ss
      for j = 1:m(i), Y = A{i}*Y; end
    end
  end

end
