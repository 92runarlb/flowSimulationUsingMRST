function [state] = fluxApprox(G, state, rock,fluid, bc,varargin)


   opt = struct('src', [], 'wells', [], ...
                'bcp',[],...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false, ...
                'Verbose',      mrstVerbose,...
                'gravity',      gravity(), ...
                'condition_number',false,...
                'pc_form','nonwetting',...
                'use_trans',false);

   opt = merge_options(opt, varargin{:});
   T = computeTrans(G, rock);
   g_vec   = opt.gravity;
   % If gravity is overriden, we cannot say anything about the effects of
   % gravity on rhs.
   grav = (norm(g_vec(1:G.griddim)) > 0) || isfield(G, 'grav_pressure');

   if all([~opt.MatrixOutput , ...
           isempty(opt.src)  , ...
           isempty(opt.bcp)  , ...
           isempty(opt.wells), ~grav]),
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   %% ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Setting up linear system...\t\t\t');
   t0 = ticif (opt.Verbose);

   % Preliminaries
   [neighborship, n_isnnc] = getNeighbourship(G, 'Topological', true);
   [cellNo, cf, cn_isnnc] = getCellNoFaces(G);
   nif    = size(neighborship, 1);
   ncf    = size(cf, 1);
   nc     = G.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;

   [mob, omega, rho] = dynamic_quantities(state, fluid);
   totmob = sum(mob, 2);

   % Compute effective (mobility-weighted) transmissibilities.
   [T, ft] = compute_trans(G, T, cellNo, cf, neighborship, totmob, opt);

   % Identify internal faces
   i  = all(neighborship ~= 0, 2);

   % Boundary conditions and source terms.
   hh = zeros(nif, 1);
   dF = false(nif, 1);
   [grav, ff] = deal(zeros(ncf, 1));
   [ff(~cn_isnnc), gg, hh(~n_isnnc), grav(~cn_isnnc), dF(~n_isnnc), dC] = ...
      computePressureRHS(G, omega, bc, opt.src);

   clear mob

   sgn = 2*(neighborship(cf, 1) == cellNo) - 1;
   j   = i(cf) | dF(cf);
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nif, 1]);
   

   p = state.cellMoments;

   %% ---------------------------------------------------------------------
   t0 = ticif (opt.Verbose);

   % Reconstruct face pressures and fluxes.
   fpress     =  ...
          accumarray(cf, (p(cellNo)+grav).*T, [nif, 1])./ ...
          accumarray(cf(:,1), T, [nif,1]);


   % Neumann faces
   b         = any(G.faces.neighbors==0, 2);
   fpress(b) = fpress(b) - hh(b)./ft(b);


   % Contribution from gravity
   %fg         = accumarray(cf, grav.*sgn, [nf, 1]);
   %fpress(~i) = fpress(~i) + fg(~i);

   % Dirichlet faces
   fpress(dF) = dC;


   % Sign for boundary faces
   sgn  = 2*(G.faces.neighbors(~i,2)==0)-1;
   ni   = neighborship(i,:);
   flux = -accumarray(find(i),  ft(i) .*(p(ni(:,2))-p(ni(:,1))-fg(i)), [nif, 1]);
   c    = sum(G.faces.neighbors(~i,:),2) ;
   fg  = accumarray(cf, grav, [nif, 1]);
   flux(~i) = -sgn.*ft(~i).*( fpress(~i) - p(c) - fg(~i) );
   %flux = -sgn.*ft((fpress(~i)-p(c)-grav));
   state.pressure(1 : nc) = p(1 : nc);
   state.flux             = flux;
   state.facePressure     = fpress;

   for k = 1 : nw,
      wc       = W(k).cells;
      dp       = norm(gravity()) * W(k).dZ*sum(rho .* W(k).compi, 2);
      state.wellSol(k).flux     = W(k).WI.*totmob(wc).*(p(nc+k) + dp - p(wc));
      state.wellSol(k).pressure = p(nc + k);
   end

   tocif(opt.Verbose, t0);
   
   
%    sol.flux = zeros(G.faces.num,1);
% neuF = strcmp(bc.type, 'flux');
% dirF = bc.face(~neuF);
% 
% 
% 
% neigh = G.faces.neighbors;
% intNeigh = all(neigh~=0,2);
% 
% %sgn      = 2*(neigh(dirF,1)==0) -1;
% dirNeigh = diag(neigh(dirF, 1 + (neigh(dirF,1)==0)));
% 
% [~,IA, hFdir] = intersect(dirNeigh, G.cells.faces);
% 
% 
% %hf = G.cells.faces(:,1);
% %hT = computeTrans(G, rock);
% %T = 1./accumarray(hf, 1./hT, [G.faces.num,1]);
% 
% % Fix normal sign
% sol.flux(bc.face(neuF)) = bc.value(neuF);
% sol.flux(dirF(IA))      = sgn(IA).*hT(hFdir).*(sol.pressure(dirNeigh(IA)) - bc.value(IA));
% sol.flux(intNeigh)      = T(intNeigh).*diff(sol.pressure(neigh(intNeigh,:)),1,2);

end



function [mob, omega, rho] = dynamic_quantities(state, fluid)
   [mu, rho] = fluid.properties(state);
   s         = fluid.saturation(state);
   kr        = fluid.relperm(s, state);

   mob    = bsxfun(@rdivide, kr, mu);
   totmob = sum(mob, 2);
   omega  = sum(bsxfun(@times, mob, rho), 2) ./ totmob;
end

%--------------------------------------------------------------------------

function [T, ft] = compute_trans(G, T, cellNo, cellFaces, neighborship, totmob, opt)
    niface = size(neighborship, 1);
    if opt.use_trans,  
      neighborcount = sum(neighborship > 0, 2);
      assert (numel(T) == niface, ...
             ['Expected one transmissibility for each interface ', ...
              '(=%d) but got %d.'], niface, numel(T));

      fmob = accumarray(cellFaces, totmob(cellNo), ...
                        [niface, 1]);
  
      fmob = fmob ./ neighborcount;
      ft   = T .* fmob;

      % Synthetic one-sided transmissibilities.
      th = ft .* neighborcount;
      T  = th(cellFaces(:,1));

   else

      % Define face transmissibility as harmonic average of mobility
      % weighted one-sided transmissibilities.
      %
      assert (numel(T) == numel(cellNo), ...
             ['Expected one one-sided transmissibility for each ', ...
              'half face (=%d), but got %d.'], numel(cellNo), numel(T));

      T  = T .* totmob(cellNo);
      ft = 1 ./ accumarray(cellFaces, 1 ./ T, [niface, 1]);

   end
end