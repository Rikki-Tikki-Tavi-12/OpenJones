classdef BallBearing < handle
    %BALLBEARING Class approximating the behaviour of a ball bearing
    %   Uses general Jones-Harris theory, assumes rotating inner race
    %   and stationary outer race.
    
%     Copyright (C) 2015 Samuel Sudhof
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.    
    
    properties
       geometry;
       % Struct containing the geometry relevant to the bearing dynamics
       % calculation. Needs to be set using the method setGeometry
       % (inputs)
       % D           ball diameter                         [m]
       % d_m         pitchline diameter                    [m]
       % alpha_free  loadless contact angle aka alpha_0    [deg]
       % r_i         inner raceway groove radius           [m]
       % r_o         outer raceway groove radius           [m]
       % z           number of balls in the bearing        [-]
       % Delta_P_D   ball clearance according to Jones not used by Harris
       %                can safely be set to 0;
       % [psi_0]     angle location of the first ball      [deg]
       % (calculated values)
       % psi         vector of the balls angular positions [deg]
       %      - setGeometry requires only the position of the first ball
       %                psi_0 to generate this vector -
       % f_i         r_i/D                                 [-]
       % f_o         r_o/D                                 [-]
       % R_i         see Harris                            [m]
       % gamma_tick  see Harris (gamma' in print)          [-]
       
       % set         variable to check if the struct has been properly
       %             initialized
       
       physical;
       % Struct containing the phsical properties of the bearing, other
       % than geometry, that are relevant to the bearing dynamics. Needs to
       % be set using the method setPhysical after calling setGeometry.
       % (inputs)
       % rho_ball    mass density of the bearing balls     [kg/m^3]
       % E_I         Young's modulus of the balls          [N/m]
       % E_II        Young's modulus of the raceways       [N/m]
       % xi_I        Poisson's number of the balls         [-]
       % Xi_II       Poisson's number of the raceways      [-]
       
       % (calculated values)
       % m_ball      mass of a ball                        [kg]
       % J_ball      moment of inertia of a ball           [kg*m^2]
       
       % set         variable to check if the struct has been propperly
       %             initialized       

       % loadcases to start a simulation are arrays with the following
       % contents (in order)
       
       % F_a       static axial load                        [N]
       % F_ry      static vertical radial load              [N]
       % F_rz      static horizontal radial load            [N]
       % Theta_y   axial alignment error (rotation around y-axis) [deg]
       % Theta_z   axial alignment error (rotation around z-axis) [deg]
       % n         rotational speed                         [RPM] (non-SI!)

    end
    
    methods
        function obj = BallBearing()
            obj.geometry=struct(...
                'D',0,...
                'd_m',0,...
                'alpha_free',0,...
                'r_i',0,'r_o',0,...
                'z',0,...
                'psi',0,...
                'f_i',0,'f_o',0,...
                'gamma_tick',0,...
                'R_i',0,...
                'set',0);
            obj.physical=struct(...
                'rho_ball',0,...
                'E_I',0,'E_II',0,...
                'xi_I',0,'xi_II',0,...
                'm_ball',0,...
                'J_ball',0,...
                'set',0);
        end
        
        function res=calcStiffness(obj, loads, init, h_force)
            %for loads and init see calcDisplacement. h_force is the
            %force step width for the numerical derivation.
            [a,init]=obj.calcDisplacement(loads, 'outer', init);
            res=ones(3,3);
            if isempty(h_force)
                h_force=sqrt(eps);
            end
            for i=1:3;
                load_pre=loads;
                load_pre(i)=load_pre(i)-h_force/2;
                load_aft=loads;
                load_aft(i)=load_aft(i)+h_force/2;
                [~,~,pre_disp]=obj.calcDisplacement(load_pre, a.raceControl, init);
                [~,~,aft_disp]=obj.calcDisplacement(load_aft, a.raceControl, init);
                res(:,i)=h_force./(aft_disp(1:3)-pre_disp(1:3));
            end
            %res(abs(res)>1e17)=inf; %For better readability, stiffnesses
            %greater than the calculation precision can be set to inf
        end
        
%         Function for debugging analytical jacobians.

%         function [res,anaJak]=calcNumJakobi(obj, loadcase, rc_scenario,...
%                 init_conditions)
%             
%             if isempty(init_conditions)
%                 valVec=obj.solvInit(obj.geometry.z);
%             else
%                 valVec=init_conditions;
%             end
%             loadStruct = obj.makeLoad(loadcase);
%             
%             [base_line,anaJak] = obj.solvStitcher(obj.geometry,...
%                 obj.physical, loadStruct, rc_scenario, valVec);
%             
%             jac=zeros(length(valVec));
%             h=1e-13;
%             for i=1:length(valVec)
%                 valAft=valVec;
%                 valAft(i)=valAft(i)+h;
%                 res_aft=obj.solvStitcher(obj.geometry, obj.physical,...
%                 loadStruct, rc_scenario, valAft);
%                 jac(:,i)=(res_aft-base_line)/h;
%             end
%             res=jac;
%         end
                
                                     
        
        function setGeometry(obj, D, d_m, alpha_free, r_i, r_o, z,...
                Delta_P_D, psi_0)
            % setGeometry Sets the BallBearing object geometry used in
            % displacement and stiffness calculations
            obj.geometry.D=D;
            obj.geometry.d_m=d_m;
            obj.geometry.alpha_free=alpha_free;
            obj.geometry.r_i=r_i;
            obj.geometry.r_o=r_o;
            obj.geometry.z=z;
            obj.geometry.Delta_P_D=Delta_P_D;
            obj.geometry.psi=obj.genPsi(z, psi_0);
            obj.geometry.f_i=r_i/D;
            obj.geometry.f_o=r_o/D;
            obj.geometry.R_i=obj.calcRadLoc(d_m, obj.geometry.f_i, D,...
                alpha_free);
            obj.geometry.gamma_tick=D/d_m; %p.308
            obj.geometry.set=1;
            
            obj.physical.set=0; %obj.physical needs updating after this
        end
        
        function setPhysical(obj, rho_ball, E_I, E_II, xi_I, xi_II)
            if obj.geometry.set==0
                error(['Bearing geometry not set. Use method'...
                    ' setGeometry before calling this method']);
            else
                obj.physical.rho_ball=rho_ball;
                obj.physical.m_ball=...
                    obj.ballMass(obj.geometry.D,obj.physical.rho_ball);
                obj.physical.J_ball=obj.ballInertiaMoment...
                    (obj.geometry.D,obj.physical.m_ball);
                obj.physical.E_I=E_I;
                obj.physical.E_II=E_II;
                obj.physical.xi_I=xi_I;
                obj.physical.xi_II=xi_II;
                obj.physical.set=1;
            end
        end
        
        function phys = getPhysical(obj) %no pun intended
            %generates a cell array of the input values for setPhysical
            %so the properties can be updated by calling
            %BallBearing.setPhysiscal(Ballbearing.getPhysical{:})
            phys = cell(1,5);
            phys{1} = obj.physical.rho_ball;
            phys{2} = obj.physical.E_I;
            phys{3} = obj.physical.E_II;
            phys{4} = obj.physical.xi_I;
            phys{5} = obj.physical.xi_II;       
        end
        
        function [resStruct, solVec, dispArray]=calcDisplacement(obj,...
                                 loadcase, rc_assumption, init_conditions)
            %Main method to run a simulation and obtain a displacement.
            %Retruns a formated struct (resStruct) of the results as well
            %as raw data. Extracting the required informaiton from
            %resStruct is the most user friendly way of using this
            %function.
            %solVec can directly be used as init_conditions for a similar
            %load case
            
            %loadcase contains the particular load applied for the
            %simulation and must be an array of the following shape:
            %[Axial load, radial Y force, radial Z force, angular error
            %around the y-axis, angular error around the z-axis, rotational
            %speed in RPM]
            
            rc_scenario=rc_assumption;
            if (obj.geometry.set==0)
                error('Bearing geometry not set. Use method setGeometry.');
            end
            if (obj.physical.set==0)
                error(['Bearing physical properties not set. Use method'...
                    ' setPhysical.']);
            end

            bea_geo=obj.geometry;
            bea_phys=obj.physical;
            
            bea_load=obj.makeLoad(loadcase);
            
            if isempty(init_conditions)
                init=obj.solvInit(bea_geo.z);
            else
                init=init_conditions;
            end
            
            options = optimoptions('fsolve','Display','iter', ...
                                'MaxFunEvals', 1e9,...
                                'MaxIter', 1e9, ...
                                'Jacobian', 'off',...
                                'FunValCheck','on',...
                                'TolFun',1e-18,...
                                'TolX', 1e-18, ...
                                'Algorithm','trust-region-dogleg');
                               % 'JacobPattern',obj.jPat(bea_geo.z),... % only used by trust-region reflective
                               % levenberg-marquardt
                               % trust-region-reflective
                               % trust-region-dogleg
                               %'FinDiffRelStep', 1e-9,...
                               %  'FinDiffType', 'central',...
                               
            loadSolv=@(val)obj.solvStitcher(bea_geo, bea_phys, bea_load,...
                rc_scenario, val);
            
            [solVec,fval,exitflag,output] = ... 
                 fsolve(loadSolv,init,options);                            
            
            resStruct = obj.genSolvStruct(solVec, bea_geo);
            
            [resStruct,rc_inner] = obj.addSigmaMaxAlpha(resStruct, bea_geo,...
                bea_phys, bea_load);
            
            %checks if the race control assumption was correct. If it was
            %not, the solver is re-run. This only checks for one ball
            %(with median stress) and decides the race control scenario for
            %the entire bearing. It is a know limitation of this class,
            %that race control is not determined for every ball separately.
            
            if rc_inner
                if strcmp(rc_scenario,'outer')
                display('Inner race control detected. Re-running solver.');
                rc_scenario='inner';
                loadSolv=@(val)obj.solvStitcher(bea_geo, bea_phys, bea_load,...
                    rc_scenario, val);
                [solVec,fval,exitflag,output]=fsolve(loadSolv,init,...
                    options);
                resStruct = obj.addSigmaMaxAlpha(resStruct, bea_geo,...
                bea_phys, bea_load);
                end
            else if strcmp(rc_scenario,'inner')
                    display('Outer race control detected. Re-running solver.');
                    rc_scenario='outer';
                    loadSolv=@(val)obj.solvStitcher(bea_geo, bea_phys, bea_load,...
                        rc_scenario, val);
                    [solVec,fval,exitflag,output]=fsolve(loadSolv,init,...
                        options);
                    resStruct = obj.addSigmaMaxAlpha(resStruct, bea_geo,...
                        bea_phys, bea_load);
                end
            end
            
            resStruct.raceControl=rc_scenario;
            
            dispArray=[resStruct.delta_a, resStruct.delta_ry...
                resStruct.delta_rz resStruct.M_y...
                resStruct.M_z];
        end
  
            
        function [res, jac] = solvStitcher(obj, geo, phys, load, scenario, valVec)
            %Main function that calls the methods to tell the solver how
            %far his current guess is from a correct solution.
            
            solv = obj.genSolvStruct(valVec, geo);
            
            [rc_id, lambda_i, lambda_o]=obj.makeRaceControl(scenario);
            omega = (load.n/60)*2*pi;
            
            
            A_1 = obj.raceCenterDistAx(solv.delta_a, geo, load);
            A_2 = obj.raceCenterDistRad(solv.delta_ry, solv.delta_rz, geo);
            
            [cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o]=...
                obj.trigSincosd(A_1, A_2, solv.X_1, solv.X_2,...
                solv.delta_i, solv.delta_o, geo.D, geo.f_i, geo.f_o);
            
            %rho_ci = obj.curvInner(cosd(geo.alpha_free), geo.f_i, geo.D,...
            %    geo.d_m);
            %rho_co = obj.curvOuter(cosd(geo.alpha_free), geo.f_o, geo.D,...
            %    geo.d_m);
            
            %More precise solution taking into account the change in attack
            %angle. Very hard to implement into an analytical Jacobi
            %matrix. Error seems negligible in tests.
            rho_ci = obj.curvInner(cos_alp_i, geo.f_i, geo.D, geo.d_m);
            rho_co = obj.curvOuter(cos_alp_o, geo.f_o, geo.D, geo.d_m);
            
            [~, ~, ~, ~, K_i] = obj.calcStress(solv.delta_i, rho_ci, phys);
            [~, ~, ~, ~, K_o] = obj.calcStress(solv.delta_o, rho_co, phys);

            M_g = obj.gyroMoment(phys.J_ball, omega, sin_alp_i,...
                cos_alp_i, sin_alp_o, cos_alp_o, geo.gamma_tick, rc_id);
            
            Geo_i_err = obj.innerGeometryError(solv.X_1, solv.X_2, A_1,...
                A_2, solv.delta_i, geo.f_i, geo.D);
            Geo_o_err = obj.outerGeometryError(solv.X_1, solv.X_2,...
                solv.delta_o, geo.f_o, geo.D);
            F_b_ax_err = obj.ballForceAxial(solv.X_1, solv.X_2, A_1,...
                A_2, solv.delta_i, solv.delta_o, geo.f_i, geo.f_o,...
                lambda_o, lambda_i, geo.D, K_o, K_i, M_g);
            F_b_r_err = obj.ballForceRadial(solv.X_1, solv.X_2, A_1,...
                A_2, solv.delta_i, solv.delta_o, geo.f_i, geo.f_o,...
                lambda_o, lambda_i, geo.D, K_o, K_i, geo.d_m, sin_alp_i,...
                cos_alp_i, sin_alp_o, cos_alp_o, M_g, geo.gamma_tick,...
                omega, phys.m_ball, rc_id);
            F_ax_err = obj.forceAxial(solv.X_1, solv.X_2, A_1, A_2,...
                solv.delta_i, K_i, lambda_i, geo.D, geo.f_i, load.F_a,M_g);
            F_rz_err = obj.forceRadial(solv.X_1, solv.X_2, A_1, A_2,...
                solv.delta_i, K_i, lambda_i, geo.D, geo.f_i, load.F_rz,...
                M_g, cosd(geo.psi));
            F_ry_err = obj.forceRadial(solv.X_1, solv.X_2, A_1, A_2,...
                solv.delta_i, K_i, lambda_i, geo.D, geo.f_i, load.F_ry,...
                M_g, sind(geo.psi));            
            M_by_err = obj.bendingMoment(solv.X_1, solv.X_2, A_1, A_2,...
                solv.delta_i, M_g, geo.D, K_i, geo.f_i, lambda_i,...
                solv.M_y, cosd(geo.psi), geo.R_i);
            M_bz_err = obj.bendingMoment(solv.X_1, solv.X_2, A_1, A_2,...
                solv.delta_i, M_g, geo.D, K_i, geo.f_i, lambda_i,...
                solv.M_z, sind(geo.psi), geo.R_i);

            %jac= double(obj.makeJacobi(geo, phys, load, solv, K_i, K_o,...
            %    lambda_i, lambda_o, M_g, A_1, A_2));
            
            %I have experimented (a lot) with analytical jacobians. I found
            %the resutls underwhelming.
            
            jac=0;
            res = double([F_ax_err F_ry_err F_rz_err M_by_err M_bz_err ...
                Geo_i_err' Geo_o_err' F_b_ax_err' F_b_r_err']);
            
        end
        
        function [Q, sigma_max, a, E, K, b] = calcStress(obj, delta, rho, phys)
            %This function calculates everything about a Herzian ellipsoid
            %that the class needs in various places. In order: force,
            %peak stress, semi-major axis, elliptical integral, stiffness,
            %and semi-minor axis.
            persistent kappa
            sum_rho = obj.curvSum(rho);
            
            %precise (see way down) determines the way that the kappa value
            % of the ellipse is calculated:
            % 0: By correlation as reccomended by Harris
            % 1: Numerically
            % Numerical calculation is currently very slow and generally
            % not worth the extra calculation time
            if precise
                if isempty(kappa)
                    R_y = 1./(rho(:,1)+rho(:,3));
                    R_x = 1./(rho(:,2)+rho(:,4));
                    kappaInit = obj.ellipseCorr(R_x, R_y);
                else kappaInit=kappa;
                end
                
                kapSolv=@(kap)obj.hertzError(kap,rho);
                options = optimoptions('fsolve','Display','none', ...
                    'MaxFunEvals', 1e9,...
                    'MaxIter', 1e9, ...
                    'Jacobian', 'on',...
                    'FunValCheck','on',...
                    'Algorithm','levenberg-marquardt',...
                    'TolFun',1e-10,...
                    'TolX', 1e-10 ...
                    );
                % trust-region-reflective 'FinDiffRelStep', 1e-6,...
                % trust-region-dogleg
                
                [kappa,~,exitflag]=fsolve(kapSolv, kappaInit, options);
                if ~exitflag
                    error('Geometry Error: Hertzian ellipsoid did not converge');
                end
                if ~isreal(kappa)
                    error('Kappa not real')
                end
            else
%                 F_rho = curvDiff(rho);
                R_y = 1./(rho(:,1)+rho(:,3));
                R_x = 1./(rho(:,2)+rho(:,4));
                kappa = obj.ellipseCorr(R_x, R_y);
            end
            [F, E] = ellipke(1-1./kappa.^2); %Non-standard definition 
            %of the input parameter:
            % m{matlab}=k^2{wikipedia}=1-1/kappa^2{Jones/Harris}
            
            delta_star=2*F/pi.*(pi./(2*kappa.^2.*E)).^(1/3);
            a_star=((2*kappa.^2.*E)/pi).^(1/3);
            b_star=((2*E)./(pi*kappa)).^(1/3);
            %Gupta has 2 instead of pi in denominator in his book,
            %but pi is correct
            one_over_E_tick=obj.elasticity(phys.xi_I, phys.xi_II,...
                phys.E_I, phys.E_II);
            Q = engPower(delta./delta_star.*2./sum_rho,(3/2))...
                .*(2*sum_rho)./(3*one_over_E_tick);
            a = a_star.*engPower(3/2*Q./sum_rho*one_over_E_tick,(1/3));
            b = b_star.*engPower(3/2*Q./sum_rho*one_over_E_tick,(1/3));
            sigma_max=(3*Q)./(2*pi*a.*b);
            K=pi*kappa*2./(one_over_E_tick).*engPower(2*E./(9*sum_rho.*F.^3),0.5);
        end
        
        function [aftStruct,rc_inner] = addSigmaMaxAlpha(obj, befStruct, geo, phys,...
                load)
            aftStruct=befStruct;
            A_1 = obj.raceCenterDistAx(aftStruct.delta_a,...
                geo, load);
            A_2 = obj.raceCenterDistRad(aftStruct.delta_ry,...
                aftStruct.delta_rz, geo);
            
            [cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o] =...
                obj.trigSincosd(A_1, A_2,...
                aftStruct.X_1, aftStruct.X_2,...
                aftStruct.delta_i, aftStruct.delta_o,...
                geo.D, geo.f_o, geo.f_i);
            
            rho_ci = obj.curvInner(cos_alp_i, geo.f_i, geo.D,...
                geo.d_m);
            rho_co = obj.curvOuter(cos_alp_o, geo.f_o, geo.D,...
                geo.d_m);
            
            [Q_i, aftStruct.sigma_max_i, a_i, E_i] = obj.calcStress(...
                aftStruct.delta_i, rho_ci, phys);
            [Q_o, aftStruct.sigma_max_o, a_o, E_o] = obj.calcStress(...
                aftStruct.delta_o, rho_co, phys);
              
            aftStruct.alpha_o=atan2d(sin_alp_o,cos_alp_o);
            aftStruct.alpha_i=atan2d(sin_alp_i,cos_alp_i);
            
            aftStruct.Q_i=Q_i;
            aftStruct.Q_o=Q_o;
            rc_trig=cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_i;
            
            rc_inner=obj.rcCheck(Q_o, a_o, E_o, rc_trig, Q_i, a_i, E_i);
            
            omega=(load.n/60)*2*pi;
            rc_id=not(rc_inner);
            switch rc_id
            case 0 %inner race control
                beta=atan2d(sin_alp_i,cos_alp_i-geo.gamma_tick);
            case 1 %outer race control
                beta=atan2d(sin_alp_o,cos_alp_o+geo.gamma_tick);
            end
            aftStruct.omega_m = obj.ballOrbitSpeed(omega, sin_alp_i,...
                cos_alp_i, sin_alp_o, cos_alp_o, geo.gamma_tick, rc_id);
            aftStruct.omega_R = obj.ballSpinSpeed(sin_alp_i, cos_alp_i,...
                sin_alp_o, cos_alp_o, omega, geo.gamma_tick, beta);
            aftStruct.beta = beta;
            
            %This is good for debugging:
            %clipboard('copy', [cos_alp_o(1) cos_alp_i(1);sin_alp_o(1) sin_alp_i(1); a_o(1) a_i(1)]);
        end
        
        function F_c = centrifugalForce(obj, d_m, omega, sin_alp_i,...
                cos_alp_i, sin_alp_o, cos_alp_o, gamma_tick, m_ball, rc_id)
            omega_m=obj.ballOrbitSpeed(omega, sin_alp_i, cos_alp_i,...
                sin_alp_o, cos_alp_o, gamma_tick, rc_id);
            F_c=0.5*m_ball*d_m*omega_m.^2;
        end
        function [cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o]=...
              trigSincosd(obj, A_1, A_2, X_1, X_2, delta_i, delta_o, D,...
              f_i, f_o)
            cos_alp_i = obj.cosICAngle(A_2, X_2, f_i, D, delta_i);
            cos_alp_o = obj.cosOCAngle(X_2, f_o, D, delta_o);
            sin_alp_i = obj.sinICAngle(A_1, X_1, f_i, D, delta_i);
            sin_alp_o = obj.sinOCAngle(X_1, f_o, D, delta_o);
        end
        
        
        function M_g = gyroMoment(obj, J, omega, sin_alp_i, cos_alp_i,...
                sin_alp_o, cos_alp_o, gamma_tick, rc_id)
            switch rc_id
                case 0 %inner race control
                    beta=atan2d(sin_alp_i,cos_alp_i-gamma_tick);
                case 1 %outer race control
                    beta=atan2d(sin_alp_o,cos_alp_o+gamma_tick);
            end
            omega_m =  obj.ballOrbitSpeed(omega, sin_alp_i, cos_alp_i,...
                sin_alp_o, cos_alp_o, gamma_tick, rc_id);
            omega_R = obj.ballSpinSpeed(sin_alp_i, cos_alp_i, sin_alp_o,...
                cos_alp_o, omega, gamma_tick, beta);
            M_g=J*omega_R.*omega_m.*sind(beta);
        end
        
        function res = ballForceRadial(obj, X_1, X_2, A_1, A_2, delta_i,...
                delta_o, f_i, f_o, lambda_o, lambda_i, D, ...
                K_o, K_i, d_m, sin_alp_i, cos_alp_i,...
                sin_alp_o, cos_alp_o, M_g, gamma_tick, omega,...
                m_ball, rc_id)            
            F_c = obj.centrifugalForce(d_m, omega, sin_alp_i, cos_alp_i,...
                sin_alp_o, cos_alp_o, gamma_tick, m_ball, rc_id);
            res = (...
                    (K_o.*engPower(delta_o,1.5).*X_2)...
                    + ((lambda_o*M_g.*X_1)/D)...
                  )...
                  ./((f_o-0.5)*D+delta_o)...
                  -...
                  ( ...
                     (K_i.*engPower(delta_i,1.5).*(A_2-X_2))...
                     +((lambda_i*M_g)/D.*(A_1-X_1)) ...
                  )...
                    ./((f_i-0.5)*D+delta_i)...
                  -F_c;
        end
        
    end
 
    
    methods (Static = true)
        
        function res = rcCheck(Q_o, a_o, E_o, trig, Q_i, a_i, E_i)
            Q_o_sort=sort(Q_o);
            rep_value=Q_o_sort(round(length(Q_o)/2));
            rep_element=find(Q_o==rep_value);
            rep_element=rep_element(1);
            Q_o_med=Q_o(rep_element);
            Q_i_med=Q_i(rep_element);
            a_o_med=a_o(rep_element);
            a_i_med=a_i(rep_element);
            E_i_med=E_i(rep_element);
            E_o_med=E_o(rep_element);
            trig_med=trig(rep_element);
            res=(Q_o_med*a_o_med*E_o_med*trig_med<=Q_i_med*a_i_med*E_i_med);
        end
        
        function res = solvInit(z)
            delta_X_1 = ones(z,1)*1e-5;
            delta_X_2 = ones(z,1)*1e-5;
            delta_i = ones(z,1)*1e-5;
            delta_o = ones(z,1)*1e-5;
            delta_a = 1e-5;
            delta_ry = 0;
            delta_rz = 0;
            M_y = 0;
            M_z = 0;
            res = [delta_X_1' delta_X_2' delta_i' delta_o' delta_a...
                delta_ry delta_rz M_y M_z];
        end
        
        function loadStruct = makeLoad(loadArray)
            loadStruct=struct(...
                'F_a',loadArray(1), ...
                'F_ry',loadArray(2), ...
                'F_rz',loadArray(3), ...
                'Theta_y',loadArray(4), ...
                'Theta_z',loadArray(5), ...
                'n', loadArray(6)...
                );
        end

        function solvStruct = genSolvStruct(valueVec, geo)
            valMat = vec2mat(valueVec, geo.z);
            solvStruct = struct(...
                'X_1', valMat(1,:)'+sind(geo.alpha_free)*(geo.f_o-0.5)*...
                    geo.D*ones(geo.z,1),...
                'X_2', valMat(2,:)'+cosd(geo.alpha_free)*(geo.f_o-0.5)*...
                    geo.D*ones(geo.z,1),...
                'delta_i', valMat(3,:)',...
                'delta_o', valMat(4,:)',...
                'delta_a', valueVec(end-4),...
                'delta_ry', valueVec(end-3),...
                'delta_rz', valueVec(end-2),...
                'M_y', valueVec(end-1)*1e4,...
                'M_z', valueVec(end)*1e4);
        end        
        
        function [rc_id, lambda_i, lambda_o]=makeRaceControl(scenario)
            switch scenario
                case 'inner'
                    rc_id=0;
                    lambda_i=1;
                    lambda_o=1;
                case 'outer'
                    rc_id=1;
                    lambda_i=0;
                    lambda_o=2;
                otherwise
                    error(['Invalid Race control mode.'...
                        ' Valid modes: inner, outer'])
            end
        end
        
        function res = jPat(z)
        %Renders a Jacobi matrix sparsity pattern for the problem not
        %currently in use, because the only Matlab solver that uses this
        %input is ill suited to this problem.
            ball=eye(z);
            whole=ones(z,5);
            whole_force=ones(5,4*z+5);
            res = [whole_force;...
                  ball ball ball ball whole; ...
                  ball ball ball ball whole; ...
                  ball ball ball ball whole; ...
                  ball ball ball ball whole];
        end
        
        function [err, jac] = hertzError(kappa, rho)
            %Funciton to nummerically solve the herzian contact ellipsoid.
            %This is a lot slower than using the correlation, which is also
            %included in this class.
            F_rho_dir = ((rho(:,1)-rho(:,2))+(rho(:,3)-rho(:,4)))...
                    ./sum(rho,2);

            [F, E] = ellipke(1-1./kappa.^2); 
            %Non-Standard definition of the input parameter:
            %   m{matlab}=k^2{wikipedia}=1-1/kappa^2{Harris}            
            F_rho_indir=((kappa.^2+1).*E-2*F)./((kappa.^2-1).*E);
            err=F_rho_dir-F_rho_indir;
            jac=diag((2*(F.^2 + 3*kappa.^2.*E.^2 - 2*E.*F -...
                2*kappa.^2.*E.*F))./(kappa.*E.^2.*(kappa.^2 - 1).^2));
%           Jacobi matrix generated by the following symmath:
%             clear;
%             syms kappa
%             [F, E] = ellipke(1-1./kappa.^2);
%             F_rho_indir=((kappa.^2+1).*E-2*F)./((kappa.^2-1).*E);
%             simplify(diff(-F_rho_indir, kappa))
        end
        
        function rho_ci = curvInner(cos_alp_i, f_i, D, d_m)
            %calculates the curvature at the inner raceway contact based on
            %Jones' derivation.
            gamma = D*cos_alp_i/d_m;
            %I: Ball
            %II: Raceway
            %1: Ellipsoid axis in raceway curvature direction (typically
            %minor axis)
            %2: Ellipsoid axis in grove curvature direction (typically
            %major axis)
            rho_I1 = 2/D+0*gamma;
            rho_I2 = 2/D+0*gamma;
            rho_II1 = (2/D)*(gamma./(1-gamma));
            rho_II2 = -1/(f_i*D)+0*gamma;
            rho_ci=[rho_I1 rho_I2 rho_II1 rho_II2];
        end
            
        function rho_co = curvOuter(cos_alp_o, f_o, D, d_m)
            %calculates the curvature at the outer raceway contact based on
            %Jones' derivation.
            gamma = D.*cos_alp_o/d_m;
            %I: Ball
            %II: Raceway
            %1: Ellipsoid axis in raceway curvature direction (typically
            %minor axis)
            %2: Ellipsoid axis in grove curvature direction (typically
            %major axis)
            rho_I1 = 2/D+0*gamma;
            rho_I2 = 2/D+0*gamma;
            rho_II1 = -(2/D)*(gamma./(1+gamma));
            rho_II2 = -1/(f_o*D)+0*gamma;
            rho_co=[rho_I1 rho_I2 rho_II1 rho_II2];
        end
        
        function one_over_E_tick = elasticity(xi_I, xi_II,E_I, E_II)
            one_over_E_tick = (1-xi_I^2)/E_I+(1-xi_II^2)/E_II;
        end
        
        function sum_rho = curvSum(rho)
            %p.63 Harris vol.2 5th edition
            sum_rho = sum(rho,2);
        end
        
        function F_rho = curvDiff(rho)
            F_rho = ((rho(:,1)-rho(:,2))+(rho(:,3)-rho(:,4)))./sum(rho,2);
        end
        
        function m = ballMass(D, rho)
            m = rho*(4/3)*pi*(D/2)^3;
        end
        
        function J = ballInertiaMoment(D, m)
            J = 2/5*m*(D/2)^2;
        end
        
        function kappa = ellipseCorr(R_x, R_y)
           kappa = 1.0339*(R_x./R_y).^0.636;
        end
        
        function psi=genPsi(z,psi_0)
            %devide a 360 deg circle into z equal parts starting at psi_0
            psi=((0:(360/z):360-(360/z))+psi_0)';
        end;   
         
        
        function res = forceAxial(X_1, X_2, A_1, A_2, delta_i,...
                K_i, lambda_i, D, f_i, F_a, M_g)
            if precise %switch to activate rounding after addition.
                %May help with some numerical problems.
                res = F_a - ballsum(vpa(...
                (...
                    K_i.*(A_1-X_1).*engPower(delta_i,1.5) ...
                    -(lambda_i*M_g)/D.*(A_2-X_2)...
                )...
                ./((f_i-0.5)*D+delta_i) ...
                ));
            else
                res = F_a - ballsum(...
                (...
                    K_i.*(A_1-X_1).*engPower(delta_i,1.5) ...
                    -(lambda_i*M_g)/D.*(A_2-X_2)...
                )...
                ./((f_i-0.5)*D+delta_i) ...
                );
            end
        end
        
        function res = forceRadial(X_1, X_2, A_1, A_2, delta_i,...
                K_i, lambda_i, D, f_i, F_r, M_g, trig_psi)
            if precise
            res = F_r - ballsum(vpa((...
                (...
                    K_i.*(A_2-X_2).*engPower(delta_i,1.5) ...
                    +(lambda_i*M_g)/D.*(A_1-X_1)...
                )...
                ./((f_i-0.5)*D+delta_i) ...
                ).*trig_psi));
            else
                res = F_r - ballsum((...
                    (...
                        K_i.*(A_2-X_2).*engPower(delta_i,1.5) ...
                        +(lambda_i*M_g)/D.*(A_1-X_1)...
                    )...
                    ./((f_i-0.5)*D+delta_i) ...
                    ).*trig_psi);
            end
        end
        
        function res = bendingMoment(X_1, X_2, A_1, A_2, delta_i, M_g,...
                D, K_i, f_i, lambda_i, M, trig_psi, R_i)
            if precise           
            res= M - ballsum(vpa(...
                (...
                    (...
                        (...
                            (...
                                (K_i.*(A_1-X_1).*engPower(delta_i,1.5))...
                                -...
                                ((lambda_i*M_g)/D.*(A_2-X_2))...
                            )*R_i...
                        )...
                        ./...
                        (...
                            (f_i-0.5)*D+delta_i...
                        )...
                    )...
                + lambda_i*f_i*M_g...
                )...
                .*trig_psi));
            else
            res= M - ballsum(...
                (...
                    (...
                        (...
                            (...
                                (K_i.*(A_1-X_1).*engPower(delta_i,1.5))...
                                -...
                                ((lambda_i*M_g)/D.*(A_2-X_2))...
                            )*R_i...
                        )...
                        ./...
                        (...
                            (f_i-0.5)*D+delta_i...
                        )...
                    )...
                + lambda_i*f_i*M_g...
                )...
                .*trig_psi);   
            end
        end
            
        
        function res = innerGeometryError(X_1, X_2, A_1,...
                A_2, delta_i, f_i, D)
            res = (A_1-X_1).^2+(A_2-X_2).^2-((f_i-0.5)*D+delta_i).^2;
        end

        function res = outerGeometryError(X_1, X_2, delta_o, f_o, D)
            res = X_1.^2+X_2.^2-((f_o-0.5)*D+delta_o).^2;
        end
        
        function res = ballForceAxial(X_1, X_2, A_1, A_2, delta_i,...
                delta_o, f_i, f_o, lambda_o, lambda_i, D, K_o, K_i, M_g)
            res = (((lambda_o*M_g.*X_2)/D) - (K_o.*engPower(delta_o,1.5).*X_1))...
                    ./((f_o-0.5)*D+delta_o)...
                 + ( ...
                    (K_i.*engPower(delta_i,1.5).*(A_1-X_1))...
                    -((lambda_i*M_g)/D.*(A_2-X_2)) ...
                 )...
                    ./((f_i-0.5)*D+delta_i);
        end

        function  sin_alp_i = sinICAngle(A_1, X_1, f_i, D, delta_i)
            sin_alp_i=(A_1 - X_1)./((f_i - 0.5)*D + delta_i);
        end
        function  sin_alp_o = sinOCAngle(X_1, f_o, D, delta_o)
            sin_alp_o = X_1./((f_o - 0.5)*D + delta_o);
        end
        function  cos_alp_i = cosICAngle(A_2, X_2, f_i, D, delta_i)
            cos_alp_i = (A_2 - X_2)./((f_i - 0.5)*D + delta_i);
        end
        function  cos_alp_o = cosOCAngle(X_2, f_o, D, delta_o)
            cos_alp_o = X_2./((f_o - 0.5)*D + delta_o);
        end
       
        function R_i = calcRadLoc(d_m, f_i, D ,alpha_free)
            R_i = 0.5*d_m+(f_i-0.5)*D*cosd(alpha_free);
        end
        
        function A_1 = raceCenterDistAx(delta_a, geo, load)
            A_1 = (geo.f_i+geo.f_o-1)*geo.D*sind(geo.alpha_free)...
                + delta_a + ...
                load.Theta_y*geo.R_i*cosd(geo.psi)+...
                load.Theta_z*geo.R_i*sind(geo.psi)+...
                geo.Delta_P_D;
            if ~isreal(A_1)
                disp('A_1 complex')
            end
        end
        
        function A_2 = raceCenterDistRad(delta_ry, delta_rz, geo)
            %Jones includes a radial clearance here, but Harris doesn't.
            %Can be added without a problem.
            A_2 = (geo.f_i+geo.f_o-1)*geo.D*cosd(geo.alpha_free) + ...
                delta_rz*cosd(geo.psi)+delta_ry*sind(geo.psi);
            if ~isreal(A_2)
                disp('A_2 complex')
            end
        end

    function omega_m = ballOrbitSpeed(omega, sin_alp_i, cos_alp_i,...
        sin_alp_o, cos_alp_o, gamma_tick, rc_id)
        %Calculates the orbital speed of the balls.
        switch rc_id
            case 0 %inner race control
                omega_m=omega*(((cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_o)...
                    -gamma_tick.*cos_alp_o)./...
                    (1+(cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_o)));
            case 1 %outer race control
                omega_m=omega*((1-gamma_tick.*cos_alp_i)./...
                    (1+(cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_o)));
        end
    %trigonometric identity: cosd(alpha_i-alpha_o) =
    %                  cosd(alpha_i)*cosd(alpha_o)+sind(alpha_i)*sind(alpha_o)
    end

    function omega_R = ballSpinSpeed(sin_alp_i, cos_alp_i, sin_alp_o,...
            cos_alp_o, omega, gamma_tick, beta)
        %Calculates the spinning speed of the balls.
        omega_R=-omega./...
            (...
                (...
                    (cos_alp_o+tand(beta).*sin_alp_o)...
                    ./(1+gamma_tick*cos_alp_o)...
                    +...
                    (cos_alp_i+tand(beta).*sin_alp_i)...
                    ./(1-gamma_tick*cos_alp_i)...
                )...
                *gamma_tick.*cosd(beta)...
            );
    end

    end    
end

function res = engPower(b, e)
    %This is a hack that vastly improves the convergence of the model:
    %Since the force model is F~delta^(3/2), if the solver ever guesses a
    %negative delta, it will crash. With this hack, it will be led back to
    %a correct solution. CAUTION, however: This class cannot desctibe cases
    %where one or more balls are truely out of contact. Those cases will be
    %indicated by a negative delta, but the results produced here are 
    %useless.
    res = sign(b).*(abs(b).^e);
end


function res = ballsum(a)
    %this helps a bit with precision, without costing too much in computing
    %res=sum(a);
    res=sort(sum(a(a>0)))+sort(sum(a(a<0)),'descend');
end

function res=precise()
    res=0;
    %switch that chooses between the precise and approximate solutions
    % in some functions. 0 means aproximate, 1 means precise. Precise mode
    % takes orders of magnitude longer to compute, and generally does not
    % change the solution siginificantly
end

% function res = sind(x)
%     res=sin(vpa(pi)*x/180);
% end
% 
% function res = cosd(x)
%     res=cos(vpa(pi)*x/180);
% end
% 
% 
% function res = tand(x)
%     res=tan(vpa(pi)*x/180);
% end
% 
% function res = atan2d(x,y)
%     res=atan2(vpa(pi)*x/180,vpa(pi)*y/180);
% end
