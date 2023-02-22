"""
Function generates 3rd order approximation to periodic solution about collinear
libration point in the CR3BP
"""


struct Struct_out_halo_analytical
   x0
   period::Float64
   fullstates
end


"""
   halo_analytical_construct(mu::Float64, lp::Int, Az_km::Float64, lstar::Float64, northsouth::Int, phi::Real=0.0)

Generate analytical approx. of halo about collinear libration point

# Arguments
      mu
      lp: (1, 2, or 3)
      Az_km
      northsouth (int): (1 or 3)
      phi: 0.0 or pi

# Returns
      T (float): 3rd order period of halo
      X0
      x_analytic_synodic
      y_analytic_synodic
      z_analytic_synodic
      xdot_analytic_synodic
      ydot_analytic_synodic
      zdot_analytic_synodic
"""
function halo_analytical_construct(
   mu::Float64, lp::Int, Az_km::Real, lstar::Real, northsouth::Int, phi::Real=0.0
)
   # function to define lagrange points
   LP = lagrangePoints(mu)

   # choose Lagrange point of interest
   if lp == 1
      gammaL = abs((1-mu) - LP[lp,1]);
   elseif lp == 2
      gammaL = abs((1-mu) - LP[lp,1]);
   elseif lp ==3
      gammaL = abs(mu - LP[lp,1]);
   end

   # computation of Legendre function coefficients
   function _legendre_cn(n::Int,mu::Float64,gammaL::Float64,Lpoint::Int)
      """Function computes Legendre function nth coefficient"""
      if Lpoint == 1
          cn = (1/gammaL^3) * (mu + (-1)^n * (1-mu)*gammaL^(n+1) / (1 - gammaL)^(n+1));
      elseif Lpoint == 2
          cn = (1/gammaL^3) * ((-1)^n*mu + (-1)^n * (1-mu)*gammaL^(n+1) / (1 + gammaL)^(n+1));
      elseif Lpoint == 3
          cn = (1/gammaL^3) * (1 -mu + mu*gammaL^(n+1) / (1 + gammaL)^(n+1));
      end
      return cn
   end

   # compute coefficients c2, c3, c4
   c2 = _legendre_cn(2,mu,gammaL,lp);
   c3 = _legendre_cn(3,mu,gammaL,lp);
   c4 = _legendre_cn(4,mu,gammaL,lp);

   # function computes characteristic root for in-plane linearized EoM
   function _compute_cr3bp_lambda(c2)
      """Function computes roots of characteristic coupled
      x- and y-EoM"""
      # find root
      f(lmb) = lmb^4 + (c2-2)*lmb^2 - (c2 - 1)*(1 + 2*c2);
      lambda = find_zero(f, 3);
      return lambda
   end

   # compute lambda
   lambda = _compute_cr3bp_lambda(c2);
   # compute Delta (frequency coefficient from 1st order to 3rd order)
   Delta = lambda^2 - c2;
   # compute k (xy-amplitude factor)
   k = 2*lambda/(lambda^2 + 1 - c2);

   # function for coefficients necessary in third-order solution
   function _richardsonCoeff(c2,c3,c4,lambda,k)
      """Function computes coefficients for Richardson (1979) solution"""
      # coefficient d1, d2
      d1 = ( 3*(lambda^2)/k ) * (k*(6*lambda^2 - 1) - 2*lambda);
      d2 = ( 8*(lambda^2)/k ) * (k*(11*lambda^2 - 1) - 2*lambda);
      # coefficient d21
      d21 = -c3/(2*lambda^2);
      # coefficient b21, b22
      b21 = ( -3*c3*lambda/(2*d1) )*(3*k*lambda - 4);
      b22 = 3*c3*lambda/d1;
      # coefficient a21 ~ a24
      a21 = 3*c3*(k^2 - 2)/(4*(1+2*c2));
      a22 = 3*c3/(4*(1+2*c2));
      a23 = (-3*c3*lambda/(4*k*d1)) * (3*k^3*lambda - 6*k*(k - lambda) + 4);
      a24 = (-3*c3*lambda/(4*k*d1)) * (2 + 3*k*lambda);
      # b31 ~ b32
      b31 = (3/(8*d2))*(8*lambda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2)) +
          (9*lambda^2 + 1 + 2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2)));
      b32 = (1/d2)*(9*lambda*(c3*(k*b22+d21-2*a24)-c4) + (3/8)*(9*lambda^2
              + 1 + 2*c2)*(4*c3*(k*a24 - b22) + k*c4));
      # a31 ~ a32
      a31 = (-9*lambda/(4*d2)) * (4*c3*(k*a23 - b21) + k*c4*(4+k^2)) +
              ((9*lambda^2 + 1 - c2)/(2*d2))*(3*c3*(2*a23 - k*b21) + c4*(2+3*k^2));
      a32 = (-1/d2)*((9*lambda/4)*(4*c3*(k*a24 - b22) + k*c4) +
              (3/2)*(9*lambda^2 + 1 - c2)*(c3*(k*b22 + d21 - 2*a24) - c4));
      # coefficient d31 ~ d32
      d31 = (3/(64*lambda^2))*(4*c3*a24 + c4);
      d32 = (3/(64*lambda^2))*(4*c3*(a23 - d21) + c4*(4+k^2));
      # frequency correction s1, s2
      s1 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*((3/2)*c3*(2*a21*(k^2-2) -
              a23*(k^2+2)-2*k*b21) - (3/8)*c4*(3*k^4-8*k^2+8));
      s2 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*((3/2)*c3*(2*a22*(k^2-2) +
              a24*(k^2+2)+2*k*b22+5*d21) + (3/8)*c4*(12-k^2));
      # amplitude constraint param. a1, a2, l1, l2
      a1 = -(3/2)*c3*(2*a21 + a23 + 5*d21) - (3/8)*c4*(12 - k^2);
      a2 = (3/2)*c3*(a24 - 2*a22) + (9/8)*c4;
      l1 = a1 + 2*lambda^2 * s1;
      l2 = a2 + 2*lambda^2 * s2;
      return a21,a22,a23,a24,a31,a32,b21,b22,b31,b32,d21,d31,d32,d1,d2,a1,a2,s1,s2,l1,l2
   end

   # compute coefficients from Richardson 1979
   a21,a22,a23,a24,a31,a32,b21,b22,b31,b32,d21,d31,d32,d1,d2,a1,a2,s1,s2,l1,l2 = _richardsonCoeff(c2,c3,c4,lambda,k);

   # amplitudes
   Az = Az_km / (lstar*gammaL);
   Ax = sqrt((- Delta - l2*Az^2)/l1);
   Ay = k*Ax;

   # period of halo orbit
   omega1 = 0;
   omega2 = s1*(Ax)^2 + s2*(Az)^2;
   omega = 1 + omega1 + omega2;
   T = 2*pi/(lambda*omega);

   # non-dimensional time array
   time_tau = omega*range(0,stop=T,length=500);
   # switching function
   deltan = 2 - northsouth;

   # construct third-order solution
   # initialize
   x_analytic_Lframe = zeros(length(time_tau))
   y_analytic_Lframe = zeros(length(time_tau))
   z_analytic_Lframe = zeros(length(time_tau))
   xdot_analytic_Lframe = zeros(length(time_tau))
   ydot_analytic_Lframe = zeros(length(time_tau))
   zdot_analytic_Lframe = zeros(length(time_tau))

   # third order solution in L-frame (centered around libration point)
   for i in 1:length(time_tau)
      # positions in Lframe
      x_analytic_Lframe[i] = a21 * Ax^2 + a22 * Az^2 - Ax*cos(lambda*time_tau[i] + phi) + (a23*Ax^2 - a24*Az^2)*cos(2*(lambda*time_tau[i] + phi)) + (a31*Ax^3 - a32*Ax*Az^2)*cos(3*(lambda*time_tau[i] + phi));
      y_analytic_Lframe[i] = k*Ax*sin(lambda*time_tau[i] + phi) + (b21*Ax^2 - b22*Az^2)*sin(2*(lambda*time_tau[i] + phi)) + (b31*Ax^3 - b32*Ax*Az^2) * sin(3*(lambda*time_tau[i] + phi));
      z_analytic_Lframe[i] = deltan * Az*cos(lambda*time_tau[i] + phi) + deltan*d21*Ax*Az*(cos(2*(lambda*time_tau[i] + phi)) - 3) + deltan*(d32*Az*Ax^2 - d31*Az^3)*cos(3*(lambda*time_tau[i] + phi));

      # velocities
      xdot_analytic_Lframe[i] = ( lambda*Ax*sin(lambda*time_tau[i] + phi) + (a23*Ax^2 - a24*Az^2)*2*lambda*sin(2*(lambda*time_tau[i] + phi)) - (a31*Ax^3 - a32*Ax*Az^2)*3*lambda*sin(3*(lambda*time_tau[i] + phi)) );
      ydot_analytic_Lframe[i] = k*Ax*lambda*omega*cos(lambda*time_tau[i] + phi) + (b21*Ax^2 - b22*Az^2)*2*lambda*omega*cos(2*(lambda*time_tau[i] + phi)) + (b31*Ax^3 - b32*Ax*Az^2) * 3*lambda*omega*cos(3*(lambda*time_tau[i] + phi));
      zdot_analytic_Lframe[i] = ( - deltan * Az*lambda*sin(lambda*time_tau[i] + phi) - deltan*d21*Ax*Az*2*lambda*sin(2*(lambda*time_tau[i] + phi)) - deltan*(d32*Az*Ax^2 - d31*Az^3)*3*lambda*sin(3*(lambda*time_tau[i] + phi)) );
   end

   # shift solution to synodic frame (centered around barycenter of CR3BP primary masses)
   x_analytic_synodic = gammaL * x_analytic_Lframe;
   for i in 1:length(x_analytic_synodic)
      x_analytic_synodic[i] = x_analytic_synodic[i] + LP[lp,1];
   end
   y_analytic_synodic = gammaL * y_analytic_Lframe;
   z_analytic_synodic = gammaL * z_analytic_Lframe;
   xdot_analytic_synodic = gammaL * xdot_analytic_Lframe;
   ydot_analytic_synodic = gammaL * ydot_analytic_Lframe;
   zdot_analytic_synodic = gammaL * zdot_analytic_Lframe;

   # prepare state-vector at t=0 as initial guess for differential corrector
   X0 = [x_analytic_synodic[1], y_analytic_synodic[1], z_analytic_synodic[1],
         xdot_analytic_synodic[1], ydot_analytic_synodic[1], zdot_analytic_synodic[1]]

   # full states
   fullstates = hcat(x_analytic_synodic, y_analytic_synodic, z_analytic_synodic,
                     xdot_analytic_synodic, ydot_analytic_synodic, zdot_analytic_synodic, time_tau)

   # prepare function output
   return Struct_out_halo_analytical(X0, T, fullstates)
   #T, X0, x_analytic_synodic, y_analytic_synodic, z_analytic_synodic,
   #   xdot_analytic_synodic, ydot_analytic_synodic, zdot_analytic_synodic
end
